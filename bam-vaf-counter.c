/*
 * bam-vaf-counter: Count variant allele frequencies directly from BAM positions
 * 
 * This program reads aligned sequencing data from BAM files and calculates
 * VAF by directly examining bases at SNP positions, without k-mer extraction.
 * 
 * OPTIMIZATION: This version directly counts ref/alt bases at SNP positions
 * instead of extracting and looking up k-mers. This is much faster because:
 * 1. No k-mer extraction needed (major speed improvement)
 * 2. Direct position-based counting using BAM alignment coordinates
 * 3. Indexed BAM access to fetch only reads overlapping SNP positions
 * 4. Multithreaded processing with kt_pipeline
 * Falls back to sequential reading if BAM index (.bai) is not available.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <zlib.h>
#include "ketopt.h"
#include "kthread.h"

#include "htslib/htslib/sam.h"
#include "htslib/htslib/hts.h"

#include "khashl.h"
// Hash map for SNP position lookup: key is (chr_tid << 32 | position), value is pattern index
KHASHL_MAP_INIT(, snp_map_t, snp_map, uint64_t, uint32_t, kh_hash_uint64, kh_eq_generic)

typedef struct {
	char chr[256];
	int tid;  // BAM chromosome ID
	int pos;  // 0-based position
	char rsid[256];
	char ref;
	char alt;
	uint32_t ref_count;
	uint32_t alt_count;
} pattern_t;

typedef struct {
	int n, m;
	pattern_t *a;
} pattern_db_t;

// Region for targeted BAM reading
typedef struct {
	char chr[256];
	int start;
	int end;
} region_t;

typedef struct {
	int n, m;
	region_t *a;
} region_list_t;

// Load pattern file - simplified for position-based counting
pattern_db_t *load_patterns(const char *fn)
{
	FILE *fp;
	pattern_db_t *db;
	pattern_t pat;
	char ref_kmer[128], alt_kmer[128];
	int start, end;
	
	fp = fopen(fn, "r");
	if (!fp) return 0;
	
	db = (pattern_db_t*)calloc(1, sizeof(pattern_db_t));
	if (!db) {
		fclose(fp);
		return 0;
	}
	
	// Pattern file format: chr start end rsid ref alt ref_kmer alt_kmer
	// BED format: start is 0-based, end is 1-based (half-open interval [start, end))
	// The SNP position is at the 0-based 'start' coordinate
	while (fscanf(fp, "%255s%d%d%255s %c %c%127s%127s", 
	              pat.chr, &start, &end, pat.rsid, 
	              &pat.ref, &pat.alt, ref_kmer, alt_kmer) == 8) {
		if (db->n == db->m) {
			pattern_t *tmp;
			db->m = db->m ? db->m << 1 : 16;
			tmp = (pattern_t*)realloc(db->a, db->m * sizeof(pattern_t));
			if (!tmp) {
				fclose(fp);
				return db;
			}
			db->a = tmp;
		}
		pat.pos = start;  // Use 0-based position (start is already 0-based from BED)
		pat.tid = -1;  // Will be set when we have the BAM header
		pat.ref_count = 0;
		pat.alt_count = 0;
		db->a[db->n++] = pat;
	}
	
	fclose(fp);
	return db;
}

void pattern_db_destroy(pattern_db_t *db)
{
	if (db) {
		free(db->a);
		free(db);
	}
}

// Comparison function for qsort
static int region_cmp(const void *a, const void *b)
{
	const region_t *ra = (const region_t*)a;
	const region_t *rb = (const region_t*)b;
	int cmp = strcmp(ra->chr, rb->chr);
	if (cmp != 0) return cmp;
	return ra->start - rb->start;
}

// Build region list from pattern database
// Merges overlapping regions to minimize redundant reads
region_list_t *build_regions(pattern_db_t *db)
{
	region_list_t *regions;
	int i;
	
	if (!db || db->n <= 0) return 0;
	
	regions = (region_list_t*)calloc(1, sizeof(region_list_t));
	if (!regions) return 0;
	
	// Allocate space for regions
	regions->m = db->n;
	regions->a = (region_t*)malloc(regions->m * sizeof(region_t));
	if (!regions->a) {
		free(regions);
		return 0;
	}
	
	// Create regions for each SNP position
	// Since we're checking positions directly, we don't need k-mer flanking
	for (i = 0; i < db->n; ++i) {
		strncpy(regions->a[i].chr, db->a[i].chr, sizeof(regions->a[i].chr) - 1);
		regions->a[i].chr[sizeof(regions->a[i].chr) - 1] = '\0';
		regions->a[i].start = db->a[i].pos;
		regions->a[i].end = db->a[i].pos + 1;  // Single position
	}
	regions->n = db->n;
	
	// Sort regions by chromosome and start position using qsort
	qsort(regions->a, regions->n, sizeof(region_t), region_cmp);
	
	// Merge overlapping or adjacent regions on same chromosome
	int write_idx = 0;
	for (i = 1; i < regions->n; ++i) {
		if (strcmp(regions->a[write_idx].chr, regions->a[i].chr) == 0 &&
		    regions->a[write_idx].end >= regions->a[i].start) {
			// Merge overlapping or adjacent regions
			if (regions->a[i].end > regions->a[write_idx].end) {
				regions->a[write_idx].end = regions->a[i].end;
			}
		} else {
			// Non-overlapping, move to next position
			++write_idx;
			if (write_idx != i) {
				regions->a[write_idx] = regions->a[i];
			}
		}
	}
	regions->n = write_idx + 1;
	
	return regions;
}

void region_list_destroy(region_list_t *regions)
{
	if (regions) {
		free(regions->a);
		free(regions);
	}
}

// Create hash table mapping SNP position to pattern index
// Key: (tid << 32) | position, Value: pattern index
snp_map_t *create_snp_map(pattern_db_t *db, sam_hdr_t *hdr)
{
	snp_map_t *h;
	int i, absent;
	khint_t itr;
	
	h = snp_map_init();
	
	for (i = 0; i < db->n; ++i) {
		// Get BAM tid for chromosome
		int tid = sam_hdr_name2tid(hdr, db->a[i].chr);
		if (tid < 0) {
			fprintf(stderr, "Warning: chromosome %s not found in BAM header\n", db->a[i].chr);
			continue;
		}
		db->a[i].tid = tid;
		
		// Create key: (tid << 32) | position
		uint64_t key = ((uint64_t)tid << 32) | (uint32_t)db->a[i].pos;
		itr = snp_map_put(h, key, &absent);
		if (absent) {
			kh_val(h, itr) = i; // store pattern index
		} else {
			fprintf(stderr, "Warning: duplicate SNP at %s:%d\n", db->a[i].chr, db->a[i].pos);
		}
	}
	
	return h;
}

// Global data for pipeline threading
typedef struct {
	int n_thread;
	samFile *fp;
	sam_hdr_t *hdr;
	hts_idx_t *idx;
	hts_itr_t *itr;
	snp_map_t *snp_map;
	pattern_db_t *db;
	region_list_t *regions;
	int curr_region;
} pipeline_t;

// Data for each pipeline step
typedef struct {
	pipeline_t *p;
	int n, m;
	bam1_t **bam_recs;
} step_t;

// Check base at position in a BAM read and count ref/alt
static inline void count_base_at_position(bam1_t *b, int ref_pos, char ref_base, char alt_base, 
                                          uint32_t *ref_count, uint32_t *alt_count)
{
	uint32_t *cigar = bam_get_cigar(b);
	uint8_t *seq = bam_get_seq(b);
	int read_pos = 0;
	int current_ref_pos = b->core.pos;  // 0-based
	
	// Walk through CIGAR to find the read position corresponding to ref_pos
	for (uint32_t i = 0; i < b->core.n_cigar; ++i) {
		int op = bam_cigar_op(cigar[i]);
		int len = bam_cigar_oplen(cigar[i]);
		
		if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
			// Match/mismatch: both read and ref advance
			if (ref_pos >= current_ref_pos && ref_pos < current_ref_pos + len) {
				// Position is in this block
				int offset = ref_pos - current_ref_pos;
				read_pos += offset;
				
				// Get base from read
				char base = seq_nt16_str[bam_seqi(seq, read_pos)];
				
				// Count ref or alt
				if (base == ref_base) {
					__sync_fetch_and_add(ref_count, 1);
				} else if (base == alt_base) {
					__sync_fetch_and_add(alt_count, 1);
				}
				// Ignore other bases (N, or other variants)
				return;
			}
			read_pos += len;
			current_ref_pos += len;
		} else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) {
			// Insertion/soft clip: only read advances
			read_pos += len;
		} else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
			// Deletion/skip: only ref advances
			if (ref_pos >= current_ref_pos && ref_pos < current_ref_pos + len) {
				// Position is in deletion, no base to count
				// Continue to next SNP position rather than returning
				return;
			}
			current_ref_pos += len;
		} else if (op == BAM_CHARD_CLIP) {
			// Hard clip: nothing to do
		}
	}
}

// Worker function for parallel position checking
static void worker_check_positions(void *data, long i, int tid)
{
	step_t *s = (step_t*)data;
	bam1_t *b = s->bam_recs[i];
	
	// Skip unmapped or low quality reads
	if (b->core.flag & (BAM_FUNMAP | BAM_FQCFAIL | BAM_FDUP)) return;
	
	int ref_tid = b->core.tid;
	int ref_start = b->core.pos;
	int ref_end = bam_endpos(b);
	
	// Check only SNPs on the same chromosome and within read boundaries
	// This is more efficient than checking all SNPs
	int j;
	for (j = 0; j < s->p->db->n; ++j) {
		pattern_t *pat = &s->p->db->a[j];
		
		// Skip SNPs on different chromosomes
		if (pat->tid != ref_tid) continue;
		
		// Skip SNPs outside read boundaries
		if (pat->pos < ref_start || pat->pos >= ref_end) continue;
		
		// This SNP overlaps the read, check the base
		count_base_at_position(b, pat->pos, pat->ref, pat->alt, 
		                       &pat->ref_count, &pat->alt_count);
	}
}

// Pipeline worker function
static void *worker_pipeline(void *data, int step, void *in)
{
	pipeline_t *p = (pipeline_t*)data;
	
	if (step == 0) { // Step 1: read BAM records from regions
		step_t *s;
		bam1_t *b;
		int ret, batch_size = 1000;
		
		s = (step_t*)calloc(1, sizeof(step_t));
		s->p = p;
		s->m = batch_size;
		s->bam_recs = (bam1_t**)malloc(batch_size * sizeof(bam1_t*));
		
		// If regions available, use indexed access; otherwise fall back to sequential
		if (p->regions && p->idx) {
			// Read from current iterator or move to next region
			for (s->n = 0; s->n < batch_size; ++s->n) {
				b = bam_init1();
				
				// Try to read from current iterator
				while (1) {
					if (p->itr) {
						ret = sam_itr_next(p->fp, p->itr, b);
						if (ret >= 0) break; // Got a read
						// Current region exhausted, move to next
						hts_itr_destroy(p->itr);
						p->itr = NULL;
					}
					
					// Move to next region
					if (p->curr_region >= p->regions->n) {
						// All regions processed
						ret = -1;
						break;
					}
					
					region_t *r = &p->regions->a[p->curr_region++];
					int tid = sam_hdr_name2tid(p->hdr, r->chr);
					if (tid < 0) {
						fprintf(stderr, "Warning: chromosome %s not found in BAM\n", r->chr);
						continue;
					}
					
					p->itr = sam_itr_queryi(p->idx, tid, r->start, r->end);
					if (!p->itr) {
						fprintf(stderr, "Warning: failed to create iterator for %s:%d-%d\n",
						        r->chr, r->start, r->end);
						continue;
					}
				}
				
				if (ret < 0) {
					bam_destroy1(b);
					break;
				}
				s->bam_recs[s->n] = b;
			}
		} else {
			// Fall back to sequential reading (original behavior)
			for (s->n = 0; s->n < batch_size; ++s->n) {
				b = bam_init1();
				ret = sam_read1(p->fp, p->hdr, b);
				if (ret < 0) {
					bam_destroy1(b);
					break;
				}
				s->bam_recs[s->n] = b;
			}
		}
		
		if (s->n == 0) {
			free(s->bam_recs);
			free(s);
			return 0;
		}
		return s;
		
	} else if (step == 1) { // Step 2: check positions in parallel
		step_t *s = (step_t*)in;
		
		// Process all reads in parallel
		kt_for(p->n_thread, worker_check_positions, s, s->n);
		
		// Clean up
		for (int i = 0; i < s->n; ++i) {
			bam_destroy1(s->bam_recs[i]);
		}
		free(s->bam_recs);
		free(s);
	}
	
	return 0;
}

// Count variants in BAM file with multi-threading
void count_bam_variants(const char *fn, int n_thread,
                        snp_map_t *snp_map, pattern_db_t *db,
                        region_list_t *regions)
{
	pipeline_t pl;
	
	pl.fp = sam_open(fn, "r");
	if (!pl.fp) {
		fprintf(stderr, "Error: failed to open BAM file: %s\n", fn);
		return;
	}
	
	pl.hdr = sam_hdr_read(pl.fp);
	if (!pl.hdr) {
		fprintf(stderr, "Error: failed to read BAM header\n");
		sam_close(pl.fp);
		return;
	}
	
	// Load BAM index for random access
	pl.idx = sam_index_load(pl.fp, fn);
	if (!pl.idx || !regions) {
		if (!pl.idx) {
			fprintf(stderr, "[M::%s] Warning: failed to load BAM index for %s, processing all reads\n", __func__, fn);
		}
		if (!regions) {
			fprintf(stderr, "[M::%s] Warning: no target regions available, processing all reads\n", __func__);
		}
		// Fall back to sequential processing without index
		pl.regions = NULL;
		pl.itr = NULL;
		pl.curr_region = 0;
		if (pl.idx) {
			hts_idx_destroy(pl.idx);
			pl.idx = NULL;
		}
	} else {
		fprintf(stderr, "[M::%s] Using indexed access to fetch reads from %d target regions\n", __func__, regions->n);
		pl.regions = regions;
		pl.itr = NULL;
		pl.curr_region = 0;
	}
	
	pl.n_thread = n_thread;
	pl.snp_map = snp_map;
	pl.db = db;
	
	kt_pipeline(2, worker_pipeline, &pl, 2);
	
	if (pl.itr) hts_itr_destroy(pl.itr);
	if (pl.idx) hts_idx_destroy(pl.idx);
	sam_hdr_destroy(pl.hdr);
	sam_close(pl.fp);
}

int main(int argc, char *argv[])
{
	int c, i, n_thread = 4;
	char *pattern_fn = 0, *out_fn = 0;
	pattern_db_t *db;
	snp_map_t *snp_map = 0;
	FILE *out_fp;
	uint64_t total_ref = 0, total_alt = 0;
	double avg_depth;
	ketopt_t o = KETOPT_INIT;
	samFile *fp;
	sam_hdr_t *hdr;
	
	while ((c = ketopt(&o, argc, argv, 1, "p:o:t:", 0)) >= 0) {
		if (c == 'p') pattern_fn = o.arg;
		else if (c == 'o') out_fn = o.arg;
		else if (c == 't') n_thread = atoi(o.arg);
	}
	
	if (!pattern_fn || !out_fn || argc - o.ind < 1) {
		fprintf(stderr, "Usage: bam-vaf-counter [options] -p <patterns.txt> -o <output.vaf> <reads.bam> [reads2.bam ...]\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -p FILE   input pattern file\n");
		fprintf(stderr, "  -o FILE   output VAF file\n");
		fprintf(stderr, "  -t INT    number of threads [%d]\n", n_thread);
		fprintf(stderr, "\nNote: This version directly counts ref/alt bases at SNP positions (no k-mer extraction).\n");
		fprintf(stderr, "      It is much faster than k-mer-based counting.\n");
		return 1;
	}
	
	fprintf(stderr, "[M::%s] Loading patterns...\n", __func__);
	db = load_patterns(pattern_fn);
	if (!db) {
		fprintf(stderr, "Error: failed to load pattern file\n");
		return 1;
	}
	fprintf(stderr, "[M::%s] Loaded %d patterns\n", __func__, db->n);
	
	// Open first BAM to get header for chromosome mapping
	fprintf(stderr, "[M::%s] Reading BAM header...\n", __func__);
	fp = sam_open(argv[o.ind], "r");
	if (!fp) {
		fprintf(stderr, "Error: failed to open BAM file: %s\n", argv[o.ind]);
		return 1;
	}
	hdr = sam_hdr_read(fp);
	if (!hdr) {
		fprintf(stderr, "Error: failed to read BAM header\n");
		sam_close(fp);
		return 1;
	}
	
	fprintf(stderr, "[M::%s] Creating SNP position map...\n", __func__);
	snp_map = create_snp_map(db, hdr);
	
	sam_hdr_destroy(hdr);
	sam_close(fp);
	
	fprintf(stderr, "[M::%s] Building target regions from patterns...\n", __func__);
	region_list_t *regions = build_regions(db);
	if (!regions) {
		fprintf(stderr, "Error: failed to build regions\n");
		return 1;
	}
	fprintf(stderr, "[M::%s] Built %d target regions (merged from %d patterns)\n", __func__, regions->n, db->n);
	
	fprintf(stderr, "[M::%s] Counting variants in BAM files with %d threads...\n", __func__, n_thread);
	for (i = o.ind; i < argc; ++i) {
		fprintf(stderr, "[M::%s] Processing %s...\n", __func__, argv[i]);
		count_bam_variants(argv[i], n_thread, snp_map, db, regions);
	}
	
	// Calculate total counts
	for (i = 0; i < db->n; ++i) {
		total_ref += db->a[i].ref_count;
		total_alt += db->a[i].alt_count;
	}
	avg_depth = (double)(total_ref + total_alt) / (db->n > 0 ? db->n : 1);
	
	fprintf(stderr, "[M::%s] Writing VAF file...\n", __func__);
	out_fp = fopen(out_fn, "w");
	if (!out_fp) {
		fprintf(stderr, "Error: failed to open output file\n");
		return 1;
	}
	
	fprintf(out_fp, "# Average depth: %.2f\n", avg_depth);
	fprintf(out_fp, "CHR\tPOS\tRSID\tREF\tALT\tREF_COUNT\tALT_COUNT\tTOTAL_COUNT\tVAF\n");
	
	for (i = 0; i < db->n; ++i) {
		pattern_t *p = &db->a[i];
		uint32_t total = p->ref_count + p->alt_count;
		double vaf = total > 0 ? (double)p->alt_count / total : 0.0;
		fprintf(out_fp, "%s\t%d\t%s\t%c\t%c\t%u\t%u\t%u\t%.4f\n",
		        p->chr, p->pos, p->rsid, p->ref, p->alt,
		        p->ref_count, p->alt_count, total, vaf);
	}
	
	fclose(out_fp);
	fprintf(stderr, "[M::%s] Done. Average depth: %.2f\n", __func__, avg_depth);
	
	region_list_destroy(regions);
	snp_map_destroy(snp_map);
	pattern_db_destroy(db);
	
	return 0;
}
