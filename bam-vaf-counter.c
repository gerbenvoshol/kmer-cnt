/*
 * bam-vaf-counter: Count k-mer occurrences from BAM files and calculate VAF
 * 
 * This program reads aligned sequencing data from BAM files and counts
 * reference and alternative k-mers from a pattern file to generate VAF files.
 * 
 * OPTIMIZATION: Uses indexed BAM access to fetch only reads overlapping SNP
 * positions, dramatically reducing the number of reads processed compared to
 * sequential scanning. Regions are extracted from the pattern file and merged
 * to minimize redundant reads. Falls back to sequential reading if BAM index
 * (.bai) is not available.
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
KHASHL_MAP_INIT(, kmer_cnt_t, kmer_cnt, uint64_t, uint32_t, kh_hash_uint64, kh_eq_generic)

#define KMER_BUF_GROWTH_FACTOR 1.2

const unsigned char seq_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

typedef struct {
	char chr[256];
	int start;
	int end;
	char rsid[256];
	char ref;
	char alt;
	char ref_kmer[128];
	char alt_kmer[128];
	uint32_t ref_count;
	uint32_t alt_count;
} pattern_t;

typedef struct {
	int n, m;
	pattern_t *a;
} pattern_db_t;

// Buffer for storing k-mers
typedef struct {
	int n, m;
	uint64_t *a;
} kmer_buf_t;

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

// Encode k-mer to 64-bit integer
uint64_t encode_kmer(const char *seq, int k)
{
	int i;
	uint64_t kmer = 0;
	for (i = 0; i < k; ++i) {
		int c = seq_nt4_table[(uint8_t)seq[i]];
		if (c >= 4) return UINT64_MAX;
		kmer = (kmer << 2) | c;
	}
	return kmer;
}

// Reverse complement
uint64_t revcomp_kmer(uint64_t kmer, int k)
{
	int i;
	uint64_t rc = 0;
	for (i = 0; i < k; ++i) {
		rc = (rc << 2) | (3 - (kmer & 3));
		kmer >>= 2;
	}
	return rc;
}

// Get canonical k-mer
uint64_t canonical_kmer(uint64_t kmer, int k)
{
	uint64_t rc = revcomp_kmer(kmer, k);
	return kmer < rc ? kmer : rc;
}

// Load pattern file
pattern_db_t *load_patterns(const char *fn, int k)
{
	FILE *fp;
	pattern_db_t *db;
	pattern_t pat;
	
	fp = fopen(fn, "r");
	if (!fp) return 0;
	
	db = (pattern_db_t*)calloc(1, sizeof(pattern_db_t));
	if (!db) {
		fclose(fp);
		return 0;
	}
	
	while (fscanf(fp, "%255s%d%d%255s %c %c%127s%127s", 
	              pat.chr, &pat.start, &pat.end, pat.rsid, 
	              &pat.ref, &pat.alt, pat.ref_kmer, pat.alt_kmer) == 8) {
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
region_list_t *build_regions(pattern_db_t *db, int k)
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
	
	// Create initial regions with flanking bases
	// Each k-mer spans k bases, so we need to read k-1 bases on each side of SNP
	for (i = 0; i < db->n; ++i) {
		strcpy(regions->a[i].chr, db->a[i].chr);
		regions->a[i].start = db->a[i].start - k + 1;
		if (regions->a[i].start < 0) regions->a[i].start = 0;
		regions->a[i].end = db->a[i].end + k - 1;
	}
	regions->n = db->n;
	
	// Sort regions by chromosome and start position using qsort
	qsort(regions->a, regions->n, sizeof(region_t), region_cmp);
	
	// Merge overlapping or adjacent regions on same chromosome
	int write_idx = 0;
	for (i = 1; i < regions->n; ++i) {
		if (strcmp(regions->a[write_idx].chr, regions->a[i].chr) == 0 &&
		    regions->a[write_idx].end + 1 >= regions->a[i].start) {
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

// Create hash table mapping k-mer to pattern index
kmer_cnt_t *create_kmer_map(pattern_db_t *db, int k, int is_ref)
{
	kmer_cnt_t *h;
	int i, absent;
	khint_t itr;
	
	h = kmer_cnt_init();
	
	for (i = 0; i < db->n; ++i) {
		char *kmer_str = is_ref ? db->a[i].ref_kmer : db->a[i].alt_kmer;
		uint64_t kmer = encode_kmer(kmer_str, k);
		if (kmer != UINT64_MAX) {
			uint64_t can = canonical_kmer(kmer, k);
			itr = kmer_cnt_put(h, can, &absent);
			if (absent) kh_val(h, itr) = i; // store pattern index
		}
	}
	
	return h;
}

// Extract k-mers from a sequence into a buffer
static void extract_kmers_to_buf(kmer_buf_t *buf, int k, int len, const char *seq)
{
	int i, l;
	uint64_t x[2], mask = (1ULL << k*2) - 1, shift = (k - 1) * 2;
	
	for (i = l = 0, x[0] = x[1] = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) {
			x[0] = (x[0] << 2 | c) & mask;
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;
			if (++l >= k) {
				uint64_t y = x[0] < x[1] ? x[0] : x[1];
				if (buf->n == buf->m) {
					buf->m = buf->m < 8 ? 8 : buf->m + (buf->m >> 1);
					buf->a = (uint64_t*)realloc(buf->a, buf->m * sizeof(uint64_t));
				}
				buf->a[buf->n++] = y;
			}
		} else {
			l = 0;
			x[0] = x[1] = 0;
		}
	}
}

// Global data for pipeline threading
typedef struct {
	int k, n_thread;
	samFile *fp;
	sam_hdr_t *hdr;
	hts_idx_t *idx;
	hts_itr_t *itr;
	kmer_cnt_t *ref_map;
	kmer_cnt_t *alt_map;
	pattern_db_t *db;
	region_list_t *regions;
	int curr_region;
} pipeline_t;

// Data for each pipeline step
typedef struct {
	pipeline_t *p;
	int n, m;
	bam1_t **bam_recs;
	kmer_buf_t buf;
} step_t;

// Worker function for parallel k-mer lookup
static void worker_lookup(void *data, long i, int tid)
{
	step_t *s = (step_t*)data;
	uint64_t kmer = s->buf.a[i];
	khint_t itr;
	
	// Check reference k-mers
	itr = kmer_cnt_get(s->p->ref_map, kmer);
	if (itr != kh_end(s->p->ref_map)) {
		int idx = kh_val(s->p->ref_map, itr);
		__sync_fetch_and_add(&s->p->db->a[idx].ref_count, 1);
	}
	
	// Check alternative k-mers
	itr = kmer_cnt_get(s->p->alt_map, kmer);
	if (itr != kh_end(s->p->alt_map)) {
		int idx = kh_val(s->p->alt_map, itr);
		__sync_fetch_and_add(&s->p->db->a[idx].alt_count, 1);
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
		
	} else if (step == 1) { // Step 2: extract k-mers
		step_t *s = (step_t*)in;
		int i, estimated_kmers = 0;
		
		// Estimate total k-mers
		for (i = 0; i < s->n; ++i) {
			int len = s->bam_recs[i]->core.l_qseq;
			if (len >= p->k) estimated_kmers += len - p->k + 1;
		}
		
		s->buf.m = (int)(estimated_kmers * KMER_BUF_GROWTH_FACTOR) + 1;
		s->buf.a = (uint64_t*)malloc(s->buf.m * sizeof(uint64_t));
		s->buf.n = 0;
		
		for (i = 0; i < s->n; ++i) {
			bam1_t *b = s->bam_recs[i];
			int len = b->core.l_qseq;
			if (len < p->k) continue;
			
			// Get sequence from BAM record
			uint8_t *seq = bam_get_seq(b);
			char *seq_str = (char*)malloc(len + 1);
			int j;
			for (j = 0; j < len; ++j) {
				seq_str[j] = seq_nt16_str[bam_seqi(seq, j)];
			}
			seq_str[len] = '\0';
			
			extract_kmers_to_buf(&s->buf, p->k, len, seq_str);
			free(seq_str);
			bam_destroy1(b);
		}
		free(s->bam_recs);
		
		return s;
		
	} else if (step == 2) { // Step 3: lookup k-mers
		step_t *s = (step_t*)in;
		
		kt_for(p->n_thread, worker_lookup, s, s->buf.n);
		
		free(s->buf.a);
		free(s);
	}
	
	return 0;
}

// Count k-mers in BAM file with multi-threading
void count_bam_kmers(const char *fn, int k, int n_thread,
                     kmer_cnt_t *ref_map, kmer_cnt_t *alt_map, pattern_db_t *db,
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
	
	pl.k = k;
	pl.n_thread = n_thread;
	pl.ref_map = ref_map;
	pl.alt_map = alt_map;
	pl.db = db;
	
	kt_pipeline(3, worker_pipeline, &pl, 3);
	
	if (pl.itr) hts_itr_destroy(pl.itr);
	if (pl.idx) hts_idx_destroy(pl.idx);
	sam_hdr_destroy(pl.hdr);
	sam_close(pl.fp);
}

int main(int argc, char *argv[])
{
	int c, k = 21, i, n_thread = 4;
	char *pattern_fn = 0, *out_fn = 0;
	pattern_db_t *db;
	kmer_cnt_t *ref_map, *alt_map;
	FILE *out_fp;
	uint64_t total_ref = 0, total_alt = 0;
	double avg_depth;
	ketopt_t o = KETOPT_INIT;
	
	while ((c = ketopt(&o, argc, argv, 1, "k:p:o:t:", 0)) >= 0) {
		if (c == 'k') k = atoi(o.arg);
		else if (c == 'p') pattern_fn = o.arg;
		else if (c == 'o') out_fn = o.arg;
		else if (c == 't') n_thread = atoi(o.arg);
	}
	
	if (!pattern_fn || !out_fn || argc - o.ind < 1) {
		fprintf(stderr, "Usage: bam-vaf-counter [options] -p <patterns.txt> -o <output.vaf> <reads.bam> [reads2.bam ...]\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -k INT    k-mer length [%d]\n", k);
		fprintf(stderr, "  -p FILE   input pattern file\n");
		fprintf(stderr, "  -o FILE   output VAF file\n");
		fprintf(stderr, "  -t INT    number of threads [%d]\n", n_thread);
		return 1;
	}
	
	fprintf(stderr, "[M::%s] Loading patterns...\n", __func__);
	db = load_patterns(pattern_fn, k);
	if (!db) {
		fprintf(stderr, "Error: failed to load pattern file\n");
		return 1;
	}
	fprintf(stderr, "[M::%s] Loaded %d patterns\n", __func__, db->n);
	
	fprintf(stderr, "[M::%s] Creating k-mer maps...\n", __func__);
	ref_map = create_kmer_map(db, k, 1);
	alt_map = create_kmer_map(db, k, 0);
	
	fprintf(stderr, "[M::%s] Building target regions from patterns...\n", __func__);
	region_list_t *regions = build_regions(db, k);
	if (!regions) {
		fprintf(stderr, "Error: failed to build regions\n");
		return 1;
	}
	fprintf(stderr, "[M::%s] Built %d target regions (merged from %d patterns)\n", __func__, regions->n, db->n);
	
	fprintf(stderr, "[M::%s] Counting k-mers in BAM files with %d threads...\n", __func__, n_thread);
	for (i = o.ind; i < argc; ++i) {
		fprintf(stderr, "[M::%s] Processing %s...\n", __func__, argv[i]);
		count_bam_kmers(argv[i], k, n_thread, ref_map, alt_map, db, regions);
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
		        p->chr, p->start, p->rsid, p->ref, p->alt,
		        p->ref_count, p->alt_count, total, vaf);
	}
	
	fclose(out_fp);
	fprintf(stderr, "[M::%s] Done. Average depth: %.2f\n", __func__, avg_depth);
	
	region_list_destroy(regions);
	kmer_cnt_destroy(ref_map);
	kmer_cnt_destroy(alt_map);
	pattern_db_destroy(db);
	
	return 0;
}
