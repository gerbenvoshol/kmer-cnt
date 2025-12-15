#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <zlib.h>
#include "ketopt.h"
#include "kthread.h"

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#include "khashl.h"

// Fast hash function for k-mers
// K-mers are already well-distributed due to canonical representation,
// so we can use a simple multiplicative hash instead of the complex kh_hash_uint64
static inline khint_t kmer_hash(uint64_t key) {
	// Simple multiplicative hash with good avalanche properties
	// Uses MurmurHash3-inspired finalizer for fast mixing
	key ^= key >> 33;
	key *= 0xff51afd7ed558ccdULL;
	key ^= key >> 33;
	return (khint_t)key;
}

// Use cached hash variant for better performance during collision resolution
// This stores the hash value in each bucket to avoid recomputing during probing
KHASHL_CMAP_INIT(, kmer_cnt_t, kmer_cnt, uint64_t, uint32_t, kmer_hash, kh_eq_generic)

#define KMER_BUF_GROWTH_FACTOR 1.2
#define KMER_REF_FLAG 0  // Flag for reference k-mer in combined map
#define KMER_ALT_FLAG 1  // Flag for alternative k-mer in combined map

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

// Create combined hash table mapping k-mer to (pattern_index << 1) | is_alt
// This allows both ref and alt k-mers in a single map
// Note: Assumes k-mers are unique (each k-mer appears in only one pattern)
// This should be guaranteed by snp-pattern-gen which filters for uniqueness
kmer_cnt_t *create_combined_kmer_map(pattern_db_t *db, int k)
{
	kmer_cnt_t *h;
	int i, absent, n_collisions = 0;
	khint_t itr;
	
	// Check for pattern index overflow (max index that can be safely encoded)
	if (db->n > (INT32_MAX >> 1)) {
		fprintf(stderr, "Error: too many patterns (%d), maximum is %d\n", 
		        db->n, INT32_MAX >> 1);
		return NULL;
	}
	
	h = kmer_cnt_init();
	// Pre-allocate hash table to avoid frequent resizing
	// We'll have 2 k-mers per pattern (ref and alt)
	// Allocate db->n * 3 to provide sufficient space with typical hash table load factors
	kmer_cnt_cm_resize(h, db->n * 3);
	
	for (i = 0; i < db->n; ++i) {
		uint64_t kmer;
		
		// Add reference k-mer
		kmer = encode_kmer(db->a[i].ref_kmer, k);
		if (kmer != UINT64_MAX) {
			uint64_t can = canonical_kmer(kmer, k);
			itr = kmer_cnt_put(h, can, &absent);
			if (absent) {
				kh_val(h, itr) = (i << 1) | KMER_REF_FLAG; // Encode pattern index and ref flag
			} else {
				++n_collisions;
			}
		}
		
		// Add alternative k-mer
		kmer = encode_kmer(db->a[i].alt_kmer, k);
		if (kmer != UINT64_MAX) {
			uint64_t can = canonical_kmer(kmer, k);
			itr = kmer_cnt_put(h, can, &absent);
			if (absent) {
				kh_val(h, itr) = (i << 1) | KMER_ALT_FLAG; // Encode pattern index and alt flag
			} else {
				++n_collisions;
			}
		}
	}
	
	if (n_collisions > 0) {
		fprintf(stderr, "[W::%s] Warning: %d k-mer collisions detected. "
		        "Some patterns may have overlapping k-mers.\n", __func__, n_collisions);
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
				// Ensure buffer has space; use branch hint since reallocation is rare
				if (__builtin_expect(buf->n == buf->m, 0)) {
					buf->m = buf->m < 8 ? 8 : buf->m + (buf->m >> 1);
					buf->a = (uint64_t*)realloc(buf->a, buf->m * sizeof(uint64_t));
				}
				uint64_t y = x[0] < x[1] ? x[0] : x[1];
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
	int k, block_len, n_thread;
	kseq_t *ks;
	kmer_cnt_t *kmer_map;  // Combined map for both ref and alt k-mers
	pattern_db_t *db;
} pipeline_t;

// Data for each pipeline step
typedef struct {
	pipeline_t *p;
	int n, m, sum_len, nk;
	int *len;
	char **seq;
	kmer_buf_t buf;
} step_t;

// Worker function for parallel k-mer lookup
static void worker_lookup(void *data, long i, int tid)
{
	step_t *s = (step_t*)data;
	uint64_t kmer = s->buf.a[i];
	khint_t itr;
	
	// Check combined k-mer map (single lookup instead of two)
	itr = kmer_cnt_get(s->p->kmer_map, kmer);
	if (itr != kh_end(s->p->kmer_map)) {
		uint32_t val = kh_val(s->p->kmer_map, itr);
		int idx = val >> 1;        // Extract pattern index
		int is_alt = val & 1;      // Extract ref/alt flag
		
		if (is_alt) {
			__sync_fetch_and_add(&s->p->db->a[idx].alt_count, 1);
		} else {
			__sync_fetch_and_add(&s->p->db->a[idx].ref_count, 1);
		}
	}
}

// Pipeline worker function
static void *worker_pipeline(void *data, int step, void *in)
{
	pipeline_t *p = (pipeline_t*)data;
	
	if (step == 0) { // Step 1: read sequences
		int ret;
		step_t *s;
		s = (step_t*)calloc(1, sizeof(step_t));
		s->p = p;
		
		while ((ret = kseq_read(p->ks)) >= 0) {
			int l = p->ks->seq.l;
			if (l < p->k) continue;
			
			if (s->n == s->m) {
				s->m = s->m < 16 ? 16 : s->m + (s->m >> 1);
				s->len = (int*)realloc(s->len, s->m * sizeof(int));
				s->seq = (char**)realloc(s->seq, s->m * sizeof(char*));
			}
			s->seq[s->n] = (char*)malloc(l);
			memcpy(s->seq[s->n], p->ks->seq.s, l);
			s->len[s->n++] = l;
			s->sum_len += l;
			s->nk += l - p->k + 1;
			
			if (s->sum_len >= p->block_len)
				break;
		}
		
		if (s->sum_len == 0) {
			free(s);
			return 0;
		}
		return s;
		
	} else if (step == 1) { // Step 2: extract k-mers
		step_t *s = (step_t*)in;
		int i;
		
		s->buf.m = (int)(s->nk * KMER_BUF_GROWTH_FACTOR) + 1;
		s->buf.a = (uint64_t*)malloc(s->buf.m * sizeof(uint64_t));
		s->buf.n = 0;
		
		for (i = 0; i < s->n; ++i) {
			extract_kmers_to_buf(&s->buf, p->k, s->len[i], s->seq[i]);
			free(s->seq[i]);
		}
		free(s->seq);
		free(s->len);
		
		return s;
		
	} else if (step == 2) { // Step 3: lookup k-mers
		step_t *s = (step_t*)in;
		
		kt_for(p->n_thread, worker_lookup, s, s->buf.n);
		
		free(s->buf.a);
		free(s);
	}
	
	return 0;
}

// Count k-mers in FASTQ files with multi-threading
void count_fastq_kmers(const char *fn, int k, int n_thread, int block_size, 
                       kmer_cnt_t *kmer_map, pattern_db_t *db)
{
	pipeline_t pl;
	gzFile fp;
	
	if ((fp = gzopen(fn, "r")) == 0) return;
	
	pl.ks = kseq_init(fp);
	pl.k = k;
	pl.n_thread = n_thread;
	pl.block_len = block_size;
	pl.kmer_map = kmer_map;
	pl.db = db;
	
	kt_pipeline(3, worker_pipeline, &pl, 3);
	
	kseq_destroy(pl.ks);
	gzclose(fp);
}

int main(int argc, char *argv[])
{
	int c, k = 21, i, n_thread = 4, block_size = 10000000;
	char *pattern_fn = 0, *out_fn = 0;
	pattern_db_t *db;
	kmer_cnt_t *kmer_map;
	FILE *out_fp;
	uint64_t total_ref = 0, total_alt = 0;
	double avg_depth;
	ketopt_t o = KETOPT_INIT;
	
	while ((c = ketopt(&o, argc, argv, 1, "k:p:o:t:b:", 0)) >= 0) {
		if (c == 'k') k = atoi(o.arg);
		else if (c == 'p') pattern_fn = o.arg;
		else if (c == 'o') out_fn = o.arg;
		else if (c == 't') n_thread = atoi(o.arg);
		else if (c == 'b') block_size = atoi(o.arg);
	}
	
	if (!pattern_fn || !out_fn || argc - o.ind < 1) {
		fprintf(stderr, "Usage: vaf-counter [options] -p <patterns.txt> -o <output.vaf> <reads.fq> [reads2.fq ...]\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -k INT    k-mer length [%d]\n", k);
		fprintf(stderr, "  -p FILE   input pattern file\n");
		fprintf(stderr, "  -o FILE   output VAF file\n");
		fprintf(stderr, "  -t INT    number of threads [%d]\n", n_thread);
		fprintf(stderr, "  -b INT    block size [%d]\n", block_size);
		return 1;
	}
	
	fprintf(stderr, "[M::%s] Loading patterns...\n", __func__);
	db = load_patterns(pattern_fn, k);
	if (!db) {
		fprintf(stderr, "Error: failed to load pattern file\n");
		return 1;
	}
	fprintf(stderr, "[M::%s] Loaded %d patterns\n", __func__, db->n);
	
	fprintf(stderr, "[M::%s] Creating k-mer map...\n", __func__);
	kmer_map = create_combined_kmer_map(db, k);
	if (!kmer_map) {
		fprintf(stderr, "Error: failed to create k-mer map\n");
		pattern_db_destroy(db);
		return 1;
	}
	
	fprintf(stderr, "[M::%s] Counting k-mers in FASTQ files with %d threads...\n", __func__, n_thread);
	for (i = o.ind; i < argc; ++i) {
		fprintf(stderr, "[M::%s] Processing %s...\n", __func__, argv[i]);
		count_fastq_kmers(argv[i], k, n_thread, block_size, kmer_map, db);
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
	
	kmer_cnt_destroy(kmer_map);
	pattern_db_destroy(db);
	
	return 0;
}
