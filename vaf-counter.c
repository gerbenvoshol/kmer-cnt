#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <zlib.h>
#include <time.h>
#include <sys/time.h>
#include "ketopt.h"
#include "kthread.h"

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#include "khashl.h"

// SIMD headers
#ifdef __x86_64__
#include <immintrin.h>
#endif

/*
 * Performance optimizations applied:
 * 1. Buffer growth factor increased to 1.5 (from 1.2) to reduce realloc frequency
 * 2. Pre-allocate exact buffer size in pipeline step 2 to avoid reallocation
 * 3. Use __ATOMIC_RELAXED for counter updates (reduces memory ordering overhead)
 * 4. SIMD optimizations for k-mer extraction using SSSE3/SSE4.1
 * 5. Software prefetching to reduce cache misses
 * 6. Verbose performance reporting option
 */

// Performance statistics structure
typedef struct {
	double time_pattern_load;
	double time_kmer_map_create;
	double time_kmer_counting;
	double time_output_write;
	uint64_t total_kmers_extracted;
	uint64_t total_sequences_processed;
	uint64_t total_bases_processed;
	int verbose;
} perf_stats_t;

// Global performance stats
static perf_stats_t g_perf_stats = {0};

// Timing helper functions
static inline double get_time(void) {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + tv.tv_usec * 1e-6;
}

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

#define KMER_BUF_GROWTH_FACTOR 1.5  // Increased from 1.2 for fewer reallocations
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
	// This reduces the number of rehashing operations during insertion
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

// SIMD-optimized sequence encoding using PSHUFB-based lookup table
// This is significantly faster than the comparison-based approach
#ifdef __x86_64__

#ifdef __SSSE3__
// Fast SIMD encoding using PSHUFB lookup table (SSSE3+)
// This approach is much faster than comparisons + blendv
static inline void encode_seq_simd_ssse3(const char *seq, int len, uint8_t *encoded) {
	int i;
	
	// Create lookup table for PSHUFB
	// Maps ASCII values to nucleotide codes (0-3) or 4 (invalid)
	// We use the low 4 bits of each ASCII character as index
	// A/a (0x41/0x61) -> low nibble 0x1 -> maps to 0
	// C/c (0x43/0x63) -> low nibble 0x3 -> maps to 1  
	// T/t (0x54/0x74) -> low nibble 0x4 -> maps to 3
	// U/u (0x55/0x75) -> low nibble 0x5 -> maps to 3
	// G/g (0x47/0x67) -> low nibble 0x7 -> maps to 2
	const __m128i lut = _mm_setr_epi8(
		4, 0, 4, 1,  3, 3, 4, 2,  // 0x0-0x7: invalid, A, invalid, C, T, U, invalid, G
		4, 4, 4, 4,  4, 4, 4, 4   // 0x8-0xF: all invalid
	);
	const __m128i mask = _mm_set1_epi8(0x0F);
	
	for (i = 0; i + 16 <= len; i += 16) {
		__m128i seq_vec = _mm_loadu_si128((__m128i*)(seq + i));
		// Extract low nibble for lookup
		__m128i nibble = _mm_and_si128(seq_vec, mask);
		// Lookup encoding
		__m128i result = _mm_shuffle_epi8(lut, nibble);
		_mm_storeu_si128((__m128i*)(encoded + i), result);
	}
	
	// Handle remaining bytes
	for (; i < len; ++i) {
		encoded[i] = seq_nt4_table[(uint8_t)seq[i]];
	}
}
#endif

#ifdef __SSE4_1__
// SSE4.1 version using blendv (fallback for non-SSSE3)
static inline void encode_seq_simd_sse41(const char *seq, int len, uint8_t *encoded) {
	int i;
	__m128i lut_A = _mm_set1_epi8('A');
	__m128i lut_C = _mm_set1_epi8('C');
	__m128i lut_G = _mm_set1_epi8('G');
	__m128i lut_T = _mm_set1_epi8('T');
	__m128i lut_a = _mm_set1_epi8('a');
	__m128i lut_c = _mm_set1_epi8('c');
	__m128i lut_g = _mm_set1_epi8('g');
	__m128i lut_t = _mm_set1_epi8('t');
	
	for (i = 0; i + 16 <= len; i += 16) {
		__m128i seq_vec = _mm_loadu_si128((__m128i*)(seq + i));
		__m128i result = _mm_set1_epi8(4); // Default to 4 (invalid)
		
		// Compare with each nucleotide
		__m128i cmp_A = _mm_or_si128(_mm_cmpeq_epi8(seq_vec, lut_A), _mm_cmpeq_epi8(seq_vec, lut_a));
		__m128i cmp_C = _mm_or_si128(_mm_cmpeq_epi8(seq_vec, lut_C), _mm_cmpeq_epi8(seq_vec, lut_c));
		__m128i cmp_G = _mm_or_si128(_mm_cmpeq_epi8(seq_vec, lut_G), _mm_cmpeq_epi8(seq_vec, lut_g));
		__m128i cmp_T = _mm_or_si128(_mm_cmpeq_epi8(seq_vec, lut_T), _mm_cmpeq_epi8(seq_vec, lut_t));
		
		// Set encoded values
		result = _mm_blendv_epi8(result, _mm_set1_epi8(0), cmp_A);
		result = _mm_blendv_epi8(result, _mm_set1_epi8(1), cmp_C);
		result = _mm_blendv_epi8(result, _mm_set1_epi8(2), cmp_G);
		result = _mm_blendv_epi8(result, _mm_set1_epi8(3), cmp_T);
		
		_mm_storeu_si128((__m128i*)(encoded + i), result);
	}
	
	// Handle remaining bytes
	for (; i < len; ++i) {
		encoded[i] = seq_nt4_table[(uint8_t)seq[i]];
	}
}
#endif

// Dispatcher function that chooses best available SIMD implementation
static inline void encode_seq_simd(const char *seq, int len, uint8_t *encoded) {
#ifdef __SSSE3__
	encode_seq_simd_ssse3(seq, len, encoded);
#elif defined(__SSE4_1__)
	encode_seq_simd_sse41(seq, len, encoded);
#else
	// Scalar fallback
	for (int i = 0; i < len; ++i) {
		encoded[i] = seq_nt4_table[(uint8_t)seq[i]];
	}
#endif
}
#endif

// Extract k-mers from a sequence into a buffer
static void extract_kmers_to_buf(kmer_buf_t *buf, int k, int len, const char *seq)
{
	int i, l;
	uint64_t x[2], mask = (1ULL << k*2) - 1, shift = (k - 1) * 2;
	
#ifdef __x86_64__
	// Use SIMD for encoding sequence on x86_64
	uint8_t *encoded = (uint8_t*)malloc(len);
	if (encoded) {
		encode_seq_simd(seq, len, encoded);
	} else {
		// Memory allocation failed, log warning and fall back to scalar
		if (g_perf_stats.verbose) {
			fprintf(stderr, "[W::%s] Failed to allocate encoding buffer (%d bytes), using scalar fallback\n", 
			        __func__, len);
		}
	}
	if (encoded) {
		
		for (i = l = 0, x[0] = x[1] = 0; i < len; ++i) {
			int c = encoded[i];
			if (c < 4) {
				x[0] = (x[0] << 2 | c) & mask;
				x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;
				if (++l >= k) {
					uint64_t y = x[0] < x[1] ? x[0] : x[1];
					// Ensure buffer has space; use branch hint since reallocation is rare
					if (__builtin_expect(buf->n == buf->m, 0)) {
						int new_m = buf->m < 8 ? 8 : buf->m + (buf->m >> 1);
						uint64_t *new_a = (uint64_t*)realloc(buf->a, new_m * sizeof(uint64_t));
						if (!new_a) {
							// Realloc failed; keep existing buffer but stop adding k-mers
							free(encoded);
							return;
						}
						buf->a = new_a;
						buf->m = new_m;
					}
					buf->a[buf->n++] = y;
					g_perf_stats.total_kmers_extracted++;
				}
			} else {
				l = 0;
				x[0] = x[1] = 0;
			}
		}
		free(encoded);
		return;
	}
#endif
	
	// Fallback to scalar implementation
	for (i = l = 0, x[0] = x[1] = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) {
			x[0] = (x[0] << 2 | c) & mask;
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;
			if (++l >= k) {
				uint64_t y = x[0] < x[1] ? x[0] : x[1];
				// Ensure buffer has space; use branch hint since reallocation is rare
				if (__builtin_expect(buf->n == buf->m, 0)) {
					int new_m = buf->m < 8 ? 8 : buf->m + (buf->m >> 1);
					uint64_t *new_a = (uint64_t*)realloc(buf->a, new_m * sizeof(uint64_t));
					if (!new_a) {
						// Realloc failed; keep existing buffer but stop adding k-mers
						return;
					}
					buf->a = new_a;
					buf->m = new_m;
				}
				buf->a[buf->n++] = y;
				g_perf_stats.total_kmers_extracted++;
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
	uint64_t total_bases;  // Track total bases for performance stats
	uint64_t total_seqs;   // Track total sequences for performance stats
} pipeline_t;

// Data for each pipeline step
typedef struct {
	pipeline_t *p;
	int n, m, sum_len, nk;
	int *len;
	char **seq;
	kmer_buf_t buf;
} step_t;

// Worker function for parallel k-mer lookup with prefetching
static void worker_lookup(void *data, long i, int tid)
{
	step_t *s = (step_t*)data;
	uint64_t kmer = s->buf.a[i];
	khint_t itr;
	
	// Prefetch next k-mers to improve cache utilization
	// This helps reduce cache misses in the tight loop
	if (i + 4 < s->buf.n) {
		__builtin_prefetch(&s->buf.a[i + 4], 0, 1);
	}
	
	// Check combined k-mer map (single lookup instead of two)
	itr = kmer_cnt_get(s->p->kmer_map, kmer);
	if (itr != kh_end(s->p->kmer_map)) {
		uint32_t val = kh_val(s->p->kmer_map, itr);
		int idx = val >> 1;        // Extract pattern index
		int is_alt = val & 1;      // Extract ref/alt flag
		
		// Prefetch the pattern data we're about to update
		__builtin_prefetch(&s->p->db->a[idx], 1, 1);
		
		// Use relaxed atomic operations for better performance
		// Memory ordering doesn't matter here since we only care about final counts
		if (is_alt) {
			__atomic_fetch_add(&s->p->db->a[idx].alt_count, 1, __ATOMIC_RELAXED);
		} else {
			__atomic_fetch_add(&s->p->db->a[idx].ref_count, 1, __ATOMIC_RELAXED);
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
			p->total_bases += l;
			p->total_seqs++;
			
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
		
		// Pre-allocate exact size needed to avoid reallocation during extraction
		s->buf.m = s->nk;
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
	double t0 = get_time();
	
	if ((fp = gzopen(fn, "r")) == 0) return;
	
	pl.ks = kseq_init(fp);
	pl.k = k;
	pl.n_thread = n_thread;
	pl.block_len = block_size;
	pl.kmer_map = kmer_map;
	pl.db = db;
	pl.total_bases = 0;
	pl.total_seqs = 0;
	
	kt_pipeline(3, worker_pipeline, &pl, 3);
	
	g_perf_stats.total_bases_processed += pl.total_bases;
	g_perf_stats.total_sequences_processed += pl.total_seqs;
	
	if (g_perf_stats.verbose) {
		double elapsed = get_time() - t0;
		fprintf(stderr, "[V::%s] Processed %s: %llu sequences, %llu bases in %.2f sec (%.2f Mbases/sec)\n",
		        __func__, fn, (unsigned long long)pl.total_seqs, (unsigned long long)pl.total_bases, 
		        elapsed, pl.total_bases / elapsed / 1e6);
	}
	
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
	double t_start, t_end;
	ketopt_t o = KETOPT_INIT;
	
	g_perf_stats.verbose = 0;
	
	while ((c = ketopt(&o, argc, argv, 1, "k:p:o:t:b:v", 0)) >= 0) {
		if (c == 'k') k = atoi(o.arg);
		else if (c == 'p') pattern_fn = o.arg;
		else if (c == 'o') out_fn = o.arg;
		else if (c == 't') n_thread = atoi(o.arg);
		else if (c == 'b') block_size = atoi(o.arg);
		else if (c == 'v') g_perf_stats.verbose = 1;
	}
	
	if (!pattern_fn || !out_fn || argc - o.ind < 1) {
		fprintf(stderr, "Usage: vaf-counter [options] -p <patterns.txt> -o <output.vaf> <reads.fq> [reads2.fq ...]\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -k INT    k-mer length [%d]\n", k);
		fprintf(stderr, "  -p FILE   input pattern file\n");
		fprintf(stderr, "  -o FILE   output VAF file\n");
		fprintf(stderr, "  -t INT    number of threads [%d]\n", n_thread);
		fprintf(stderr, "  -b INT    block size [%d]\n", block_size);
		fprintf(stderr, "  -v        verbose mode (report performance statistics)\n");
		return 1;
	}
	
	t_start = get_time();
	
	fprintf(stderr, "[M::%s] Loading patterns...\n", __func__);
	t_end = get_time();
	db = load_patterns(pattern_fn, k);
	if (!db) {
		fprintf(stderr, "Error: failed to load pattern file\n");
		return 1;
	}
	g_perf_stats.time_pattern_load = get_time() - t_end;
	fprintf(stderr, "[M::%s] Loaded %d patterns in %.3f sec\n", __func__, db->n, g_perf_stats.time_pattern_load);
	
	fprintf(stderr, "[M::%s] Creating k-mer map...\n", __func__);
	t_end = get_time();
	kmer_map = create_combined_kmer_map(db, k);
	if (!kmer_map) {
		fprintf(stderr, "Error: failed to create k-mer map\n");
		pattern_db_destroy(db);
		return 1;
	}
	g_perf_stats.time_kmer_map_create = get_time() - t_end;
	if (g_perf_stats.verbose) {
		fprintf(stderr, "[V::%s] Created k-mer map with %d entries in %.3f sec\n", 
		        __func__, kh_size(kmer_map), g_perf_stats.time_kmer_map_create);
	}
	
	fprintf(stderr, "[M::%s] Counting k-mers in FASTQ files with %d threads...\n", __func__, n_thread);
	t_end = get_time();
	for (i = o.ind; i < argc; ++i) {
		fprintf(stderr, "[M::%s] Processing %s...\n", __func__, argv[i]);
		count_fastq_kmers(argv[i], k, n_thread, block_size, kmer_map, db);
	}
	g_perf_stats.time_kmer_counting = get_time() - t_end;
	
	// Calculate total counts
	for (i = 0; i < db->n; ++i) {
		total_ref += db->a[i].ref_count;
		total_alt += db->a[i].alt_count;
	}
	avg_depth = (double)(total_ref + total_alt) / (db->n > 0 ? db->n : 1);
	
	fprintf(stderr, "[M::%s] Writing VAF file...\n", __func__);
	t_end = get_time();
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
	g_perf_stats.time_output_write = get_time() - t_end;
	
	fprintf(stderr, "[M::%s] Done. Average depth: %.2f\n", __func__, avg_depth);
	
	// Print performance statistics if verbose mode enabled
	if (g_perf_stats.verbose) {
		double total_time = get_time() - t_start;
		fprintf(stderr, "\n=== Performance Statistics ===\n");
		fprintf(stderr, "Total runtime:           %.3f sec\n", total_time);
		fprintf(stderr, "  Pattern loading:       %.3f sec (%.1f%%)\n", 
		        g_perf_stats.time_pattern_load, 100.0 * g_perf_stats.time_pattern_load / total_time);
		fprintf(stderr, "  K-mer map creation:    %.3f sec (%.1f%%)\n", 
		        g_perf_stats.time_kmer_map_create, 100.0 * g_perf_stats.time_kmer_map_create / total_time);
		fprintf(stderr, "  K-mer counting:        %.3f sec (%.1f%%)\n", 
		        g_perf_stats.time_kmer_counting, 100.0 * g_perf_stats.time_kmer_counting / total_time);
		fprintf(stderr, "  Output writing:        %.3f sec (%.1f%%)\n", 
		        g_perf_stats.time_output_write, 100.0 * g_perf_stats.time_output_write / total_time);
		fprintf(stderr, "\nThroughput:\n");
		fprintf(stderr, "  Sequences processed:   %llu\n", (unsigned long long)g_perf_stats.total_sequences_processed);
		fprintf(stderr, "  Bases processed:       %llu (%.2f Mbases)\n", 
		        (unsigned long long)g_perf_stats.total_bases_processed,
		        g_perf_stats.total_bases_processed / 1e6);
		fprintf(stderr, "  K-mers extracted:      %llu (%.2f million)\n",
		        (unsigned long long)g_perf_stats.total_kmers_extracted,
		        g_perf_stats.total_kmers_extracted / 1e6);
		if (g_perf_stats.time_kmer_counting > 0) {
			fprintf(stderr, "  Speed:                 %.2f Mbases/sec\n",
			        g_perf_stats.total_bases_processed / g_perf_stats.time_kmer_counting / 1e6);
			fprintf(stderr, "  K-mer throughput:      %.2f million k-mers/sec\n",
			        g_perf_stats.total_kmers_extracted / g_perf_stats.time_kmer_counting / 1e6);
		}
		fprintf(stderr, "\nMemory:\n");
		fprintf(stderr, "  Patterns:              %d\n", db->n);
		fprintf(stderr, "  Hash table entries:    %d\n", kh_size(kmer_map));
		fprintf(stderr, "  Hash table capacity:   %d\n", kh_capacity(kmer_map));
		fprintf(stderr, "  Hash table load:       %.1f%%\n", 
		        100.0 * kh_size(kmer_map) / kh_capacity(kmer_map));
		
		// Report SIMD capabilities
		fprintf(stderr, "\nOptimizations:\n");
#ifdef __SSSE3__
		fprintf(stderr, "  SIMD:                  SSSE3 enabled (optimized PSHUFB)\n");
#elif defined(__SSE4_1__)
		fprintf(stderr, "  SIMD:                  SSE4.1 enabled (blendv)\n");
#elif defined(__SSE2__)
		fprintf(stderr, "  SIMD:                  SSE2 enabled (basic)\n");
#else
		fprintf(stderr, "  SIMD:                  Not available (scalar fallback)\n");
#endif
		fprintf(stderr, "  Threads:               %d workers\n", n_thread);
		fprintf(stderr, "==============================\n");
	}
	
	kmer_cnt_destroy(kmer_map);
	pattern_db_destroy(db);
	
	return 0;
}
