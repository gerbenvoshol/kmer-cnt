#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <zlib.h>
#include "ketopt.h"

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#include "khashl.h"
KHASHL_MAP_INIT(, kmer_cnt_t, kmer_cnt, uint64_t, uint32_t, kh_hash_uint64, kh_eq_generic)

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

// Count k-mers in FASTQ files
void count_fastq_kmers(const char *fn, int k, kmer_cnt_t *ref_map, kmer_cnt_t *alt_map, pattern_db_t *db)
{
	gzFile fp;
	kseq_t *ks;
	int l;
	uint64_t x[2], mask = (1ULL << k*2) - 1, shift = (k - 1) * 2;
	
	if ((fp = gzopen(fn, "r")) == 0) return;
	ks = kseq_init(fp);
	
	while (kseq_read(ks) >= 0) {
		int i;
		char *seq = ks->seq.s;
		int len = ks->seq.l;
		
		for (i = l = 0, x[0] = x[1] = 0; i < len; ++i) {
			int c = seq_nt4_table[(uint8_t)seq[i]];
			if (c < 4) {
				x[0] = (x[0] << 2 | c) & mask;
				x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;
				if (++l >= k) {
					uint64_t y = x[0] < x[1] ? x[0] : x[1];
					khint_t itr;
					
					// Check reference k-mers
					itr = kmer_cnt_get(ref_map, y);
					if (itr != kh_end(ref_map)) {
						int idx = kh_val(ref_map, itr);
						++db->a[idx].ref_count;
					}
					
					// Check alternative k-mers
					itr = kmer_cnt_get(alt_map, y);
					if (itr != kh_end(alt_map)) {
						int idx = kh_val(alt_map, itr);
						++db->a[idx].alt_count;
					}
				}
			} else {
				l = 0;
				x[0] = x[1] = 0;
			}
		}
	}
	
	kseq_destroy(ks);
	gzclose(fp);
}

int main(int argc, char *argv[])
{
	int c, k = 21, i;
	char *pattern_fn = 0, *out_fn = 0;
	pattern_db_t *db;
	kmer_cnt_t *ref_map, *alt_map;
	FILE *out_fp;
	uint64_t total_ref = 0, total_alt = 0;
	double avg_depth;
	ketopt_t o = KETOPT_INIT;
	
	while ((c = ketopt(&o, argc, argv, 1, "k:p:o:", 0)) >= 0) {
		if (c == 'k') k = atoi(o.arg);
		else if (c == 'p') pattern_fn = o.arg;
		else if (c == 'o') out_fn = o.arg;
	}
	
	if (!pattern_fn || !out_fn || argc - o.ind < 1) {
		fprintf(stderr, "Usage: vaf-counter -k %d -p <patterns.txt> -o <output.vaf> <reads.fq> [reads2.fq ...]\n", k);
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -k INT    k-mer length [%d]\n", k);
		fprintf(stderr, "  -p FILE   input pattern file\n");
		fprintf(stderr, "  -o FILE   output VAF file\n");
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
	
	fprintf(stderr, "[M::%s] Counting k-mers in FASTQ files...\n", __func__);
	for (i = o.ind; i < argc; ++i) {
		fprintf(stderr, "[M::%s] Processing %s...\n", __func__, argv[i]);
		count_fastq_kmers(argv[i], k, ref_map, alt_map, db);
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
	
	kmer_cnt_destroy(ref_map);
	kmer_cnt_destroy(alt_map);
	pattern_db_destroy(db);
	
	return 0;
}
