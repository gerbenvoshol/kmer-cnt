#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <zlib.h>
#include "ketopt.h"

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#include "khashl.h"
KHASHL_MAP_INIT(, kmer_t, kmer, uint64_t, uint32_t, kh_hash_uint64, kh_eq_generic)

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
	char *name;
	char *seq;
	int len;
} fasta_seq_t;

typedef struct {
	int n, m;
	fasta_seq_t *a;
} fasta_db_t;

typedef struct {
	char chr[256];
	int start;
	int end;
	char rsid[256];
	char ref;
	char alt;
} snp_t;

// Load entire FASTA file into memory
fasta_db_t *load_fasta(const char *fn)
{
	gzFile fp;
	kseq_t *ks;
	fasta_db_t *db;
	
	if ((fp = gzopen(fn, "r")) == 0) return 0;
	ks = kseq_init(fp);
	db = (fasta_db_t*)calloc(1, sizeof(fasta_db_t));
	if (!db) {
		kseq_destroy(ks);
		gzclose(fp);
		return 0;
	}
	
	while (kseq_read(ks) >= 0) {
		if (db->n == db->m) {
			fasta_seq_t *tmp;
			db->m = db->m? db->m << 1 : 16;
			tmp = (fasta_seq_t*)realloc(db->a, db->m * sizeof(fasta_seq_t));
			if (!tmp) break;
			db->a = tmp;
		}
		db->a[db->n].name = strdup(ks->name.s);
		db->a[db->n].seq = strdup(ks->seq.s);
		if (!db->a[db->n].name || !db->a[db->n].seq) break;
		db->a[db->n].len = ks->seq.l;
		++db->n;
	}
	
	kseq_destroy(ks);
	gzclose(fp);
	return db;
}

void fasta_db_destroy(fasta_db_t *db)
{
	int i;
	if (db == 0) return;
	for (i = 0; i < db->n; ++i) {
		free(db->a[i].name);
		free(db->a[i].seq);
	}
	free(db->a);
	free(db);
}

// Find sequence by chromosome name
fasta_seq_t *find_seq(fasta_db_t *db, const char *chr)
{
	int i;
	for (i = 0; i < db->n; ++i) {
		if (strcmp(db->a[i].name, chr) == 0)
			return &db->a[i];
	}
	return 0;
}

// Encode k-mer to 64-bit integer
uint64_t encode_kmer(const char *seq, int k)
{
	int i;
	uint64_t kmer = 0;
	for (i = 0; i < k; ++i) {
		int c = seq_nt4_table[(uint8_t)seq[i]];
		if (c >= 4) return UINT64_MAX; // has N
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

// Count k-mer occurrences in entire genome
void count_genome_kmers(fasta_db_t *db, int k, kmer_t *h)
{
	int i, j, l;
	uint64_t x[2], mask = (1ULL << k*2) - 1, shift = (k - 1) * 2;
	
	for (i = 0; i < db->n; ++i) {
		char *seq = db->a[i].seq;
		int len = db->a[i].len;
		
		for (j = l = 0, x[0] = x[1] = 0; j < len; ++j) {
			int c = seq_nt4_table[(uint8_t)seq[j]];
			if (c < 4) {
				x[0] = (x[0] << 2 | c) & mask;
				x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;
				if (++l >= k) {
					khint_t itr;
					int absent;
					uint64_t y = x[0] < x[1] ? x[0] : x[1];
					itr = kmer_put(h, y, &absent);
					if (absent) kh_val(h, itr) = 0;
					++kh_val(h, itr);
				}
			} else {
				l = 0;
				x[0] = x[1] = 0;
			}
		}
	}
}

// Extract k-mer around SNP position
int extract_snp_kmer(fasta_seq_t *seq, int pos, char alt, int k, char *ref_kmer, char *alt_kmer)
{
	int flank = k / 2;
	int start = pos - flank;
	int i;
	
	if (start < 0 || start + k > seq->len) return 0;
	
	// Check for N bases
	for (i = 0; i < k; ++i) {
		if (seq_nt4_table[(uint8_t)seq->seq[start + i]] >= 4)
			return 0;
	}
	
	// Extract reference k-mer
	memcpy(ref_kmer, seq->seq + start, k);
	ref_kmer[k] = '\0';
	
	// Extract alternative k-mer
	memcpy(alt_kmer, seq->seq + start, k);
	alt_kmer[flank] = alt;
	alt_kmer[k] = '\0';
	
	return 1;
}

int main(int argc, char *argv[])
{
	int c, k = 21;
	char *bed_fn = 0, *fasta_fn = 0, *out_fn = 0;
	FILE *bed_fp, *out_fp;
	fasta_db_t *db;
	kmer_t *h;
	snp_t snp;
	char ref_kmer[128], alt_kmer[128];
	int n_total = 0, n_unique = 0;
	ketopt_t o = KETOPT_INIT;
	
	while ((c = ketopt(&o, argc, argv, 1, "k:b:f:o:", 0)) >= 0) {
		if (c == 'k') k = atoi(o.arg);
		else if (c == 'b') bed_fn = o.arg;
		else if (c == 'f') fasta_fn = o.arg;
		else if (c == 'o') out_fn = o.arg;
	}
	
	if (k % 2 == 0) {
		fprintf(stderr, "Error: k must be odd\n");
		return 1;
	}
	
	if (!bed_fn || !fasta_fn || !out_fn) {
		fprintf(stderr, "Usage: snp-pattern-gen -k %d -b <snps.bed> -f <ref.fa> -o <patterns.txt>\n", k);
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -k INT    k-mer length (must be odd) [%d]\n", k);
		fprintf(stderr, "  -b FILE   input BED file with SNPs\n");
		fprintf(stderr, "  -f FILE   input reference genome FASTA file\n");
		fprintf(stderr, "  -o FILE   output pattern file\n");
		return 1;
	}
	
	fprintf(stderr, "[M::%s] Loading reference genome...\n", __func__);
	db = load_fasta(fasta_fn);
	if (!db) {
		fprintf(stderr, "Error: failed to load FASTA file\n");
		return 1;
	}
	fprintf(stderr, "[M::%s] Loaded %d sequences\n", __func__, db->n);
	
	fprintf(stderr, "[M::%s] Counting genome k-mers...\n", __func__);
	h = kmer_init();
	count_genome_kmers(db, k, h);
	fprintf(stderr, "[M::%s] Found %ld unique k-mers in genome\n", __func__, (long)kh_size(h));
	
	bed_fp = fopen(bed_fn, "r");
	if (!bed_fp) {
		fprintf(stderr, "Error: failed to open BED file\n");
		return 1;
	}
	
	out_fp = fopen(out_fn, "w");
	if (!out_fp) {
		fprintf(stderr, "Error: failed to open output file\n");
		return 1;
	}
	
	fprintf(stderr, "[M::%s] Processing SNPs...\n", __func__);
	while (fscanf(bed_fp, "%255s%d%d%255s %c %c", snp.chr, &snp.start, &snp.end, snp.rsid, &snp.ref, &snp.alt) == 6) {
		fasta_seq_t *seq = find_seq(db, snp.chr);
		++n_total;
		
		if (!seq) {
			fprintf(stderr, "Warning: chromosome %s not found\n", snp.chr);
			continue;
		}
		
		if (extract_snp_kmer(seq, snp.start, snp.alt, k, ref_kmer, alt_kmer)) {
			uint64_t ref_enc = encode_kmer(ref_kmer, k);
			uint64_t alt_enc = encode_kmer(alt_kmer, k);
			
			if (ref_enc == UINT64_MAX || alt_enc == UINT64_MAX) continue;
			
			uint64_t ref_can = canonical_kmer(ref_enc, k);
			uint64_t alt_can = canonical_kmer(alt_enc, k);
			
			khint_t ref_itr = kmer_get(h, ref_can);
			khint_t alt_itr = kmer_get(h, alt_can);
			
			// Check if both k-mers are unique (occur exactly once)
			if (ref_itr != kh_end(h) && kh_val(h, ref_itr) == 1 &&
			    alt_itr == kh_end(h)) {
				fprintf(out_fp, "%s\t%d\t%d\t%s\t%c\t%c\t%s\t%s\n",
				        snp.chr, snp.start, snp.end, snp.rsid, snp.ref, snp.alt,
				        ref_kmer, alt_kmer);
				++n_unique;
			}
		}
	}
	
	fprintf(stderr, "[M::%s] Total SNPs: %d, Unique k-mer pairs: %d\n", __func__, n_total, n_unique);
	
	fclose(bed_fp);
	fclose(out_fp);
	kmer_destroy(h);
	fasta_db_destroy(db);
	
	return 0;
}
