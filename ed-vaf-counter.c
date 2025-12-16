/*
 * ed-vaf-counter: Variant Allele Frequency counter using edlib for approximate matching
 * 
 * This program takes a different approach from vaf-counter:
 * Instead of extracting all k-mers from FASTQs and looking them up in a hash table,
 * it searches for each pattern k-mer in the FASTQ reads using edlib's approximate
 * string matching. This allows for mismatches and can be more efficient when the
 * number of patterns is small relative to the total number of k-mers in the reads.
 * 
 * Usage: ed-vaf-counter [options] -p <patterns.txt> -o <output.vaf> <reads.fq> [...]
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <zlib.h>
#include "ketopt.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#include "edlib.h"

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

// Load pattern file
pattern_db_t *load_patterns(const char *fn)
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

// Search for k-mer pattern in a sequence using edlib
// Returns number of matches found
static int search_kmer_in_seq(const char *kmer, const char *seq, int seq_len, int max_edit_dist)
{
	EdlibAlignResult result;
	EdlibAlignConfig config;
	int count = 0;
	int kmer_len = strlen(kmer);
	
	// Configure edlib for infix alignment (k-mer anywhere in sequence)
	// with specified maximum edit distance
	config = edlibNewAlignConfig(max_edit_dist, EDLIB_MODE_HW, EDLIB_TASK_LOC, NULL, 0);
	
	// Perform alignment
	result = edlibAlign(kmer, kmer_len, seq, seq_len, config);
	
	if (result.status == EDLIB_STATUS_OK && result.editDistance >= 0) {
		// Count all locations where alignment was found
		// For exact matches (editDistance == 0), we count each occurrence
		// For approximate matches, we also count based on number of locations found
		count = result.numLocations;
	}
	
	// Free the result
	edlibFreeAlignResult(result);
	
	return count;
}

// Count k-mers in a FASTQ file by searching for pattern k-mers
static void count_fastq_kmers(const char *fn, pattern_db_t *db, int max_edit_dist)
{
	gzFile fp;
	kseq_t *ks;
	int ret;
	int i;
	
	if ((fp = gzopen(fn, "r")) == 0) {
		fprintf(stderr, "Warning: failed to open %s\n", fn);
		return;
	}
	
	ks = kseq_init(fp);
	
	// Read each sequence from the FASTQ file
	while ((ret = kseq_read(ks)) >= 0) {
		const char *seq = ks->seq.s;
		int seq_len = ks->seq.l;
		
		// For each pattern, search for both ref and alt k-mers
		for (i = 0; i < db->n; ++i) {
			int ref_matches = search_kmer_in_seq(db->a[i].ref_kmer, seq, seq_len, max_edit_dist);
			int alt_matches = search_kmer_in_seq(db->a[i].alt_kmer, seq, seq_len, max_edit_dist);
			
			// Update counts (using atomic add to be thread-safe if we add threading later)
			__sync_fetch_and_add(&db->a[i].ref_count, ref_matches);
			__sync_fetch_and_add(&db->a[i].alt_count, alt_matches);
		}
	}
	
	kseq_destroy(ks);
	gzclose(fp);
}

int main(int argc, char *argv[])
{
	int c, i;
	int max_edit_dist = 0;  // Default: exact match only
	char *pattern_fn = 0, *out_fn = 0;
	pattern_db_t *db;
	FILE *out_fp;
	uint64_t total_ref = 0, total_alt = 0;
	double avg_depth;
	ketopt_t o = KETOPT_INIT;
	
	while ((c = ketopt(&o, argc, argv, 1, "p:o:e:", 0)) >= 0) {
		if (c == 'p') pattern_fn = o.arg;
		else if (c == 'o') out_fn = o.arg;
		else if (c == 'e') max_edit_dist = atoi(o.arg);
	}
	
	if (!pattern_fn || !out_fn || argc - o.ind < 1) {
		fprintf(stderr, "Usage: ed-vaf-counter [options] -p <patterns.txt> -o <output.vaf> <reads.fq> [reads2.fq ...]\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -p FILE   input pattern file\n");
		fprintf(stderr, "  -o FILE   output VAF file\n");
		fprintf(stderr, "  -e INT    maximum edit distance for approximate matching [%d]\n", max_edit_dist);
		fprintf(stderr, "\nDescription:\n");
		fprintf(stderr, "  This program uses edlib to search for pattern k-mers in FASTQ reads.\n");
		fprintf(stderr, "  Unlike vaf-counter which extracts all k-mers from reads and looks them up,\n");
		fprintf(stderr, "  ed-vaf-counter searches for each pattern k-mer in the reads using approximate\n");
		fprintf(stderr, "  string matching. This can be more efficient for small pattern sets.\n");
		fprintf(stderr, "  Set -e 0 for exact matches only (default), or higher values to allow mismatches.\n");
		return 1;
	}
	
	fprintf(stderr, "[M::%s] Loading patterns...\n", __func__);
	db = load_patterns(pattern_fn);
	if (!db) {
		fprintf(stderr, "Error: failed to load pattern file\n");
		return 1;
	}
	fprintf(stderr, "[M::%s] Loaded %d patterns\n", __func__, db->n);
	
	fprintf(stderr, "[M::%s] Searching for k-mers in FASTQ files (max edit distance: %d)...\n", 
	        __func__, max_edit_dist);
	for (i = o.ind; i < argc; ++i) {
		fprintf(stderr, "[M::%s] Processing %s...\n", __func__, argv[i]);
		count_fastq_kmers(argv[i], db, max_edit_dist);
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
		pattern_db_destroy(db);
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
	
	pattern_db_destroy(db);
	
	return 0;
}
