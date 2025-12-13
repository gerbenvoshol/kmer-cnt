/*
 * vcf-vaf-counter: Generate VAF files from VCF files
 * 
 * This program reads genotype information from VCF files and generates
 * VAF files compatible with the correlation-matrix and match-classifier tools.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "ketopt.h"

#include "htslib/htslib/vcf.h"
#include "htslib/htslib/hts.h"

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

// Find pattern by position
int find_pattern(pattern_db_t *db, const char *chr, int pos)
{
	int i;
	for (i = 0; i < db->n; ++i) {
		if (strcmp(db->a[i].chr, chr) == 0 && db->a[i].start == pos) {
			return i;
		}
	}
	return -1;
}

// Process VCF file and extract genotype information
void process_vcf(const char *fn, pattern_db_t *db, int sample_idx, int min_depth)
{
	htsFile *fp;
	bcf_hdr_t *hdr;
	bcf1_t *rec;
	int ngt_arr = 0, ndp_arr = 0, nad_arr = 0;
	int *gt_arr = NULL, *dp_arr = NULL, *ad_arr = NULL;
	
	fp = hts_open(fn, "r");
	if (!fp) {
		fprintf(stderr, "Error: failed to open VCF file: %s\n", fn);
		return;
	}
	
	hdr = bcf_hdr_read(fp);
	if (!hdr) {
		fprintf(stderr, "Error: failed to read VCF header\n");
		hts_close(fp);
		return;
	}
	
	rec = bcf_init1();
	
	while (bcf_read(fp, hdr, rec) == 0) {
		bcf_unpack(rec, BCF_UN_ALL);
		
		// Get chromosome and position
		const char *chr = bcf_hdr_id2name(hdr, rec->rid);
		int pos = rec->pos;
		
		// Find matching pattern
		int pat_idx = find_pattern(db, chr, pos);
		if (pat_idx < 0) continue;
		
		// Check if this is a SNP
		if (rec->n_allele != 2) continue;
		if (strlen(rec->d.allele[0]) != 1 || strlen(rec->d.allele[1]) != 1) continue;
		
		// Verify alleles match pattern
		if (rec->d.allele[0][0] != db->a[pat_idx].ref ||
		    rec->d.allele[1][0] != db->a[pat_idx].alt) continue;
		
		// Get genotype (GT)
		int ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
		if (ngt <= 0) continue;
		
		int sample_offset = sample_idx * 2;
		if (sample_offset + 1 >= ngt) continue;
		
		int allele1 = bcf_gt_allele(gt_arr[sample_offset]);
		int allele2 = bcf_gt_allele(gt_arr[sample_offset + 1]);
		
		// Skip missing genotypes
		if (allele1 < 0 || allele2 < 0) continue;
		
		// Get depth (DP or AD)
		int depth = 0;
		int ref_depth = 0, alt_depth = 0;
		
		// Try to get allele depth (AD)
		int nad = bcf_get_format_int32(hdr, rec, "AD", &ad_arr, &nad_arr);
		if (nad > 0) {
			int n_alleles = rec->n_allele;
			int offset = sample_idx * n_alleles;
			if (offset + 1 < nad && ad_arr[offset] != bcf_int32_missing && 
			    ad_arr[offset + 1] != bcf_int32_missing) {
				ref_depth = ad_arr[offset];
				alt_depth = ad_arr[offset + 1];
				depth = ref_depth + alt_depth;
			}
		}
		
		// If AD not available, try DP and estimate from genotype
		if (depth == 0) {
			int ndp = bcf_get_format_int32(hdr, rec, "DP", &dp_arr, &ndp_arr);
			if (ndp > 0 && sample_idx < ndp && dp_arr[sample_idx] != bcf_int32_missing) {
				depth = dp_arr[sample_idx];
				
				// Estimate allele counts from genotype
				if (allele1 == 0 && allele2 == 0) {
					// 0/0 - homozygous reference
					ref_depth = depth;
					alt_depth = 0;
				} else if (allele1 == 1 && allele2 == 1) {
					// 1/1 - homozygous alternative
					ref_depth = 0;
					alt_depth = depth;
				} else {
					// 0/1 or 1/0 - heterozygous (assume 50/50)
					ref_depth = depth / 2;
					alt_depth = depth - ref_depth;
				}
			}
		}
		
		// Apply minimum depth filter
		if (depth < min_depth) continue;
		
		// Update pattern counts
		db->a[pat_idx].ref_count = ref_depth;
		db->a[pat_idx].alt_count = alt_depth;
	}
	
	free(gt_arr);
	free(dp_arr);
	free(ad_arr);
	bcf_destroy1(rec);
	bcf_hdr_destroy(hdr);
	hts_close(fp);
}

int main(int argc, char *argv[])
{
	int c, i, sample_idx = 0, min_depth = 1;
	char *pattern_fn = 0, *out_fn = 0, *vcf_fn = 0;
	pattern_db_t *db;
	FILE *out_fp;
	uint64_t total_ref = 0, total_alt = 0;
	double avg_depth;
	ketopt_t o = KETOPT_INIT;
	
	while ((c = ketopt(&o, argc, argv, 1, "p:o:v:s:d:", 0)) >= 0) {
		if (c == 'p') pattern_fn = o.arg;
		else if (c == 'o') out_fn = o.arg;
		else if (c == 'v') vcf_fn = o.arg;
		else if (c == 's') sample_idx = atoi(o.arg);
		else if (c == 'd') min_depth = atoi(o.arg);
	}
	
	if (!pattern_fn || !out_fn || !vcf_fn) {
		fprintf(stderr, "Usage: vcf-vaf-counter [options] -p <patterns.txt> -v <input.vcf> -o <output.vaf>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -p FILE   input pattern file\n");
		fprintf(stderr, "  -v FILE   input VCF/BCF file\n");
		fprintf(stderr, "  -o FILE   output VAF file\n");
		fprintf(stderr, "  -s INT    sample index (0-based) [%d]\n", sample_idx);
		fprintf(stderr, "  -d INT    minimum depth [%d]\n", min_depth);
		return 1;
	}
	
	fprintf(stderr, "[M::%s] Loading patterns...\n", __func__);
	db = load_patterns(pattern_fn);
	if (!db) {
		fprintf(stderr, "Error: failed to load pattern file\n");
		return 1;
	}
	fprintf(stderr, "[M::%s] Loaded %d patterns\n", __func__, db->n);
	
	fprintf(stderr, "[M::%s] Processing VCF file...\n", __func__);
	process_vcf(vcf_fn, db, sample_idx, min_depth);
	
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
	
	pattern_db_destroy(db);
	
	return 0;
}
