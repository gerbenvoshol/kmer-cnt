#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ketopt.h"

#define MAX_SAMPLES 1000
#define MAX_SNPS 100000
#define MAX_LINE 4096

typedef struct {
	char name[256];
	double *vaf;    // VAF values for each SNP
	int *depth;     // Depth (total count) for each SNP
	int n_snps;
} sample_t;

typedef struct {
	char chr[256];
	int pos;
	char rsid[256];
} snp_info_t;

// Read VAF file and populate sample data
int load_vaf_file(const char *fn, sample_t *sample, snp_info_t *snps, int *n_snps)
{
	FILE *fp;
	char line[MAX_LINE];
	int i = 0;
	
	fp = fopen(fn, "r");
	if (!fp) return 0;
	
	// Extract sample name from filename
	const char *base = strrchr(fn, '/');
	if (base) base++;
	else base = fn;
	strncpy(sample->name, base, 255);
	sample->name[255] = '\0';
	
	// Remove .vaf extension if present
	char *ext = strstr(sample->name, ".vaf");
	if (ext) *ext = '\0';
	
	sample->vaf = (double*)malloc(MAX_SNPS * sizeof(double));
	sample->depth = (int*)malloc(MAX_SNPS * sizeof(int));
	if (!sample->vaf || !sample->depth) {
		free(sample->vaf);
		free(sample->depth);
		fclose(fp);
		return 0;
	}
	sample->n_snps = 0;
	
	while (fgets(line, MAX_LINE, fp)) {
		if (line[0] == '#') continue; // skip header/comment
		if (strncmp(line, "CHR", 3) == 0) continue; // skip column header
		
		char chr[256], rsid[256], ref, alt;
		int pos, ref_count, alt_count, total_count;
		double vaf;
		
		if (sscanf(line, "%255s\t%d\t%255s\t%c\t%c\t%d\t%d\t%d\t%lf",
		           chr, &pos, rsid, &ref, &alt,
		           &ref_count, &alt_count, &total_count, &vaf) == 9) {
			
			// Check bounds
			if (i >= MAX_SNPS) {
				fprintf(stderr, "Warning: too many SNPs (max %d), truncating\n", MAX_SNPS);
				break;
			}
			
			// Store SNP info (only once, from first file)
			if (*n_snps == i) {
				strcpy(snps[i].chr, chr);
				snps[i].pos = pos;
				strcpy(snps[i].rsid, rsid);
				(*n_snps)++;
			}
			
			sample->vaf[i] = vaf;
			sample->depth[i] = total_count;
			i++;
		}
	}
	
	sample->n_snps = i;
	fclose(fp);
	return 1;
}

// Calculate depth-aware Pearson correlation coefficient
// Only include SNPs with sufficient depth (>=1) in both samples
double pearson_correlation_depth_aware(double *x, int *depth_x, double *y, int *depth_y, int n, int min_snps)
{
	int i, valid_count = 0;
	double sum_x = 0, sum_y = 0, sum_xy = 0;
	double sum_x2 = 0, sum_y2 = 0;
	double mean_x, mean_y, numerator, denom_x, denom_y;
	
	// First pass: count valid SNPs with sufficient depth
	for (i = 0; i < n; ++i) {
		if (depth_x[i] >= 1 && depth_y[i] >= 1) {
			valid_count++;
		}
	}
	
	// Require minimum number of valid SNPs (default 20, like NGSCheckMate)
	if (valid_count < min_snps) return 0.0;
	
	// Calculate means using only valid SNPs
	for (i = 0; i < n; ++i) {
		if (depth_x[i] >= 1 && depth_y[i] >= 1) {
			sum_x += x[i];
			sum_y += y[i];
		}
	}
	mean_x = sum_x / valid_count;
	mean_y = sum_y / valid_count;
	
	// Calculate correlation using only valid SNPs
	for (i = 0; i < n; ++i) {
		if (depth_x[i] >= 1 && depth_y[i] >= 1) {
			double dx = x[i] - mean_x;
			double dy = y[i] - mean_y;
			sum_xy += dx * dy;
			sum_x2 += dx * dx;
			sum_y2 += dy * dy;
		}
	}
	
	numerator = sum_xy;
	denom_x = sqrt(sum_x2);
	denom_y = sqrt(sum_y2);
	
	// Handle division by zero with small epsilon (like NGSCheckMate)
	if (denom_x < 1e-10 || denom_y < 1e-10) {
		return numerator / (sqrt(sum_x2 * sum_y2) + 0.00001);
	}
	
	return numerator / (denom_x * denom_y);
}

// Calculate distance matrix using depth-aware Pearson correlation
void calculate_correlation_matrix(sample_t *samples, int n_samples, double **corr_matrix, int min_snps)
{
	int i, j;
	
	for (i = 0; i < n_samples; ++i) {
		corr_matrix[i][i] = 1.0; // self-correlation is 1
		for (j = i + 1; j < n_samples; ++j) {
			double r = pearson_correlation_depth_aware(
				samples[i].vaf, samples[i].depth,
				samples[j].vaf, samples[j].depth,
				samples[i].n_snps, min_snps);
			corr_matrix[i][j] = r;
			corr_matrix[j][i] = r; // symmetric
		}
	}
}

// Simple UPGMA clustering for tree building
typedef struct {
	int left, right;  // -1 means leaf (sample), >=0 means internal node
	double distance;
	int sample_id;    // for leaves
} tree_node_t;

// Find minimum distance in matrix (excluding diagonal)
void find_min_distance(double **dist, int *active, int n, int *min_i, int *min_j, double *min_dist)
{
	int i, j;
	*min_dist = 1e10;
	*min_i = -1;
	*min_j = -1;
	
	for (i = 0; i < n; ++i) {
		if (!active[i]) continue;
		for (j = i + 1; j < n; ++j) {
			if (!active[j]) continue;
			if (dist[i][j] < *min_dist) {
				*min_dist = dist[i][j];
				*min_i = i;
				*min_j = j;
			}
		}
	}
}

// Build simple tree using UPGMA-like algorithm
void build_tree(sample_t *samples, int n_samples, double **corr_matrix, FILE *tree_fp)
{
	int i, j, k;
	double **dist;
	int *active;
	
	// Convert correlation to distance (1 - r)
	dist = (double**)malloc(n_samples * sizeof(double*));
	if (!dist) return;
	for (i = 0; i < n_samples; ++i) {
		dist[i] = (double*)malloc(n_samples * sizeof(double));
		if (!dist[i]) {
			for (j = 0; j < i; ++j) free(dist[j]);
			free(dist);
			return;
		}
		for (j = 0; j < n_samples; ++j) {
			dist[i][j] = 1.0 - corr_matrix[i][j];
		}
	}
	
	active = (int*)malloc(n_samples * sizeof(int));
	if (!active) {
		for (i = 0; i < n_samples; ++i) free(dist[i]);
		free(dist);
		return;
	}
	for (i = 0; i < n_samples; ++i) active[i] = 1;
	
	fprintf(tree_fp, "# Simple dendrogram (UPGMA-like clustering)\n");
	fprintf(tree_fp, "# Format: (Sample1:distance, Sample2:distance)\n");
	
	int n_active = n_samples;
	while (n_active > 1) {
		int min_i, min_j;
		double min_dist;
		
		find_min_distance(dist, active, n_samples, &min_i, &min_j, &min_dist);
		
		if (min_i == -1 || min_j == -1) break;
		
		fprintf(tree_fp, "Cluster: %s (%.4f) <-> %s (%.4f)\n",
		        samples[min_i].name, min_dist/2,
		        samples[min_j].name, min_dist/2);
		
		// Update distances (average linkage)
		for (k = 0; k < n_samples; ++k) {
			if (k == min_i || k == min_j || !active[k]) continue;
			dist[min_i][k] = (dist[min_i][k] + dist[min_j][k]) / 2.0;
			dist[k][min_i] = dist[min_i][k];
		}
		
		// Deactivate min_j
		active[min_j] = 0;
		n_active--;
	}
	
	// Free memory
	for (i = 0; i < n_samples; ++i) free(dist[i]);
	free(dist);
	free(active);
}

int main(int argc, char *argv[])
{
	int c, i, j;
	char *out_fn = 0;
	int build_tree_flag = 0;
	int min_snps = 20;  // Minimum SNPs required for correlation (like NGSCheckMate)
	sample_t *samples;
	snp_info_t *snps;
	int n_samples, n_snps = 0;
	double **corr_matrix;
	FILE *out_fp, *tree_fp = NULL;
	ketopt_t o = KETOPT_INIT;
	
	while ((c = ketopt(&o, argc, argv, 1, "o:tm:", 0)) >= 0) {
		if (c == 'o') out_fn = o.arg;
		else if (c == 't') build_tree_flag = 1;
		else if (c == 'm') min_snps = atoi(o.arg);
	}
	
	n_samples = argc - o.ind;
	
	if (!out_fn || n_samples < 2) {
		fprintf(stderr, "Usage: correlation-matrix -o <output.corr> [-t] [-m INT] <sample1.vaf> <sample2.vaf> [sample3.vaf ...]\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -o FILE   output correlation matrix file\n");
		fprintf(stderr, "  -t        build tree/dendrogram (outputs to <output.tree>)\n");
		fprintf(stderr, "  -m INT    minimum SNPs with depth >= 1 required for correlation [%d]\n", min_snps);
		return 1;
	}
	
	fprintf(stderr, "[M::%s] Loading %d VAF files...\n", __func__, n_samples);
	
	samples = (sample_t*)calloc(n_samples, sizeof(sample_t));
	snps = (snp_info_t*)calloc(MAX_SNPS, sizeof(snp_info_t));
	if (!samples || !snps) {
		fprintf(stderr, "Error: failed to allocate memory\n");
		return 1;
	}
	
	for (i = 0; i < n_samples; ++i) {
		if (!load_vaf_file(argv[o.ind + i], &samples[i], snps, &n_snps)) {
			fprintf(stderr, "Error: failed to load %s\n", argv[o.ind + i]);
			return 1;
		}
		fprintf(stderr, "[M::%s] Loaded %s: %d SNPs\n", __func__, samples[i].name, samples[i].n_snps);
	}
	
	fprintf(stderr, "[M::%s] Computing correlation matrix...\n", __func__);
	
	// Allocate correlation matrix
	corr_matrix = (double**)malloc(n_samples * sizeof(double*));
	if (!corr_matrix) {
		fprintf(stderr, "Error: failed to allocate correlation matrix\n");
		return 1;
	}
	for (i = 0; i < n_samples; ++i) {
		corr_matrix[i] = (double*)malloc(n_samples * sizeof(double));
		if (!corr_matrix[i]) {
			fprintf(stderr, "Error: failed to allocate correlation matrix row\n");
			return 1;
		}
	}
	
	calculate_correlation_matrix(samples, n_samples, corr_matrix, min_snps);
	
	// Write correlation matrix
	fprintf(stderr, "[M::%s] Writing correlation matrix...\n", __func__);
	out_fp = fopen(out_fn, "w");
	if (!out_fp) {
		fprintf(stderr, "Error: failed to open output file\n");
		return 1;
	}
	
	// Write header
	fprintf(out_fp, "Sample");
	for (i = 0; i < n_samples; ++i) {
		fprintf(out_fp, "\t%s", samples[i].name);
	}
	fprintf(out_fp, "\n");
	
	// Write matrix
	for (i = 0; i < n_samples; ++i) {
		fprintf(out_fp, "%s", samples[i].name);
		for (j = 0; j < n_samples; ++j) {
			fprintf(out_fp, "\t%.6f", corr_matrix[i][j]);
		}
		fprintf(out_fp, "\n");
	}
	
	fclose(out_fp);
	fprintf(stderr, "[M::%s] Correlation matrix written to %s\n", __func__, out_fn);
	
	// Build tree if requested
	if (build_tree_flag) {
		char tree_fn[512];
		snprintf(tree_fn, sizeof(tree_fn), "%s", out_fn);
		char *ext = strstr(tree_fn, ".corr");
		if (ext) strcpy(ext, ".tree");
		else strcat(tree_fn, ".tree");
		
		fprintf(stderr, "[M::%s] Building dendrogram...\n", __func__);
		tree_fp = fopen(tree_fn, "w");
		if (tree_fp) {
			build_tree(samples, n_samples, corr_matrix, tree_fp);
			fclose(tree_fp);
			fprintf(stderr, "[M::%s] Dendrogram written to %s\n", __func__, tree_fn);
		}
	}
	
	// Free memory
	for (i = 0; i < n_samples; ++i) {
		free(samples[i].vaf);
		free(samples[i].depth);
		free(corr_matrix[i]);
	}
	free(corr_matrix);
	free(samples);
	free(snps);
	
	return 0;
}
