#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ketopt.h"

#define MAX_SAMPLES 1000
#define MAX_LINE 4096

typedef struct {
	char name[256];
	double avg_depth;  // Average depth for this sample
} sample_info_t;

typedef struct {
	int n_samples;
	sample_info_t *samples;
	double **corr_matrix;
	int has_depth_info;
} correlation_data_t;

// Forward declarations
void free_correlation_data(correlation_data_t *data);

// Load correlation matrix file
correlation_data_t *load_correlation_matrix(const char *fn)
{
	FILE *fp;
	char line[MAX_LINE];
	correlation_data_t *data;
	int i, j;
	
	fp = fopen(fn, "r");
	if (!fp) return NULL;
	
	data = (correlation_data_t*)calloc(1, sizeof(correlation_data_t));
	if (!data) {
		fclose(fp);
		return NULL;
	}
	
	// Read header line to get sample names
	if (!fgets(line, MAX_LINE, fp)) {
		free(data);
		fclose(fp);
		return NULL;
	}
	
	// Count samples and allocate
	char *token = strtok(line, "\t\n");
	if (!token || strcmp(token, "Sample") != 0) {
		free(data);
		fclose(fp);
		return NULL;
	}
	
	// Count samples in header
	int n = 0;
	while ((token = strtok(NULL, "\t\n")) != NULL && n < MAX_SAMPLES) {
		n++;
	}
	
	data->n_samples = n;
	data->samples = (sample_info_t*)calloc(n, sizeof(sample_info_t));
	data->corr_matrix = (double**)malloc(n * sizeof(double*));
	for (i = 0; i < n; ++i) {
		data->corr_matrix[i] = (double*)malloc(n * sizeof(double));
	}
	
	if (!data->samples || !data->corr_matrix) {
		free(data->samples);
		free(data);
		fclose(fp);
		return NULL;
	}
	
	// Re-read header to get sample names
	fseek(fp, 0, SEEK_SET);
	if (!fgets(line, MAX_LINE, fp)) {
		free_correlation_data(data);
		fclose(fp);
		return NULL;
	}
	token = strtok(line, "\t\n"); // Skip "Sample"
	for (i = 0; i < n; ++i) {
		token = strtok(NULL, "\t\n");
		if (token) {
			strncpy(data->samples[i].name, token, 255);
			data->samples[i].name[255] = '\0';
		}
	}
	
	// Read correlation values
	for (i = 0; i < n; ++i) {
		if (!fgets(line, MAX_LINE, fp)) break;
		
		token = strtok(line, "\t\n"); // Sample name (skip)
		for (j = 0; j < n; ++j) {
			token = strtok(NULL, "\t\n");
			if (token) {
				data->corr_matrix[i][j] = atof(token);
			}
		}
	}
	
	fclose(fp);
	return data;
}

void free_correlation_data(correlation_data_t *data)
{
	int i;
	if (!data) return;
	
	if (data->corr_matrix) {
		for (i = 0; i < data->n_samples; ++i) {
			free(data->corr_matrix[i]);
		}
		free(data->corr_matrix);
	}
	free(data->samples);
	free(data);
}

// Load depth information from VAF files
int load_depth_info(correlation_data_t *data, int argc, char **argv, int start_idx)
{
	int i, n_loaded = 0;
	
	for (i = 0; i < data->n_samples && start_idx + i < argc; ++i) {
		FILE *fp;
		char line[MAX_LINE];
		const char *vaf_fn = argv[start_idx + i];
		
		fp = fopen(vaf_fn, "r");
		if (!fp) {
			fprintf(stderr, "Warning: could not open %s for depth info\n", vaf_fn);
			continue;
		}
		
		// Read first line to get average depth
		if (fgets(line, MAX_LINE, fp)) {
			if (line[0] == '#') {
				double depth;
				if (sscanf(line, "# Average depth: %lf", &depth) == 1) {
					data->samples[i].avg_depth = depth;
					n_loaded++;
				}
			}
		}
		fclose(fp);
	}
	
	if (n_loaded > 0) {
		fprintf(stderr, "[M::%s] Loaded depth info for %d samples\n", __func__, n_loaded);
		data->has_depth_info = 1;
		return 1;
	}
	
	return 0;
}

// NGSCheckMate predefined model parameters
typedef struct {
	double mean_matched;
	double std_matched;
	double mean_unmatched;
	double std_unmatched;
} model_params_t;

// Get predefined model parameters based on depth (NGSCheckMate approach)
model_params_t get_predefined_model(double depth, int family_mode)
{
	model_params_t params;
	
	if (family_mode) {
		// Family mode: related samples
		if (depth > 10) {
			params.mean_matched = 0.874611;
			params.std_matched = 0.022596;
			params.mean_unmatched = 0.644481;
			params.std_unmatched = 0.020908;
		} else if (depth > 5) {
			params.mean_matched = 0.785312;
			params.std_matched = 0.021318;
			params.mean_unmatched = 0.596133;
			params.std_unmatched = 0.022502;
		} else if (depth > 2) {
			params.mean_matched = 0.650299;
			params.std_matched = 0.019252;
			params.mean_unmatched = 0.5346;
			params.std_unmatched = 0.020694;
		} else if (depth > 1) {
			params.mean_matched = 0.578582;
			params.std_matched = 0.018379;
			params.mean_unmatched = 0.495017;
			params.std_unmatched = 0.021652;
		} else if (depth > 0.5) {
			params.mean_matched = 0.524757;
			params.std_matched = 0.023218;
			params.mean_unmatched = 0.465653;
			params.std_unmatched = 0.027378;
		} else {
			// Warning: depth too low
			params.mean_matched = 0.524757;
			params.std_matched = 0.023218;
			params.mean_unmatched = 0.465653;
			params.std_unmatched = 0.027378;
		}
	} else {
		// Non-family mode: unrelated samples
		if (depth > 10) {
			params.mean_matched = 0.874546;
			params.std_matched = 0.022211;
			params.mean_unmatched = 0.310549;
			params.std_unmatched = 0.060058;
		} else if (depth > 5) {
			params.mean_matched = 0.785249;
			params.std_matched = 0.021017;
			params.mean_unmatched = 0.279778;
			params.std_unmatched = 0.054104;
		} else if (depth > 2) {
			params.mean_matched = 0.650573;
			params.std_matched = 0.018699;
			params.mean_unmatched = 0.238972;
			params.std_unmatched = 0.047196;
		} else if (depth > 1) {
			params.mean_matched = 0.578386;
			params.std_matched = 0.018526;
			params.mean_unmatched = 0.222322;
			params.std_unmatched = 0.041186;
		} else if (depth > 0.5) {
			params.mean_matched = 0.529327;
			params.std_matched = 0.025785;
			params.mean_unmatched = 0.217839;
			params.std_unmatched = 0.040334;
		} else {
			// Warning: depth too low
			params.mean_matched = 0.529327;
			params.std_matched = 0.025785;
			params.mean_unmatched = 0.217839;
			params.std_unmatched = 0.040334;
		}
	}
	
	return params;
}

// Calculate depth-dependent threshold for a pair of samples
// Uses NGSCheckMate's predefined model approach
double get_depth_dependent_threshold(double depth1, double depth2, int family_mode)
{
	double min_depth = depth1 < depth2 ? depth1 : depth2;
	model_params_t params = get_predefined_model(min_depth, family_mode);
	
	// Calculate threshold as midpoint between matched and unmatched distributions
	// This is a simple approach; NGSCheckMate uses more sophisticated classification
	double threshold = (params.mean_matched + params.mean_unmatched) / 2.0;
	
	return threshold;
}

// Classify matches based on threshold (with optional depth adjustment)
void classify_matches(correlation_data_t *data, double base_threshold, FILE *out_fp, int verbose, int family_mode, int use_predefined_model)
{
	int i, j;
	int n_matches = 0;
	
	if (use_predefined_model && data->has_depth_info) {
		fprintf(out_fp, "# Match classification using NGSCheckMate predefined model (%s mode)\n", 
		        family_mode ? "family" : "non-family");
		fprintf(out_fp, "Sample1\tSample2\tDepth1\tDepth2\tCorrelation\tThreshold\tStatus\n");
	} else if (data->has_depth_info) {
		fprintf(out_fp, "# Match classification with base threshold %.4f (depth-adjusted)\n", base_threshold);
		fprintf(out_fp, "Sample1\tSample2\tDepth1\tDepth2\tCorrelation\tThreshold\tStatus\n");
	} else {
		fprintf(out_fp, "# Match classification with correlation threshold >= %.4f\n", base_threshold);
		fprintf(out_fp, "Sample1\tSample2\tCorrelation\tStatus\n");
	}
	
	for (i = 0; i < data->n_samples; ++i) {
		for (j = i + 1; j < data->n_samples; ++j) {
			double r = data->corr_matrix[i][j];
			double threshold = base_threshold;
			const char *status;
			
			// Use predefined model or simple depth-dependent threshold
			if (data->has_depth_info && use_predefined_model) {
				threshold = get_depth_dependent_threshold(
					data->samples[i].avg_depth,
					data->samples[j].avg_depth,
					family_mode);
			} else if (data->has_depth_info) {
				// Simple depth adjustment (legacy mode)
				double min_depth = data->samples[i].avg_depth < data->samples[j].avg_depth ? 
				                   data->samples[i].avg_depth : data->samples[j].avg_depth;
				if (min_depth < 5) {
					threshold = base_threshold - 0.05;
				} else if (min_depth > 15) {
					threshold = base_threshold + 0.02;
				}
			}
			
			if (r >= threshold) {
				status = "MATCHED";
				n_matches++;
			} else {
				status = "UNMATCHED";
			}
			
			// Only output matches by default, or all if verbose
			if (r >= threshold || verbose) {
				if (data->has_depth_info) {
					fprintf(out_fp, "%s\t%s\t%.2f\t%.2f\t%.6f\t%.4f\t%s\n",
					        data->samples[i].name,
					        data->samples[j].name,
					        data->samples[i].avg_depth,
					        data->samples[j].avg_depth,
					        r, threshold, status);
				} else {
					fprintf(out_fp, "%s\t%s\t%.6f\t%s\n",
					        data->samples[i].name,
					        data->samples[j].name,
					        r, status);
				}
			}
		}
	}
	
	if (use_predefined_model) {
		fprintf(stderr, "[M::%s] Found %d matched pairs using predefined model\n", 
		        __func__, n_matches);
	} else {
		fprintf(stderr, "[M::%s] Found %d matched pairs (threshold >= %.4f)\n", 
		        __func__, n_matches, base_threshold);
	}
}

// Calculate optimal threshold from training data
double calculate_optimal_threshold(correlation_data_t *data, const char *matched_list_fn)
{
	FILE *fp;
	char line[MAX_LINE];
	int n_matched_pairs = 0;
	double sum_matched = 0.0;
	double sum_unmatched = 0.0;
	int n_unmatched_pairs = 0;
	int i, j;
	char **matched_pairs;
	int n_pairs = 0, max_pairs = 100;
	
	// Load list of known matched pairs
	fp = fopen(matched_list_fn, "r");
	if (!fp) {
		fprintf(stderr, "Warning: could not open matched pairs file, using default threshold\n");
		return 0.95;
	}
	
	matched_pairs = (char**)malloc(max_pairs * sizeof(char*));
	while (fgets(line, MAX_LINE, fp) && n_pairs < max_pairs) {
		if (line[0] == '#') continue;
		
		// Remove newline
		line[strcspn(line, "\n")] = 0;
		matched_pairs[n_pairs] = strdup(line);
		n_pairs++;
	}
	fclose(fp);
	
	// Calculate average correlation for matched and unmatched pairs
	for (i = 0; i < data->n_samples; ++i) {
		for (j = i + 1; j < data->n_samples; ++j) {
			char pair_str[512];
			int is_matched = 0;
			int k;
			
			snprintf(pair_str, sizeof(pair_str), "%s\t%s", 
			         data->samples[i].name, data->samples[j].name);
			
			// Check if this is a known matched pair
			for (k = 0; k < n_pairs; ++k) {
				if (strstr(matched_pairs[k], data->samples[i].name) &&
				    strstr(matched_pairs[k], data->samples[j].name)) {
					is_matched = 1;
					break;
				}
			}
			
			if (is_matched) {
				sum_matched += data->corr_matrix[i][j];
				n_matched_pairs++;
			} else {
				sum_unmatched += data->corr_matrix[i][j];
				n_unmatched_pairs++;
			}
		}
	}
	
	// Free matched pairs
	for (i = 0; i < n_pairs; ++i) {
		free(matched_pairs[i]);
	}
	free(matched_pairs);
	
	if (n_matched_pairs == 0) {
		fprintf(stderr, "Warning: no matched pairs found in training data\n");
		return 0.95;
	}
	
	double avg_matched = sum_matched / n_matched_pairs;
	double avg_unmatched = n_unmatched_pairs > 0 ? sum_unmatched / n_unmatched_pairs : 0.0;
	
	// Set threshold midway between matched and unmatched averages
	double threshold = (avg_matched + avg_unmatched) / 2.0;
	
	fprintf(stderr, "[M::%s] Training statistics:\n", __func__);
	fprintf(stderr, "[M::%s]   Matched pairs: %d, avg correlation: %.4f\n", 
	        __func__, n_matched_pairs, avg_matched);
	fprintf(stderr, "[M::%s]   Unmatched pairs: %d, avg correlation: %.4f\n", 
	        __func__, n_unmatched_pairs, avg_unmatched);
	fprintf(stderr, "[M::%s]   Calculated threshold: %.4f\n", __func__, threshold);
	
	return threshold;
}

int main(int argc, char *argv[])
{
	int c, verbose = 0, family_mode = 0, use_predefined_model = 0;
	char *corr_fn = NULL, *out_fn = NULL, *train_fn = NULL;
	double threshold = 0.95; // Default threshold like NGSCheckMate
	correlation_data_t *data;
	FILE *out_fp;
	ketopt_t o = KETOPT_INIT;
	int vaf_start_idx = -1;
	
	while ((c = ketopt(&o, argc, argv, 1, "c:o:t:T:vFP", 0)) >= 0) {
		if (c == 'c') corr_fn = o.arg;
		else if (c == 'o') out_fn = o.arg;
		else if (c == 't') threshold = atof(o.arg);
		else if (c == 'T') train_fn = o.arg;
		else if (c == 'v') verbose = 1;
		else if (c == 'F') family_mode = 1;
		else if (c == 'P') use_predefined_model = 1;
	}
	
	vaf_start_idx = o.ind;
	
	if (!corr_fn || !out_fn) {
		fprintf(stderr, "Usage: match-classifier -c <correlation.corr> -o <matches.txt> [options] [VAF files...]\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -c FILE    input correlation matrix file\n");
		fprintf(stderr, "  -o FILE    output matches file\n");
		fprintf(stderr, "  -t FLOAT   correlation threshold for matching [%.2f]\n", threshold);
		fprintf(stderr, "  -T FILE    training file with known matched pairs (auto-calculate threshold)\n");
		fprintf(stderr, "  -P         use NGSCheckMate predefined model (requires VAF files for depth)\n");
		fprintf(stderr, "  -F         family mode (for related samples, used with -P)\n");
		fprintf(stderr, "  -v         verbose mode (output all pairs, not just matches)\n");
		fprintf(stderr, "\nDefault thresholds (NGSCheckMate-inspired):\n");
		fprintf(stderr, "  r >= 0.95  : Matched (same individual or technical replicates)\n");
		fprintf(stderr, "  r >= 0.80  : Possibly related (siblings, parent-child)\n");
		fprintf(stderr, "  r <  0.80  : Unrelated\n");
		fprintf(stderr, "\nDepth-dependent thresholds (with -P):\n");
		fprintf(stderr, "  depth > 10 : High confidence matching\n");
		fprintf(stderr, "  depth 5-10 : Medium confidence\n");
		fprintf(stderr, "  depth 2-5  : Lower confidence\n");
		fprintf(stderr, "  depth < 2  : Very low confidence\n");
		return 1;
	}
	
	fprintf(stderr, "[M::%s] Loading correlation matrix from %s...\n", __func__, corr_fn);
	data = load_correlation_matrix(corr_fn);
	if (!data) {
		fprintf(stderr, "Error: failed to load correlation matrix\n");
		return 1;
	}
	fprintf(stderr, "[M::%s] Loaded %d samples\n", __func__, data->n_samples);
	
	// Load depth information from VAF files if using predefined model
	if (use_predefined_model && vaf_start_idx < argc) {
		fprintf(stderr, "[M::%s] Loading depth information from VAF files...\n", __func__);
		if (!load_depth_info(data, argc, argv, vaf_start_idx)) {
			fprintf(stderr, "Warning: -P flag used but no depth info loaded, using simple threshold\n");
			use_predefined_model = 0;
		} else {
			fprintf(stderr, "[M::%s] Using NGSCheckMate predefined model (%s mode)\n", 
			        __func__, family_mode ? "family" : "non-family");
		}
	}
	
	// Calculate optimal threshold from training data if provided
	if (train_fn) {
		fprintf(stderr, "[M::%s] Calculating threshold from training data...\n", __func__);
		threshold = calculate_optimal_threshold(data, train_fn);
	}
	
	out_fp = fopen(out_fn, "w");
	if (!out_fp) {
		fprintf(stderr, "Error: failed to open output file\n");
		free_correlation_data(data);
		return 1;
	}
	
	if (use_predefined_model) {
		fprintf(stderr, "[M::%s] Classifying matches using predefined model...\n", __func__);
	} else {
		fprintf(stderr, "[M::%s] Classifying matches with threshold %.4f...\n", __func__, threshold);
	}
	classify_matches(data, threshold, out_fp, verbose, family_mode, use_predefined_model);
	
	fclose(out_fp);
	free_correlation_data(data);
	
	fprintf(stderr, "[M::%s] Results written to %s\n", __func__, out_fn);
	
	return 0;
}
