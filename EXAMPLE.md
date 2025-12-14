# SNP K-mer Analysis Tools - Usage Example

This document provides a step-by-step example of using the SNP k-mer analysis pipeline.

## Overview

The pipeline consists of four programs:
1. `snp-pattern-gen` - Extracts unique k-mers from SNP positions
2. `vaf-counter` - Counts k-mer occurrences in sequencing data
3. `correlation-matrix` - Computes sample correlations
4. `match-classifier` - Classifies matched samples with depth-dependent thresholds

## Step 1: Prepare Input Files

### BED File Format

Create a BED file containing SNP positions. Format: `chr start end rsID ref alt`

Example `snps.bed`:
```
chr17	46549406	46549407	rs201103889	A	C
chr1	152308305	152308306	rs2184953	T	C
chr19	4945904	4945905	rs2250981	C	T
chr3	12816528	12816529	rs3732675	C	G
chr19	41397314	41397315	rs10853751	C	T
```

### Reference Genome

You'll need a reference genome in FASTA format (e.g., hg38.fa, GRCh37.fa).

## Step 2: Generate SNP K-mer Patterns

Extract unique k-mers around SNP positions:

```bash
./snp-pattern-gen -k 21 -b snps.bed -f hg38.fa -o patterns.txt
```

Options:
- `-k 21`: Use 21-mer (10 bases on each side of SNP)
- `-b snps.bed`: Input BED file with SNPs
- `-f hg38.fa`: Reference genome FASTA file
- `-o patterns.txt`: Output pattern file

The program will:
1. Load the reference genome into memory
2. Count all k-mers in the genome
3. For each SNP, extract the reference and alternative k-mers
4. Check if both k-mers are unique in the genome
5. Output only SNPs with unique k-mer pairs

Output example (`patterns.txt`):
```
chr17	46549406	46549407	rs201103889	A	C	AAATTTGGGCCCAAATTTCCC	AAATTTGGGCCCAAATTTGCC
chr1	152308305	152308306	rs2184953	T	C	GCTAGCTAGCTAGCTAGCTA	GCTAGCTAGCAAGCTAGCTA
```

## Step 3: Count K-mers in Sequencing Data

You have three options for generating VAF files, depending on your input data:

### Option A: From FASTQ files (unaligned reads)

For each sample, count occurrences of reference and alternative k-mers:

```bash
# Sample 1
./vaf-counter -k 21 -t 4 -p patterns.txt -o sample1.vaf \
    sample1_R1.fastq.gz sample1_R2.fastq.gz

# Sample 2
./vaf-counter -k 21 -t 4 -p patterns.txt -o sample2.vaf \
    sample2_R1.fastq.gz sample2_R2.fastq.gz

# Sample 3
./vaf-counter -k 21 -t 4 -p patterns.txt -o sample3.vaf \
    sample3_R1.fastq.gz sample3_R2.fastq.gz
```

Options:
- `-k 21`: K-mer length (must match pattern generation)
- `-p patterns.txt`: Pattern file from step 2
- `-o sample1.vaf`: Output VAF file
- `-t 4`: Number of threads (default: 4, for faster processing)
- `-b INT`: Block size for batching (default: 10000000)
- FASTQ files can be gzipped or uncompressed

**Performance Note:** vaf-counter uses multi-threading to parallelize k-mer counting, similar to kc-c4.c. The `-t` option controls the number of worker threads used for k-mer lookup operations.

### Option B: From BAM files (aligned reads)

For aligned sequencing data, use bam-vaf-counter with htslib:

```bash
# Index BAM files first (if not already indexed) for optimal performance
samtools index sample1.bam
samtools index sample2.bam
samtools index sample3.bam

# Sample 1
./bam-vaf-counter -k 21 -t 4 -p patterns.txt -o sample1.vaf sample1.bam

# Sample 2
./bam-vaf-counter -k 21 -t 4 -p patterns.txt -o sample2.vaf sample2.bam

# Sample 3
./bam-vaf-counter -k 21 -t 4 -p patterns.txt -o sample3.vaf sample3.bam
```

Options:
- `-k 21`: K-mer length (must match pattern generation)
- `-p patterns.txt`: Pattern file from step 2
- `-o sample1.vaf`: Output VAF file
- `-t 4`: Number of threads (default: 4)

**Note:** This tool reads BAM/SAM/CRAM files using htslib and extracts k-mers from aligned reads. The output format is identical to vaf-counter.

**Performance:** bam-vaf-counter is optimized to use indexed BAM access (requires .bai file) to fetch only reads from SNP regions instead of processing all reads. This provides significant speedup (potentially 100x-1000x faster) when working with large BAM files. Make sure your BAM files are indexed with `samtools index` before running bam-vaf-counter for best performance.

### Option C: From VCF files (variant calls)

If you already have called variants with genotype information, use vcf-vaf-counter:

```bash
# Sample 1 (first sample in VCF)
./vcf-vaf-counter -p patterns.txt -v sample1.vcf.gz -o sample1.vaf -s 0 -d 5

# Sample 2 (second sample in VCF, if multi-sample)
./vcf-vaf-counter -p patterns.txt -v samples.vcf.gz -o sample2.vaf -s 1 -d 5

# Sample 3 (third sample in VCF)
./vcf-vaf-counter -p patterns.txt -v samples.vcf.gz -o sample3.vaf -s 2 -d 5
```

Options:
- `-p patterns.txt`: Pattern file from step 2
- `-v sample1.vcf.gz`: Input VCF/BCF file (can be compressed)
- `-o sample1.vaf`: Output VAF file
- `-s 0`: Sample index (0-based) for multi-sample VCF [default: 0]
- `-d 5`: Minimum depth filter [default: 1]

**Note:** This tool reads VCF files using htslib and extracts genotype (GT) and allele depth (AD or DP) information. If AD is available, it uses actual allele counts; otherwise, it estimates from DP and genotype.

### Output Format

All three methods produce the same VAF file format:

```
# Average depth: 25.50
CHR	POS	RSID	REF	ALT	REF_COUNT	ALT_COUNT	TOTAL_COUNT	VAF
chr17	46549406	rs201103889	A	C	45	5	50	0.1000
chr1	152308305	rs2184953	T	C	48	2	50	0.0400
```

## Step 4: Compute Sample Correlation Matrix

Compare all samples to identify relationships. **Use different modes for different scenarios:**

### For Matched Samples (Same Individual, Technical Replicates)
```bash
./correlation-matrix -M matched -o correlation.corr -t \
    sample1.vaf sample2_replicate.vaf
```
Uses **stricter** thresholds: depth ≥5, minimum 10 SNPs (higher quality requirement)

### For Unmatched/Related Samples
```bash
./correlation-matrix -M unmatched -o correlation.corr -t \
    sample1.vaf sample2.vaf sample3.vaf
```
Uses **lenient** thresholds: depth ≥1, minimum 20 SNPs (more SNPs needed for confidence)

### For High-Confidence Analysis
```bash
./correlation-matrix -M strict -o correlation.corr \
    sample1.vaf sample2.vaf
```
Uses **strictest** thresholds: depth ≥10, minimum 30 SNPs

### Custom Configuration
```bash
./correlation-matrix -d 3 -m 15 -o correlation.corr sample1.vaf sample2.vaf
```

Options:
- `-M MODE`: Preset mode (`matched`, `unmatched`, or `strict`)
- `-d INT`: Minimum depth per SNP (overrides mode preset)
- `-m INT`: Minimum SNPs with sufficient depth (overrides mode preset)
- `-o FILE`: Output correlation matrix file
- `-t`: Generate dendrogram/tree file (optional)

**Important:** Different scenarios require **different depth cutoffs**:
- **Matched samples**: Higher depth cutoff ensures confident matching (avoids false negatives)
- **Unmatched samples**: Lower depth cutoff but more SNPs for statistical confidence
- The mode presets follow NGSCheckMate's approach to depth-dependent correlation

Output example (`correlation.corr`):
```
Sample	sample1	sample2	sample3
sample1	1.000000	0.950000	0.450000
sample2	0.950000	1.000000	0.480000
sample3	0.450000	0.480000	1.000000
```

If `-t` flag is used, a tree file is also generated (`correlation.tree`):
```
# Simple dendrogram (UPGMA-like clustering)
# Format: (Sample1:distance, Sample2:distance)
Cluster: sample1 (0.0250) <-> sample2 (0.0250)
Cluster: sample1 (0.2750) <-> sample3 (0.2750)
```

## Interpretation

### VAF Values
- **VAF = 0.0**: Homozygous reference (only reference allele present)
- **VAF = 0.5**: Heterozygous (equal amounts of reference and alternative alleles)
- **VAF = 1.0**: Homozygous alternative (only alternative allele present)

### Correlation Values
- **r = 1.0**: Identical samples (perfect correlation)
- **r > 0.95**: Likely the same individual or clonal replicates
- **r = 0.8-0.95**: Possibly related individuals (siblings, parent-child)
- **r < 0.5**: Unrelated samples
- **r = 0.0**: Insufficient data (fewer than minimum SNPs with adequate depth)

**Note on Depth Dependency:** Correlation values are depth-dependent. Low-coverage samples may return 0.0 correlation even if they are from the same individual, simply because there aren't enough SNPs with sufficient depth to calculate a reliable correlation. Always check that samples have adequate depth before interpreting low correlation values.

## Step 5: Classify Matched Samples

Use the match classifier to automatically identify which samples are matched (same individual):

### Simple Threshold Mode
```bash
./match-classifier -c correlation.corr -o matches.txt -t 0.95
```

Output (`matches.txt`):
```
# Match classification with correlation threshold >= 0.9500
Sample1	Sample2	Correlation	Status
sample1	sample2	0.980000	MATCHED
```

### NGSCheckMate Predefined Model (Recommended)

For more accurate classification using depth-dependent thresholds:

```bash
# Non-family mode (for unrelated samples)
./match-classifier -c correlation.corr -o matches.txt -P \
    sample1.vaf sample2.vaf sample3.vaf

# Family mode (for related samples)
./match-classifier -c correlation.corr -o matches_family.txt -P -F \
    sample1.vaf sample2.vaf sample3.vaf
```

Output with depth-dependent thresholds (`matches.txt`):
```
# Match classification using NGSCheckMate predefined model (non-family mode)
Sample1	Sample2	Depth1	Depth2	Correlation	Threshold	Status
sample1	sample2	15.50	12.00	0.980000	0.5925	MATCHED
sample1	sample3	15.50	3.00	0.450000	0.4448	MATCHED
```

**How it works:**
- The classifier uses **depth-stratified thresholds** from NGSCheckMate
- Higher depth samples require higher correlation to be matched (more confident)
- Lower depth samples use lower thresholds (less confident, but still useful)
- Family mode adjusts thresholds for related samples vs unrelated

**Depth-dependent thresholds:**
- **depth > 10**: High confidence (threshold ~0.87 non-family, ~0.78 family)
- **depth 5-10**: Medium confidence
- **depth 2-5**: Lower confidence  
- **depth < 2**: Very low confidence (threshold ~0.53 non-family, ~0.58 family)

### Applications
1. **Sample Mix-up Detection**: Automatically identify mislabeled samples
2. **Contamination Detection**: Identify samples with unexpected matches
3. **Relatedness Estimation**: Determine family relationships between samples
4. **Quality Control**: Verify technical replicates are correctly matched

## Performance Tips

1. **K-mer Length**: 
   - Longer k-mers (21-31) provide better specificity
   - Must be odd to center on the SNP position

2. **SNP Selection**:
   - Use common SNPs (MAF > 0.1) for better coverage
   - Select SNPs with unique flanking sequences
   - Aim for 10,000-50,000 SNPs for good coverage

3. **Memory Usage**:
   - `snp-pattern-gen` loads entire genome into memory
   - For human genome (~3GB), expect 5-10GB RAM usage
   - `vaf-counter` and `correlation-matrix` use minimal memory

4. **Speed**:
   - Pattern generation is I/O bound (genome loading)
   - K-mer counting uses multi-threading (use `-t` to specify threads, default: 4)
   - More threads significantly speed up k-mer counting on large FASTQ files
   - Correlation computation is very fast (< 1 second)

## Common Issues

### No unique k-mer pairs found
- Try increasing k-mer length
- Ensure reference genome matches SNP coordinates
- Check that SNP positions are valid (not in repetitive regions)

### Low depth or coverage
- Increase sequencing depth
- Use more SNPs in the pattern file
- Check FASTQ file quality

### Unexpected correlations
- Verify sample labels
- Check for cross-contamination
- Ensure reference genome matches sequencing data

## References

This pipeline is inspired by:
- NGSCheckMate: https://github.com/parklab/NGSCheckMate
- K-mer counting techniques from this repository

## Contact

For issues or questions, please refer to the main repository README.
