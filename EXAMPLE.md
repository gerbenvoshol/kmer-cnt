# SNP K-mer Analysis Tools - Usage Example

This document provides a step-by-step example of using the SNP k-mer analysis pipeline.

## Overview

The pipeline consists of three programs:
1. `snp-pattern-gen` - Extracts unique k-mers from SNP positions
2. `vaf-counter` - Counts k-mer occurrences in sequencing data
3. `correlation-matrix` - Computes sample correlations

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

For each sample, count occurrences of reference and alternative k-mers:

```bash
# Sample 1
./vaf-counter -k 21 -p patterns.txt -o sample1.vaf \
    sample1_R1.fastq.gz sample1_R2.fastq.gz

# Sample 2
./vaf-counter -k 21 -p patterns.txt -o sample2.vaf \
    sample2_R1.fastq.gz sample2_R2.fastq.gz

# Sample 3
./vaf-counter -k 21 -p patterns.txt -o sample3.vaf \
    sample3_R1.fastq.gz sample3_R2.fastq.gz
```

Options:
- `-k 21`: K-mer length (must match pattern generation)
- `-p patterns.txt`: Pattern file from step 2
- `-o sample1.vaf`: Output VAF file
- FASTQ files can be gzipped or uncompressed

Output example (`sample1.vaf`):
```
# Average depth: 25.50
CHR	POS	RSID	REF	ALT	REF_COUNT	ALT_COUNT	TOTAL_COUNT	VAF
chr17	46549406	rs201103889	A	C	45	5	50	0.1000
chr1	152308305	rs2184953	T	C	48	2	50	0.0400
```

## Step 4: Compute Sample Correlation Matrix

Compare all samples to identify relationships:

```bash
./correlation-matrix -o correlation.corr -t \
    sample1.vaf sample2.vaf sample3.vaf
```

Options:
- `-o correlation.corr`: Output correlation matrix file
- `-t`: Generate dendrogram/tree file (optional)
- `-m INT`: Minimum SNPs with depth ≥1 required for correlation (default: 20)

**Important:** The correlation calculation is **depth-aware**:
- Only SNPs with depth ≥1 in both samples are included
- Requires minimum 20 valid SNPs by default (like NGSCheckMate)
- Returns 0.0 correlation if insufficient SNPs have adequate depth
- Prevents spurious correlations from low-coverage samples

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

### Applications
1. **Sample Mix-up Detection**: Verify that labeled samples match expected correlations
2. **Contamination Detection**: Identify samples with unexpected VAF patterns
3. **Relatedness Estimation**: Determine family relationships between samples
4. **Quality Control**: Ensure technical replicates show high correlation

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
   - K-mer counting benefits from gzipped FASTQ files
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
