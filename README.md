## SNP K-mer Analysis Tools

This repository includes seven C programs for SNP-based sample correlation analysis, inspired by [NGSCheckMate](https://github.com/parklab/NGSCheckMate):

### Overview

The pipeline consists of seven programs:

1. **snp-pattern-gen** - Extracts unique k-mers from SNP positions
2. **vaf-counter** - Counts k-mer occurrences in FASTQ files and calculates variant allele frequencies
3. **ed-vaf-counter** - Alternative FASTQ k-mer counter using edlib for approximate matching
4. **bam-vaf-counter** - Counts k-mer occurrences in BAM files (uses htslib)
5. **vcf-vaf-counter** - Generates VAF files directly from VCF files (uses htslib)
6. **correlation-matrix** - Computes Pearson correlation between samples
7. **match-classifier** - Classifies matched samples using depth-dependent thresholds

### Usage

#### 1. Generate SNP k-mer patterns

```sh
./snp-pattern-gen -k 21 -b snps.bed -f reference.fa -o patterns.txt
```

**Input:**
- `-b` BED file with SNPs (format: chr, start, end, rsID, ref, alt)
- `-f` Reference genome in FASTA format
- `-k` K-mer length (must be odd, default: 21)

**Output:** Pattern file containing unique reference and alternative k-mers for each SNP

**Example BED file:**
```
chr17	46549406	46549407	rs201103889	A	C
chr1	152308305	152308306	rs2184953	T	C
chr19	4945904	4945905	rs2250981	C	T
```

#### 2. Count k-mers in FASTQ files

```sh
./vaf-counter -k 21 -t 4 -p patterns.txt -o sample1.vaf reads1.fq reads2.fq

# With verbose performance reporting
./vaf-counter -k 21 -t 4 -p patterns.txt -o sample1.vaf -v reads1.fq reads2.fq
```

**Input:**
- `-p` Pattern file from step 1
- `-k` K-mer length (must match pattern generation)
- `-t` Number of threads (default: 4, for faster multi-threaded processing)
- `-b` Block size for batching sequences (default: 10000000)
- `-v` Verbose mode: report detailed performance statistics
- One or more FASTQ files (can be gzipped)

**Output:** VAF file with variant allele frequencies and depth information

**Performance:** vaf-counter uses multi-threaded k-mer counting similar to kc-c4.c, with a 3-stage pipeline (read sequences, extract k-mers, lookup k-mers) for optimal performance. The tool includes SIMD optimizations (SSSE3/SSE4.1) for faster sequence encoding and prefetching for reduced cache misses. Use the `-v` flag to see detailed performance metrics including throughput (Mbases/sec), k-mer extraction rate, and SIMD optimization level.

#### 2alt. Alternative: Count k-mers using approximate matching (ed-vaf-counter)

For situations where approximate matching is desired, or when the pattern set is small:

```sh
./ed-vaf-counter -p patterns.txt -o sample1.vaf -e 0 reads1.fq reads2.fq
```

**Input:**
- `-p` Pattern file from step 1
- `-e` Maximum edit distance for approximate matching (default: 0 for exact match)
- One or more FASTQ files (can be gzipped)

**Output:** VAF file with variant allele frequencies and depth information

**Approach:** Unlike vaf-counter which extracts all k-mers from reads and looks them up in a hash table, ed-vaf-counter uses [edlib](https://github.com/Martinsos/edlib) to search for each pattern k-mer in the FASTQ reads using approximate string matching. This approach:
- Reverses the search direction (pattern→reads instead of reads→patterns)
- Allows approximate matching with configurable edit distance (-e parameter)
- Can be more efficient when the number of patterns is small relative to total k-mers
- Enables fuzzy matching to account for sequencing errors

**Use cases:**
- When you want to allow mismatches/indels in k-mer matching (set -e > 0)
- When you have a very small pattern set and large FASTQ files
- For exploratory analysis with approximate matching

#### 2b. Alternative: Count variants in BAM files

For aligned sequencing data (BAM/SAM/CRAM files):

```sh
./bam-vaf-counter -t 4 -p patterns.txt -o sample1.vaf sample1.bam
```

**Input:**
- `-p` Pattern file from step 1
- `-t` Number of threads (default: 4)
- One or more BAM/SAM/CRAM files

**Output:** VAF file with variant allele frequencies and depth information

**Note:** This program uses htslib to read BAM files and directly counts ref/alt bases at SNP positions without k-mer extraction, making it significantly faster than the k-mer-based approach.

**Performance Optimization:** bam-vaf-counter uses indexed BAM access (requires .bai index file) to fetch only reads overlapping SNP positions, dramatically improving performance compared to sequential processing of all reads. The program directly examines bases at SNP positions using CIGAR string parsing, eliminating the overhead of k-mer extraction. Regions are automatically merged to minimize redundant fetching. If no BAM index is available, the program falls back to sequential processing with a warning.

#### 2c. Alternative: Generate VAF from VCF files

For variant call data (VCF/BCF files):

```sh
./vcf-vaf-counter -p patterns.txt -v sample1.vcf.gz -o sample1.vaf -s 0 -d 5
```

**Input:**
- `-p` Pattern file from step 1
- `-v` VCF or BCF file (can be compressed)
- `-s` Sample index (0-based) for multi-sample VCF [default: 0]
- `-d` Minimum depth filter [default: 1]

**Output:** VAF file with variant allele frequencies and depth information

**Note:** This program uses htslib to read VCF files and extracts genotype (GT) and allele depth (AD) or depth (DP) information. If AD is available, it uses actual allele counts; otherwise, it estimates from DP and genotype. This is useful when you already have called variants and want to skip the k-mer counting step.

#### 3. Compute correlation matrix

```sh
# Using preset mode for matched samples (same individual)
./correlation-matrix -M matched -o correlation.corr -t sample1.vaf sample2.vaf

# Using preset mode for unmatched/related samples
./correlation-matrix -M unmatched -o correlation.corr -t sample1.vaf sample2.vaf sample3.vaf

# Manual configuration
./correlation-matrix -d 5 -m 15 -o correlation.corr sample1.vaf sample2.vaf
```

**Input:**
- Multiple VAF files from step 2
- `-M MODE` Preset modes: `matched` (depth≥5, SNPs≥10), `unmatched` (depth≥1, SNPs≥20), `strict` (depth≥10, SNPs≥30)
- `-t` Optional flag to generate dendrogram/tree
- `-d INT` Minimum depth per SNP (default: 1)
- `-m INT` Minimum SNPs with sufficient depth required (default: 20)

**Output:** 
- Correlation matrix file (depth-aware Pearson correlation)
- Optional tree file (UPGMA clustering)

**Note:** The correlation calculation is depth-aware with different thresholds for different scenarios:
- **Matched samples** (same individual, replicates): Use higher depth cutoff (≥5) for confident matches
- **Unmatched samples** (related/unrelated): Use lower depth cutoff (≥1) but require more SNPs (≥20)
- Manual `-m` and `-d` flags override preset modes

#### 4. Classify matched samples

```sh
# Simple threshold-based classification
./match-classifier -c correlation.corr -o matches.txt -t 0.95

# NGSCheckMate predefined model with depth-dependent thresholds
./match-classifier -c correlation.corr -o matches.txt -P sample1.vaf sample2.vaf sample3.vaf

# Family mode (for related samples)
./match-classifier -c correlation.corr -o matches.txt -P -F sample1.vaf sample2.vaf
```

**Input:**
- `-c` Correlation matrix from step 3
- `-P` Use NGSCheckMate predefined model (requires VAF files for depth info)
- `-F` Family mode for related samples (use with `-P`)
- `-t` Manual threshold (default: 0.95)

**Output:** List of matched sample pairs with status and depth-dependent thresholds

**Note:** The classifier uses depth-dependent thresholds following NGSCheckMate:
- **depth > 10**: High confidence matching (stricter thresholds)
- **depth 5-10**: Medium confidence
- **depth 2-5**: Lower confidence  
- **depth < 2**: Very low confidence (more lenient thresholds)

### Quick Start

```sh
# Build the tools (includes htslib)
make snp-pattern-gen vaf-counter ed-vaf-counter bam-vaf-counter vcf-vaf-counter correlation-matrix match-classifier

# Run the complete pipeline with FASTQ files
./snp-pattern-gen -k 21 -b dbSNP.bed -f hg38.fa -o patterns.txt
./vaf-counter -k 21 -p patterns.txt -o sample1.vaf sample1_R1.fq.gz sample1_R2.fq.gz
./vaf-counter -k 21 -p patterns.txt -o sample2.vaf sample2_R1.fq.gz sample2_R2.fq.gz
./correlation-matrix -M matched -o correlation.corr sample1.vaf sample2.vaf
./match-classifier -c correlation.corr -o matches.txt -P sample1.vaf sample2.vaf

# Or use ed-vaf-counter for approximate matching (allows mismatches)
./ed-vaf-counter -p patterns.txt -e 1 -o sample1.vaf sample1_R1.fq.gz sample1_R2.fq.gz
./ed-vaf-counter -p patterns.txt -e 1 -o sample2.vaf sample2_R1.fq.gz sample2_R2.fq.gz

# Or use BAM files instead (position-based, no k-mer needed)
./bam-vaf-counter -p patterns.txt -o sample1.vaf sample1.bam
./bam-vaf-counter -p patterns.txt -o sample2.vaf sample2.bam

# Or use VCF files directly
./vcf-vaf-counter -p patterns.txt -v sample1.vcf.gz -o sample1.vaf -s 0
./vcf-vaf-counter -p patterns.txt -v sample2.vcf.gz -o sample2.vaf -s 0
```

### Applications

- Sample identity verification and mix-up detection
- Cross-contamination detection
- Sample relatedness estimation
- Quality control in sequencing projects

---

## K-mer Counting Tools - Getting Started

```sh
git clone https://github.com/lh3/kmer-cnt
cd kmer-cnt
make  # C++11 required to compile the two C++ implementations
wget https://github.com/lh3/kmer-cnt/releases/download/v0.1/M_abscessus_HiSeq_10M.fa.gz
./yak-count M_abscessus_HiSeq_10M.fa.gz > kc-c4.out
```

## Introduction

K-mer counting is the foundation of many mappers, assemblers and miscellaneous
tools (e.g. genotypers, metagenomics profilers, etc). It is one of the most
important classes of algorithms in Bioinformatics. Here we will implement basic
k-mer counting algorithms but with advanced engineering tricks. We will see how
far better engineering can go.

In this repo, each `{kc,yak}-*.*` file implements a standalone k-mer counter.
As to other files: ketopt.h is a command line option parser; khashl.h is a
generic hash table library in C; kseq.h is a fasta/fastq parser; kthread.{h,c}
provides two multi-threading models; robin\_hood.h is a C++11 hash table
library.

## Results

We provide eight k-mer counters, which are detailed below the result table. All
implementations count canonical k-mers, the lexicographically *smaller* k-mer
between the k-mers on the two DNA strands.

The following table shows the timing and peak memory of different
implementations for counting 31-mers from 2.5 million pairs of 100bp reads
sampled from the HiSeq *M. abscessus* 6G-0125-R dateset in [GAGE-B][gage-b].
They were run on a Linux server equipped with two EPYC 7301 CPUs and 512GB RAM.

|Implementation                 |Limitation          |Elapsed time (s)|CPU time (s)|Peak RAM (GB)|
|:------------------------------|:-------------------|---------------:|-----------:|------------:|
|[kc-py1](kc-py1.py) + Python3.7|                    |           499.6|       499.5|         8.15|
|[kc-py1](kc-py1.py) + Pypy7.3  |                    |          1220.8|      1220.8|        12.21|
|[kc-cpp1](kc-cpp1.cpp)         |                    |           528.0|       527.9|         8.27|
|[kc-cpp2](kc-cpp2.cpp)         |                    |           319.6|       319.6|         6.90|
|[kc-c1](kc-c1.c)               |<=32-mer            |            39.3|        38.3|         1.52|
|[kc-c2](kc-c2.c)               |<=32-mer; <1024 count|           38.7|        37.9|         1.05|
|[kc-c3](kc-c3.c)               |<=32-mer; <1024 count|           34.1|        38.7|         1.15|
|[kc-c4](kc-c4.c) (2+4 threads) |<=32-mer; <1024 count|            7.5|        35.1|         1.27|
|[yak-count](yak-count.c) (2+4; >=2 count)|<=32-mer; <1024 count| 14.6|        54.8|         0.47|
|[jellyfish2][jf] (16 threads)  |                    |            10.8|       163.9|         0.82|
|[KMC3][KMC] (16 thr; in-mem)   |                    |             9.2|        36.2|         5.02|

## Discussions

* [kc-py1.py](kc-py1.py) is a basic Python3 implementation. It uses string
  translate for fast complementary. Interestingly, pypy is much slower than
  python3. Perhaps the official python3 comes with a better hash table
  implementation. Just a guess. I often recommend pypy over python. I need to
  be more careful about this recommendation in future.

* [kc-cpp1.cpp](kc-cpp1.cpp) implements a basic counter in C++11 using STL's
  [unordered\_map][unordermap]. It is slower than python3. This is partly
  because STL's hash table implementation is very inefficient. C++ does not
  necessarily lead to a fast implementation.

* [kc-cpp2.cpp](kc-cpp2.cpp) replaces `std::unordered_map` with Martin Ankerl's
  [robin\_hood][rhhash] hash table library, which is [among the
  fastest][rhbench] hash table implementations. It is now faster than
  kc-py1.py, though the performance gap is small.

* [kc-c1.c](kc-c1.c) packs k-mers no longer than 32bp into 64-bit integers.
  This dramatically improves speed and reduces the peak memory. Most practical
  k-mer counters employs bit packing. Excluding library files, this counter has
  less than 100 coding lines, not much more complex than the C++ or the python
  implementations.

* [kc-c2.c](kc-c2.c) uses an ensemble of hash tables to save 8 bits for
  counter. This reduces the peak memory. The key advantage of using multiple
  hash tables is to implement multithreading. See below.

* [kc-c3.c](kc-c3.c) puts file reading and parsing into a separate thread. The
  performance improvement is minor here, but it sets the stage for the next
  multi-threaded implementation.

* [kc-c4.c](kc-c4.c) is the fastest counter in this series. Due to the use of
  an ensembl of hash tables in kc-c2, we can parallelize the insertion of a
  batch of k-mers. It is much faster than the previous versions. Notably, kc-c4
  also uses less CPU time. This is probably because batching helps data
  locality.

* [yak-count.c](yak-count.c) is adapted from [yak][yak] and uses the same kc-c4
  algorithm. Similar to [BFCounter][BFCnt], it optionally adds a bloom filter
  to filter out most singleton k-mers (k-mers occurring only once in the
  input). Yak needs to update the bloom filter, read the input twice and count
  twice. It is slower but uses less memory. Yak-count is the most complex
  example in this repo, but it is still short. Its code is also better
  organized. Command line: `-b30` (bloom filter with 1 billion bits).

* [jellyfish2][jf] is probably the fastest in-memory k-mer counter to date. It
  uses less memory and more flexible than kc-c4, but it is slower and much more
  complex. Command line: `count -m 31 -C -s 100000000 -o /dev/null -t 16`.

* [KMC3][KMC] is one of the fastest k-mer counters. It uses minimizers and
  relies on sorting. KMC3 is run in the in-memory mode here. The disk mode is
  as fast. KMC3 is optimized for counting much larger datasets. Although it
  uses more RAM here, it generally uses less RAM than jellyfish2 and other
  in-memory counters given high-coverage human data. Command line: `-k31 -t16
  -r -fa`.

## Conclusions

The k-mer counters here are fairly basic implementations only using generic
hash tables. Nonetheless, we show better engineering can carry the basic idea a
long way. If you want to implement your own k-mer counter,
[yak-count.c](yak-count.c) could be a good starting point. It is fast and
relatively simple. By the way, if you have an efficient and simple k-mer
counter (implemented in a few files), please let me know. I will be happy to add it to the table.

[jf]: http://www.genome.umd.edu/jellyfish.html
[unordermap]: http://www.cplusplus.com/reference/unordered_map/unordered_map/
[rhhash]: https://github.com/martinus/robin-hood-hashing
[rhbench]: https://martin.ankerl.com/2019/04/01/hashmap-benchmarks-01-overview/
[gage-b]: https://ccb.jhu.edu/gage_b/datasets/index.html
[yak]: https://github.com/lh3/yak
[BFCnt]: https://github.com/pmelsted/BFCounter
[KMC]: https://github.com/refresh-bio/KMC
