# back-reference-qc

A bioinformatics pipeline that detects unreliable reads with low QV. Quality values are cross-referenced from its own kmerized e.g. Illumina library or your choice of complementary technology (Illumina or long reads) for the sample in query. It's also possible to filter out the quality reads but otherwise not required. By default, this pipeine will produce kernel density estimation plots before and after defined z-score parameter to reflect the distribution pattern.

## Getting started
1. Clone the repo
2. Install Snakemake version >= 7.16.0
3. Fill in the config and manifests.
   
## Anlaysis
 * Begin with a dry-run
 ```
 ./runlocal 10 -np
 ```
 * If dry-run looks good, proceed with:
 ```
 ./runlocal 10
 ```

## Analysis Output
```commandline
   results/
   ├── overwrite_target_lists
   │   └── unique_sample_name
   ├── plots
   │   └── unique_sample_name
   ├── read_qv
   │   └── unique_sample_name
   └── reads_filtered
       └── unique_sample_name
           ├── fastq
           ├── filtered_out
           │   ├── fasta
           │   └── kraken2
           └── log
```

```commandline
$ zcat results/reads_filtered/unique_sample_name/filtered_out/kraken2/summary.tsv.gz
sample  reads_mapped    taxonomy        cell_name
unique_sample_name        54      Mycoplasmataceae (taxid 2092)   m84046_230412_163336_s3.hifi_reads
total   54  N/A     54037280
```
The row total pertains to the sum of all reads mapped to the taxa (e.g. 54), and 54037280 represents the sum of all reads in cells (e.g. m84046_230412_163336_s3.hifi_reads)

```commandline
./results/plots/unique_sample_name
├── kde-after_filter.png
└── kde-before_filter.png
```
The plots showcases the kernel density estimate distribution for qv values before and after filtering. The expectation is that the reads are normally distributed.

```commandline
./results/reads_filtered/unique_sample_name/log/
└── m84046_230412_163336_s3.hifi_reads-extract_undesirable_reads.log

$ cat
comparison_type: self for sample: unique_sample_name
z_filter: -2 for sample: unique_sample_name
QV-99	99.0	m84046_230412_163336_s3.hifi_reads	0.00
QV-0	0.0	m84046_230412_163336_s3.hifi_reads	0.00
QV-low	19.02	m84046_230412_163336_s3.hifi_reads	0.00
```

When qc analysis is performed only if new_fastq=True is set in the config file, a preparation file for overwriting is generated from the information of the original fastq file location and the newly created fastq file.

```commandline
zcat results/overwrite_target_lists/unique_sample_name.overwrite_target_lists.tab.gz
SAMPLE   CELL   ORIGINAL_PATH   ORIGINAL_SIZE   ORIGINAL_MD5   CLEANED_PATH   CLEANED_SIZE   CLEANED_MD5
unique_sample_name   m84046_230412_163336_s3.hifi_reads   fastq/m84046_230412_163336_s3.hifi_reads.fastq.gz   16255055450   ca00950d807af12a5ee0e7ca6229cc44   results/reads_filtered/unique_sample_name/fastq/m84046_230412_163336_s3.hifi_reads.fastq_target-reads_with_reference_help-subset.fastq.gz   16539432159   932d1bf2f2fc65f01fc76640a2321308
```

If comparison_type and z_filter are throwaway headers. We can read the table where lines start with Q. I describe the columns as such:
  - QV thresholds
  - Median value
  - Cell/movie name
  - Proportion of reads that are less than defined z score/filter.

QV thresholds:
  * QV-99 are perfect k-mer matches
  * QV-0 do not match k-mers whatsoever (definitely non-human reads)
  * QV-low are reads beyond the define z_filter score (e.g. -2). 19.02 is the median value of reads that fall beyond a z-score of -2. If we round up to 20, it equates to an error rate of 1 in 100 bases called or 0.01% probability that the base is incorrect.

## Overwrite(optional)
 * Begin with a dry-run
 ```
 ./runlocal overwrite 10 -np
 ```
 * If dry-run looks good, proceed with:
 ```
 ./runlocal overwrite 10
 ```

## Overwrite Output
```commandline
   results/
   └── overwrite_records
       └── unique_sample_name
```

```commandline
zcat results/overwrite_records/unique_sample_name.overwrite_records.tab.gz
SAMPLE   CELL   ORIGINAL_PATH   ORIGINAL_SIZE   ORIGINAL_MD5   CLEANED_PATH   CLEANED_SIZE   CLEANED_MD5   DATE   STATUS   INFO
unique_sample_name   m84046_230412_163336_s3.hifi_reads   fastq/m84046_230412_163336_s3.hifi_reads.fastq.gz   16255055450   ca00950d807af12a5ee0e7ca6229cc44   results/reads_filtered/unique_sample_name/fastq/m84046_230412_163336_s3.hifi_reads.fastq_target-reads_with_reference_help-subset.fastq.gz   16539432159   932d1bf2f2fc65f01fc76640a2321308   2024-01-01T12:34:56   Overwritten   
```
Due to the current limitations of the algorithm, if the back-reference-qc pipeline is run again with the replaced cleaned fastq, be aware that there will be a further reduction of reads based on the zscore -2 criterion from reads that have already been filtered.

**Before performing overwriting, please check the write permissions of the files and folders.**

## To-do
- [ ] Put in CI tests
- [ ] Add conda enviroments
- [ ] Build bare-bone container to get minimal example running
- [ ] Expand on z_filter and new_fastq options
- [ ] Change to allow execution in units of fastq rather than by sample unit.

## FAQ
1. What is an example config and manifest?
   ```commandline
   sample  query_fofn      reference_fofn  comparison_type
   test_sample   fofn/test_sample.hifi.fastq.fofn   fofn/test_sample.illumina.fofn   self
   ```
2. Where can I find the Kraken2 database?
   1. `wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20230605.tar.gz`
   2. It contains Refeq archaea, bacteria, viral, plasmid, human1, UniVec_Core.
   3. You can find other types of databases here: https://benlangmead.github.io/aws-indexes/k2

## Overview
![pipeline vector](https://github.com/youngjun0827/back-reference-qc/blob/main/dag.svg)
