# Preparing Test Data

This guide covers how to download benchmark datasets and downsample them to create manageable test data for running assembly tools locally.

---

## Table of Contents

1. [Downloading from NCBI SRA](#downloading-from-ncbi-sra)
2. [Downsampling with Rasusa](#downsampling-with-rasusa)
3. [Fastq Statistics](#fastq-statistics)
4. [Benchmark Datasets](#benchmark-datasets)
   - [Metagenomic Datasets](#metagenomic-datasets)
   - [Microbial Isolate Datasets](#microbial-isolate-datasets)

---

## Downloading from NCBI SRA

For SRA downloads you need the SRA Toolkit installed from [here](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit).

### Prefetch

```bash
prefetch SRR8073716
```

Use `--max-size` for files larger than 20 GB. Data will be downloaded to the SRA configured directory from your setup.

### Extract FASTQ

```bash
fastq-dump SRR8073716 -O <outdir> --split-files --gzip
```

`--split-files` is required to split paired-end data into R1 and R2 files, otherwise everything merges into a single file. `fasterq-dump` also works but does not support `--gzip` directly.

### Parallel Download (Recommended)

For faster downloads, use `parallel-fastq-dump`:

```bash
mamba install -c bioconda parallel-fastq-dump
```

```bash
parallel-fastq-dump \
    --sra-id SRR8073716 \
    --threads 4 \
    --outdir out/ \
    --split-files \
    --gzip \
    --tmpdir tmp/
```

Set `--tmpdir` to avoid errors with large files filling up the default temp directory.

---

## Downsampling with Rasusa

For large datasets, downsample to a smaller size using [`rasusa`](https://github.com/mbhall88/rasusa) without introducing bias.

```bash
mamba install -c bioconda rasusa
```

### Suggested Sizes

| Dataset size | Use case |
|-------------|----------|
| `0.2Mb` | Mini — quick tool testing, CI/CD |
| `10Mb` | Small — local testing, parameter exploration |
| `2Gb` | Mid — realistic benchmarking |

### Paired-End Reads

```bash
rasusa reads --bases 0.2Mb -o out1.fastq.gz -o out2.fastq.gz input_R1.fastq.gz input_R2.fastq.gz
```

### Long Reads (Single File)

```bash
rasusa reads --bases 10Mb -o downsampled.fastq.gz input.fastq.gz
```

---

## Fastq Statistics

Check your data before and after downsampling:

```bash
seqkit stats -Ta -j 14 **/*.fastq.gz | tee Full_Fastq_Stats.tsv
```

This fetches statistics for all `.fastq.gz` files in sub-directories and saves results to `Full_Fastq_Stats.tsv`.

---

## Benchmark Datasets

### Metagenomic Datasets

#### ZymoBIOMICS Microbial Community Standards

The ZymoBIOMICS mock communities are widely used benchmarks with known composition, making them ideal for evaluating assembly and binning pipelines.

##### Zymo Even

Product: [D6322](https://zymoresearch.eu/collections/zymobiomics-microbial-community-standards/products/zymobiomics-hmw-dna-standard) — HMW DNA Standard (even distribution).

**Zymo_HMW_R104-GridION-EVEN** is the Zymo even mock community deep sequenced from [this publication](https://www.nature.com/articles/s41592-022-01539-7) and can be downloaded from [ENA PRJEB48692](https://www.ebi.ac.uk/ena/browser/view/PRJEB48692) (~49 GB). This dataset can be subsetted to different sequencing depths.

##### Zymo Log

Product: [D6310](https://zymoresearch.eu/collections/zymobiomics-microbial-community-standards/products/zymobiomics-microbial-community-standard-ii-log-distribution) — Log distribution.

**Zymo-GridION-LOG-BB-SN** is a mid-depth sequencing run of the log-distributed mock community HMW DNA standard. Download from the [LomanLab mock community repo](https://github.com/LomanLab/mockcommunity?tab=readme-ov-file#data-availability) (~16 GB, with a ~146 GB deep sequencing run also available). Reference: [Nicholls et al., GigaScience 2019](https://academic.oup.com/gigascience/article/8/5/giz043/5486468).

##### Zymo Reference Genomes

Reference genome sequences for the Zymo sets can be found in the [Reference_Genomes repo](./Zymo_Reference_Genomes/).

#### ZymoBIOMICS Fecal Reference (ONT PromethION)

Product: [D6323](https://www.zymoresearch.com/collections/zymobiomics-microbial-community-standards/products/zymobiomics-fecal-reference-with-trumatrix-technology) — Fecal Reference with TruMatrix Technology.

A high-depth metagenomic dataset from Oxford Nanopore, ideal for benchmarking metagenome assembly and binning (this is the dataset used in the MAG recovery figures in the main tutorial).

| Detail | Value |
|--------|-------|
| **Platform** | PromethION |
| **Flow cell** | FLO-PRO114M |
| **Chemistry** | R10.4.1 |
| **Basecall model** | v5.0.0 |
| **Modified bases** | 4mC, 5mC, 6mA |
| **Replicates** | 2 biological, 1 flow cell each |

| Data | Size |
|------|------|
| Raw (POD5) | 4.6 TB |
| Basecalls (BAM) | 850 GB |
| Analysis | 1.2 TB |

**Download:**

```bash
# Full dataset
aws s3 sync --no-sign-request s3://ont-open-data/zymo_fecal_2025.05 zymo_fecal_2025.05

# Or browse files first
# https://42basepairs.com/browse/s3/ont-open-data/zymo_fecal_2025.05
```

The analysis was performed with metaMDBG (v1.0), SemiBin2 (v2.1.0), MetaBAT2 (v2.17), DAS Tool (v1.1.7), and CheckM2 (v1.0.2). See the [ONT metagenomics application note](https://nanoporetech.com/resource-centre/oxford-nanopore-sequencing-provides-superior-metagenome-assembled-genome-recovery-and-strain-level-resolution-from-a-complex-microbiome) for details.

> **Tip**: This dataset is very large. For local testing, download only the basecalls and downsample with rasusa to 1-5 Gb.

#### Pathogen Surveillance (ONT GridION)

A metagenomic sequencing dataset for pathogen surveillance using Oxford Nanopore's rapid metagenomic surveillance protocol.

| Detail | Value |
|--------|-------|
| **Sample** | Human sputum spiked with Zeptometrix Respiratory Panel 2.1 |
| **Targets** | DNA viruses, RNA viruses, bacteria, atypical bacteria |
| **Platform** | GridION |
| **Flow cell** | FLO-MIN114 |
| **Chemistry** | R10.4.1 |
| **Basecall model** | HAC v4.3.0 (rebasecalled with Dorado v0.9.1) |
| **Replicates** | 3 technical per panel (2 panels) |

| Data | Size |
|------|------|
| Raw (POD5) | 255 GB |
| Basecalls (BAM) | 37 GB |
| Analysis | 13 MB |

**Download:**

```bash
# Full dataset
aws s3 sync --no-sign-request s3://ont-open-data/pathogen_surveillance_2025.09 pathogen_surveillance_2025.09

# Or browse files first
# https://42basepairs.com/browse/s3/ont-open-data/pathogen_surveillance_2025.09
```

Analysis was performed with [wf-metagenomics](https://github.com/epi2me-labs/wf-metagenomics) (v2.13.0). More details at the [epi2me dataset page](https://epi2me.nanoporetech.com/pathogen_surveillance_2025.09/).

> **Tip**: For local testing, download just the basecalls (37 GB) and downsample.

#### BMock12 / MicroBench Metagenomic Datasets

[BMock12](https://github.com/Kirk3gaard/MicroBench) is a collection of mono cultures and metagenomic samples. [MicroBench](https://github.com/Kirk3gaard/MicroBench) extends this with 1000 metagenomic samples sequenced on MinION. Data deposited at ENA: [PRJEB85558](https://www.ebi.ac.uk/ena/browser/view/PRJEB85558).

---

### Microbial Isolate Datasets

#### MicroBench Isolate Datasets

[MicroBench](https://github.com/Kirk3gaard/MicroBench) newest datasets are R10 PromethION datasets with the Rapid Barcoding Kit (check the repo for details).

[Anabaena variabilis PCC 7120](https://github.com/Kirk3gaard/MicroBench#anabaena-variabilis-pcc-7120-dsm-107007) datasets are available from R10 PromethION with Rapid Barcoding Kit, including fast, hac, and sup basecalling models as well as raw pod5 files.
