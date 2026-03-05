# Microbial Genome Assembly – MSc Tutorial

A hands-on walkthrough for assembling a microbial genome from Oxford Nanopore long reads using **Flye**, managed entirely with **conda / mamba**.

---

## Table of Contents

1. [Introduction](#1-introduction)
2. [Learning Objectives](#2-learning-objectives)
3. [Prerequisites](#3-prerequisites)
4. [Environment Setup](#4-environment-setup)
5. [Obtaining the Data](#5-obtaining-the-data)
6. [Read Quality Control](#6-read-quality-control)
7. [Genome Assembly with Flye](#7-genome-assembly-with-flye)
8. [Assembly Quality Assessment](#8-assembly-quality-assessment)
9. [Visualising the Assembly Graph](#9-visualising-the-assembly-graph)
10. [Further Reading & Next Steps](#10-further-reading--next-steps)

---

## 1. Introduction

Long-read sequencing technologies (Oxford Nanopore Technologies, Pacific Biosciences) have transformed microbial genomics by producing reads that span repetitive regions and produce highly contiguous assemblies, often finishing a bacterial genome in a single contig.

In this tutorial you will assemble a small bacterial genome end-to-end:

```
Raw reads → QC → Assembly (Flye) → Assessment (QUAST / BUSCO) → Graph (Bandage)
```

All software is installed through **conda / mamba** so the environment is fully reproducible.

---

## 2. Learning Objectives

By the end of this tutorial you will be able to:

- Create and manage reproducible bioinformatics environments with conda / mamba.
- Assess the quality of Oxford Nanopore reads with NanoPlot.
- Assemble a microbial genome with Flye and understand its key parameters.
- Evaluate assembly contiguity and completeness with QUAST and BUSCO.
- Inspect the assembly graph with Bandage.

---

## 3. Prerequisites

| Requirement | Notes |
|---|---|
| Linux or macOS terminal | Windows users: use WSL 2 |
| conda **or** mamba installed | See [Miniforge](https://github.com/conda-forge/miniforge) for a fast installation |
| ~5 GB free disk space | For raw data, tools, and output |
| Basic command-line familiarity | `cd`, `ls`, `mkdir`, pipe `|` |

> **Tip – mamba vs conda:** `mamba` is a drop-in replacement for `conda` that resolves environments much faster. Simply replace `conda` with `mamba` in any command below if you have mamba available.

---

## 4. Environment Setup

### 4.1 Create the conda environment

We install all required tools in a single, named environment.

```bash
conda create -n assembly -c conda-forge -c bioconda \
    flye \
    nanoplot \
    nanostat \
    quast \
    busco \
    bandage \
    sra-tools \
    -y
```

> Using mamba? Replace `conda create` with `mamba create` for faster dependency resolution.

### 4.2 Activate the environment

```bash
conda activate assembly
```

You should see `(assembly)` at the start of your prompt. All subsequent commands assume this environment is active.

### 4.3 Verify tool installation

```bash
flye --version
NanoPlot --version
quast --version
busco --version
```

---

## 5. Obtaining the Data

We will use a publicly available Nanopore dataset of *Escherichia coli* (a small, well-characterised genome of ~5 Mb) from NCBI SRA.

### 5.1 Download the reads

```bash
mkdir -p ~/assembly_tutorial/raw_data
cd ~/assembly_tutorial/raw_data

# Download SRA accession SRR23415347 (~300 MB)
prefetch SRR23415347
fasterq-dump SRR23415347 --outdir .
```

> `prefetch` + `fasterq-dump` are both part of **sra-tools** (already installed).
> Depending on your internet connection this may take several minutes.

### 5.2 Inspect the download

```bash
ls -lh
head -8 SRR23415347.fastq
```

Expected output: FASTQ records with four lines each (header, sequence, `+`, quality).

---

## 6. Read Quality Control

Before assembly it is important to understand the characteristics of your reads: length distribution, quality scores, and any potential issues.

### 6.1 Run NanoPlot

```bash
mkdir -p ~/assembly_tutorial/qc
cd ~/assembly_tutorial

NanoPlot \
    --fastq raw_data/SRR23415347.fastq \
    --outdir qc/ \
    --threads 4 \
    --plots dot
```

### 6.2 Interpret the results

Open `qc/NanoPlot-report.html` in your browser. Pay attention to:

| Metric | What to look for |
|---|---|
| **N50 read length** | Higher is better; >10 kb is typical for a good Nanopore run |
| **Mean read quality** | Q10+ is acceptable; Q15+ is excellent |
| **Total bases** | Should be ≥50× coverage of the genome (50 × 5 Mb = 250 Mb) |

> **Quick summary without plots:**
> ```bash
> NanoStat --fastq raw_data/SRR23415347.fastq
> ```

### 6.3 Questions to consider

- What is the estimated sequencing depth (coverage)?
- Are there very short reads that might reduce assembly quality?
- What is the modal read length?

---

## 7. Genome Assembly with Flye

[Flye](https://github.com/mikolmogorov/Flye) is a de-novo assembler for single-molecule sequencing reads that builds an assembly graph and then resolves it into contigs.

### 7.1 Run Flye

```bash
mkdir -p ~/assembly_tutorial/assembly
cd ~/assembly_tutorial

flye \
    --nano-raw raw_data/SRR23415347.fastq \
    --out-dir assembly/ \
    --threads 4 \
    --genome-size 5m
```

**Key parameters explained:**

| Parameter | Meaning |
|---|---|
| `--nano-raw` | Input reads are Oxford Nanopore with standard accuracy. For reads basecalled to ≥Q20 accuracy (Guppy high-accuracy model / Dorado), use `--nano-hq` instead |
| `--out-dir` | Output directory |
| `--threads` | Number of CPU threads; adjust to your machine |
| `--genome-size` | Approximate expected genome size (helps Flye estimate coverage) |

> **High-quality basecalling?** If your reads were basecalled with a recent Dorado model (accuracy ≥ Q20), use `--nano-hq` instead of `--nano-raw` for a better assembly.

### 7.2 Monitor progress

Flye prints progress to the terminal. A typical *E. coli* run on 4 threads takes **10–30 minutes**.

Flye will iterate through several stages:
1. Overlap detection
2. Assembly graph construction
3. Graph polishing
4. Consensus generation

### 7.3 Examine the output

```bash
ls -lh assembly/
```

| File | Description |
|---|---|
| `assembly.fasta` | Final assembly sequences (contigs / scaffolds) |
| `assembly_graph.gfa` | Assembly graph in GFA format |
| `assembly_info.txt` | Per-contig statistics (length, coverage, circularity) |
| `flye.log` | Full log of the run |

```bash
# Verify assembly.fasta was created successfully
ls -lh assembly/assembly.fasta

# Quick summary of contigs
grep "^>" assembly/assembly.fasta | wc -l        # number of contigs
cat assembly/assembly_info.txt                    # detailed per-contig stats
```

**What does a good assembly look like?**
For a typical bacterial genome you should see 1–5 contigs, with the largest contig representing the chromosome. Plasmids appear as smaller, often circular, contigs.

---

## 8. Assembly Quality Assessment

### 8.1 QUAST – contiguity metrics

[QUAST](https://github.com/ablab/quast) reports standard assembly statistics.

```bash
mkdir -p ~/assembly_tutorial/quast_results
cd ~/assembly_tutorial

quast assembly/assembly.fasta \
    --output-dir quast_results/ \
    --threads 4
```

Open `quast_results/report.html` and review:

| Metric | Ideal value for a complete bacterial genome |
|---|---|
| **# contigs** | 1–5 |
| **Largest contig** | ~5 Mb (full chromosome) |
| **N50** | Close to the chromosome length |
| **L50** | 1 |
| **Total length** | ~5 Mb |

> If you have a reference genome available (e.g., downloaded from NCBI), provide it with `-r reference.fasta` for a more detailed comparison.

### 8.2 BUSCO – completeness assessment

[BUSCO](https://busco.ezlab.org/) checks how many conserved single-copy orthologous genes are present in your assembly, giving a proxy for completeness.

```bash
mkdir -p ~/assembly_tutorial/busco_results
cd ~/assembly_tutorial

busco \
    -i assembly/assembly.fasta \
    -o busco_results \
    -m genome \
    -l enterobacterales_odb10 \
    --cpu 4
```

> The lineage dataset (`enterobacterales_odb10`) will be downloaded automatically on first run.
> Choose the most specific lineage available for your organism; run `busco --list-datasets` to see all options.

Examine the short summary:

```bash
cat busco_results/short_summary.specific.enterobacterales_odb10.busco_results.txt
```

A high-quality assembly should show **≥95% Complete BUSCOs**.

---

## 9. Visualising the Assembly Graph

The assembly graph reveals how contigs are connected and whether the chromosome is fully resolved.

### 9.1 Open Bandage

[Bandage](https://rrwick.github.io/Bandage/) is a GUI tool for viewing assembly graphs.

```bash
Bandage &
```

### 9.2 Load the graph

1. Click **File → Load graph**.
2. Navigate to `~/assembly_tutorial/assembly/assembly_graph.gfa` and open it.
3. Click **Draw graph**.

### 9.3 Interpret the graph

| Pattern | Interpretation |
|---|---|
| Single circular node | Fully assembled circular chromosome – excellent! |
| Linear node with no dead ends | Complete chromosome in linear representation |
| Multiple connected nodes | Unresolved repeat regions; normal for complex genomes |
| Many disconnected fragments | Fragmented assembly; may need more coverage or polishing |

You can colour nodes by depth (coverage) to spot plasmids (often higher coverage than the chromosome).

---

## 10. Further Reading & Next Steps

### Polishing (optional)

Nanopore assemblies can be further improved by polishing with short Illumina reads if available:
- [Medaka](https://github.com/nanoporetech/medaka) – Nanopore-only polishing
- [Pilon](https://github.com/broadinstitute/pilon) – Illumina-based polishing

### Annotation

Once you are satisfied with your assembly, you can annotate it:
- [Prokka](https://github.com/tseemann/prokka) – fast prokaryotic annotation
- [Bakta](https://github.com/oschwengers/bakta) – more up-to-date alternative to Prokka

### Key references

- **Flye:** Kolmogorov M. et al. (2019). Assembly of long, error-prone reads using repeat graphs. *Nature Methods*, 16, 792–794. <https://doi.org/10.1038/s41592-019-0072-4>
- **QUAST:** Gurevich A. et al. (2013). QUAST: quality assessment tool for genome assemblies. *Bioinformatics*, 29(8), 1072–1075. <https://doi.org/10.1093/bioinformatics/btt086>
- **BUSCO:** Manni M. et al. (2021). BUSCO Update: Novel and Streamlined Workflows along with Broader and Deeper Phylogenetic Coverage. *Molecular Biology and Evolution*, 38(10), 4647–4654. <https://doi.org/10.1093/molbev/msab199>
- **Bandage:** Wick R.R. et al. (2015). Bandage: interactive visualization of de novo genome assemblies. *Bioinformatics*, 31(20), 3350–3352. <https://doi.org/10.1093/bioinformatics/btv383>

### Cheat sheet

```bash
# One-liner to check coverage before assembly
NanoStat --fastq raw_data/SRR23415347.fastq | grep -E "Total bases|Number of reads|Mean read length"

# Rerun Flye with high-quality flag
flye --nano-hq raw_data/SRR23415347.fastq --out-dir assembly_hq/ --threads 4 --genome-size 5m

# Activate / deactivate environment
conda activate assembly
conda deactivate
```

---

*Tutorial prepared for MSc Bioinformatics students. Please raise an issue in this repository if you encounter problems.*
