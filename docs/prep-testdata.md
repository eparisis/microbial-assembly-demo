# Benchmark Data Instructions

## Download from NCBI SRA

For SRA downloads you need to download the SRA Toolkit from [here](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit).

Configure your setting and run the following command to donwload a dataset:

```bash
prefetch SRR8073716
```

`--max-size` may be needed to download larger files than 20GB.

The data will be downloaded to the SRA configured directory from your setup.

Then run the following command to extraxt the fastq data:

```bash
fastq-dump SRR8073716 -O <outdir> --split-files --gzip
```

`fasterq-dump` will also split the fastq files into R1 and R2 files.
`--split-files` is required to split the fastq files into R1 and R2 files otherwise the data will be merged into a single file.

You can also use `parallel-fastq-dump` a wrapper to run `fastq-dump` in parallel. You can install it using conda:

`mamba install parallel-fastq-dump`

```bash
parallel-fastq-dump --sra-id SRR8073716 --threads 4 --outdir out/ --split-files --gzip --tmpdir tmp/
```

You will need to set a temporary directory for the download with `--tmpdir` to avoid errors with large files.

## Downsampling

For large datasets, you can downsample the data to a smaller size using [`rasusa`](https://github.com/mbhall88/rasusa) without introducing bias.

Downsampling was performed using the following command:

For the `mini datasets` a `0.2Mb` of data was chosen to be downsampled to.
For the `mid datasets` a `2Gb` of data was chosen to be downsampled to.

```bash
rasusa reads --bases 0.2Mb -o out1.fastq.gz -o out2.fastq.gz input_R1.fastq.gz input_R2.fastq.gz
```

To get `fastq` statistics you can use the following command:

```bash
seqkit stats -Ta -j 14 **/*.fastq.gz | tee Full_Fastq_Stats.tsv
```

This will fetch statistics for all `.fastq.gz` files in sub-directories and save the results to `Full_Fastq_Stats.tsv`.

# Benchmark Datasets

## Metagenomic Datasets

### DNASEQ Simulated Dataset

This Dataset is a CAMISIM simulated ONT dataset consisting of 73 strains. It is mainly used for profiling and classification benchmarks. Reference genomes can be found in the [GenomeTable](./DNASEQ_Simulated_reads/GenomeTable.tsv) while the profile composition can be found under [taxonomic_profile_0.txt](./DNASEQ_Simulated_reads/taxonomic_profile_0.txt).

### ZymoBIOMICS Microbial Community Standards

#### Zymo Even

Productname: [D6322](https://zymoresearch.eu/collections/zymobiomics-microbial-community-standards/products/zymobiomics-hmw-dna-standard).

*__Zymo_HMW_R104-GridION-EVEN__* is the Zymo even mock community deep sequenced from [this publication](https://www.nature.com/articles/s41592-022-01539-7?fromPaywallRec=false#Sec2) and can be downloaded from [this ENA link](https://www.ebi.ac.uk/ena/browser/view/PRJEB48692) (About ~49GB). This dataset can be subsetted to different sequencing depths.

The *__Zymo_HMW_EVEN_DNASEQ__* is a shallow sequencing run on the MinION from the DNASEQ AITOL_A project (Barcode 22).

#### Zymo Log

Productname: [D6310](https://zymoresearch.eu/collections/zymobiomics-microbial-community-standards/products/zymobiomics-microbial-community-standard-ii-log-distribution).

*__Zymo-GridION-LOG-BB-SN__* is a mid depth sequenced run from the Zymobiomics log distributed mock community HMW DNA standard. Dataset can be downloaded [here](https://github.com/LomanLab/mockcommunity?tab=readme-ov-file#data-availability) and is about ~16GB. An very deep sequncing run can also be downloaded from the same repo (~146GB). (Possible reference [here](https://academic.oup.com/gigascience/article/8/5/giz043/5486468?login=false#136473957))

#### Zymo Reference Genomes

Reference genome sequences for the Zymo set can be found in the [Reference_Genomes repo](./Zymo_Reference_Genomes/).

### BMock12

[Bmock12](https://github.com/Kirk3gaard/MicroBench) is a collection of mono cultures, and metagenomic samples. Data has been deposited to the ENA ([PRJEB85558](https://www.ebi.ac.uk/ena/browser/view/PRJEB85558)).

### MicroBench Metagenomic Datasets

[MicroBench](https://github.com/Kirk3gaard/MicroBench) is a dataset of 1000 metagenomic samples sequenced on the MinION platform. The dataset can be downloaded from [ENA](https://www.ebi.ac.uk/ena/browser/view/PRJEB85558).

## Microbial Isolate Datasets

### MicroBench Microbial Isolate Datasets

MicroBench newest datasets are R10 Promethion Datasets with the Rapid Barcoding kit (checkout the excel file for more details).

[Anabaena variabilis PCC 7120](https://github.com/Kirk3gaard/MicroBench#anabaena-variabilis-pcc-7120-dsm-107007) datasets can be found from R10 Promethion Datasets with Rapid Barcoding Kit and fast, hac and sup basecalling models aswell as raw pod5 files.
