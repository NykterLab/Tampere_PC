# bpnet

Reproduce the analysis of chromatin accessibility using the BPNET method to discover TF-binding motifs.
The workflow is managed with Snakemake. For convenience, each part of the analysis can also be performed running the following scripts:

### Setup

Prepare all the input data:

- Run `generate_config_files.sh` to generate all the BPNET configuration files and BPNET data specification files used to train the models.
- Run `download_datasets.sh` to download all the required sequencing data, reference genome, and the UniBindPWM database.
- Run `process_input_data.sh` to perform all the required reads alignments, peak calling, data indexing, and selection of the genomic regions to be used for BPNET models training.

### BPNET models training

Build BPNET models:

- Run `train_chip_seq_BPNET.sh` to train a BPNET model using VCaP AR Chip-seq data.
- Run `train_VCaP_atac_seq_BPNET.sh` to train BPNET models using VCaP ATAC-seq data. The training regions for each model are selected from the binding sites of a different TF among AR, FOXA, HOXB13, and ERG.
- Run `train_PC_atac_seq_BPNET.sh` to train BPNET models using prostate cancer clinical tissue samples ATAC-seq data. The training regions for each model are selected from the binding sites of a different TF among AR, FOXA, HOXB13, and ERG.

### BPNET models interpretation

Interpret trained BPNET models:

- Run `motif_discovery_chip_seq_BPNET.sh` to compute contribution scores and run TF-MoDISCo for the VCaP AR Chip-seq BPNET model.
- Run `motif_discovery_VCaP_atac_seq_BPNET.sh` to compute contribution scores and run TF-MoDISCo for the VCaP ATAC-seq BPNET models.
- Run `motif_discovery_PC_atac_seq_BPNET.sh` to compute contribution scores and run TF-MoDISCo for the clinical samples ATAC-seq BPNET models.

### Motifs discovery

Analyse the sequence motifs resulting from the interpretation of BPNET models:

- Run `analyse_chip_seq.sh` to plot the obtained PFMs, compare them to known PWMs collected in the HOCOMOCO database, and generate a report for the VCaP AR Chip-seq BPNET model.
- Run `analyse_VCaP_atac_seq.sh` to plot the obtained PFMs, compare them to known PWMs collected in the HOCOMOCO database, and generate a report for the VCaP ATAC-seq BPNET models.
- Run `analyse_PC_atac_seq.sh` to plot the obtained PFMs, compare them to known PWMs collected in the HOCOMOCO database, and generate a report for the clinical samples ATAC-seq BPNET models.

## Directories

Executing the code will create new directories. The complete directory tree up to the third level is shown below:

```
.
├── genomes ......................... reference genome files
├── meme_db ......................... HOCOMOCO human full database in MEME format
├── model_config..................... BPNET models configuration gin files 
│   ├── atac_seq
│   │   ├── tampere_pc
│   │   └── VCaP
│   └── chip_seq
│       └── VCaP_SRA012454
├── model_dataspec .................. BPNET models data specification yml files
│   ├── atac_seq
│   │   ├── tampere_pc
│   │   └── VCaP
│   └── chip_seq
│       └── VCaP_SRA012454
├── models .......................... BPNET models - results are stored under each model directory 
│   ├── atac_seq
│   │   ├── tampere_pc
│   │   └── VCaP
│   └── chip_seq
│       └── VCaP_SRA012454
├── peaks ........................... chromatin accessibility peaks files in BED format
│   ├── atac_seq
│   │   ├── tampere_pc
│   │   └── VCaP
│   └── chip_seq
│       └── VCaP_SRA012454
├── reads ........................... sequencing data files in FASTQ and BAM formats
│   ├── atac_seq
│   │   ├── tampere_pc
│   │   └── VCaP
│   └── chip_seq
│       └── VCaP_SRA012454
├── src ............................. Python scripts
├── tracks .......................... sequencing data in BigWig format used to train the models
│   ├── atac_seq
│   │   ├── tampere_pc
│   │   └── VCaP
│   └── chip_seq
│       └── VCaP_SRA012454
├── training_regions ................ peaks files in BED format used to train the models
│   ├── atac_seq
│   │   ├── tampere_pc
│   │   └── VCaP
│   └── chip_seq
│       └── VCaP_SRA012454
└── UniBind_tfbs .................... UniBindPWM TF binding sites data files in BED format
    ├── UniBindPWM
    │   ├── AHR
    │   ├── AR
    │   ├── ARID3A
    │   ├── ARNT
    │   ├── ARNTL
    │   ├── ASCL1
    │   ├── ATF1
    │   ├── [...]
    └── UniBindPWM_merged_by_TF
```

## Software requirements

- [Snakemake](https://snakemake.readthedocs.io/en/stable/)
- [SRA Toolkit](https://ncbi.github.io/sra-tools/)
- [BOWTIE2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [SAMTools](http://www.htslib.org/)
- [BEDTools](https://bedtools.readthedocs.io/en/latest/)
- [BPNET](https://github.com/kundajelab/bpnet)
- [TOMTOM](http://meme-suite.org/tools/tomtom)
