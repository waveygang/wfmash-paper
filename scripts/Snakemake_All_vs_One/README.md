# Snakemake for aligners

## Prepare
- Download fasta and HiFi reads
    - Rename all fasta into Chr1/Chr2|chr1|chr2, remove any unscaffolded sequences

## Folder structrure
```bash
.
├── Athaliana.metadata.tsv
├── Snakefile
├── aligner_parameters.tsv
├── fasta
│   ├── Athaliana
│   │   ├── Col-CC.fa.gz
│   │   ├── Cvi-0.fa.gz
│   │   └── Ler-0.fa.gz
│   └── Slycopersicum
│       ├── SL5.fa.gz
│       └── TS421.fa.gz
├── fastq
│   ├── Athaliana
│   │   ├── Cvi-0.hifi.fastq.gz
│   │   └── Ler-0.hifi.fastq.gz
│   └── Slycopersicum
│       ├── TS421.hifi.fastq.gz
│       └── TS60.hifi.fastq.gz
├── ref
│   ├── Athaliana
│   │   └── Col-CC.gene.gff3
│   └── Slycopersicum
│       └── SL5.gene.gff3
├── rules
│   └── all-to-one.smk
└── species.metadata.tsv
```
## Test
```
snakemake -npr -S Snakefile
```
You could see this output
```
Job stats:
job                   count    min threads    max threads
------------------  -------  -------------  -------------
AnchorWave_align          3              1              1
AnchorWave_maf2paf        3              1              1
AnchorWave_pre            2              1              1
all                       1              1              1
dv                        3              1              1
dv_prep                   2              1              1
faindex                   5              1              1
hifi_mapping              3              1              1
minimap2                  3              1              1
nucmer                    3              1              1
wfmash                    3              1              1
total                    31              1              1
```

## Run on Slurm

# use slurm profile
https://github.com/Snakemake-Profiles/slurm

