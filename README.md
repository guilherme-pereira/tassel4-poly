# Tassel4-Poly
Modified version of Tassel4 (v.4.3.7) for running the Tassel-GBS pipeline modified for polyploid species with high read depths

## Overview

[Tassel](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btm308) is a well-known Java program that has also implemented a genotyping-by-sequencing (GBS) analysis pipeline called [Tassel-GBS](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0090346). This pipeline has been broadly used by plant geneticists who want to use GBS-based single nucleotide polymorphism (SNP) in their genetic studies. However, one limitation for its use in polyploid species is that the pipeline provides an approximate or truncated number of read counts (up to 127), which are ultimately used to call diploidized genotypes. We modified Tassel-GBS pipeline so that it could store and return the actual read counts per individual per locus. 

## Running the pipeline

In the same way as Tassel3, Tassel4 is run from the terminal with command lines such as those found in this [tutorial](https://bytebucket.org/tasseladmin/tassel-5-source/wiki/docs/TasselPipelineGBS.pdf). Only two plugins have been directly affected by the modification: `FastqToTBT` and `DiscoverySNPCaller`. In order to get the exact read counts, one should change the flag `-y` to `-sh` so that a maximum number of 32,767 reads will be recorded (instead of up to 127 with the `-y`). This number should suffice the GBS experiments with higher depths demanded by polyploid species studies. 

**_Example commands_** (adapted from the [tutorial](https://bytebucket.org/tasseladmin/tassel-5-source/wiki/docs/TasselPipelineGBS.pdf)):

- `FastqToTBTPlugin`
```
./tassel4-poly/run_pipeline.pl -fork1 -FastqToTBTPlugin -i fastq -k myGBSProject_key.txt -e ApeKI -o tbt -sh –t mergedTagCounts/myMasterTags.cnt -endPlugin -runfork1
```
- `DiscoverySNPCallerPlugin`
```
./tassel4-poly/run_pipeline.pl -fork1 -DiscoverySNPCallerPlugin -i mergedTBT/myStudy.tbt.shrt -sh -m topm/myMasterTags.topm -mUpd topm/myMasterTagsWithVariants.topm -vcf -o vcf/raw/myGBSGenos_chr+ -mnF 0.8 -p myPedigreeFile.ped -mnMAF 0.02 -mnMAC 100000 -ref MyReferenceGenome.fa -sC 1 -eC 10 -endPlugin -runfork1
```

Remember to change the file extension `*.byte` to `*.shrt` where needed in subsequent plugins. Also, notice that VCF files returned by Tassel4-Poly still contains genotype calls as in a diploid species. In order to get quantitative (dosage) genotype calls for polyploid species as provided by [SuperMASSA](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0030906) software, you may want to use the VCF files as input data in the [VCF2SM](https://github.com/guilherme-pereira/vcf2sm) pipeline.

## Cite

If you end up using Tassel4-Poly in your study, please make sure to cite the publications below in your paper: 

- Glaubitz JC, Casstevens TM, Lu F, Harriman J, Elshire RJ, Sun Q, Buckler ES. (2014) TASSEL-GBS: A High Capacity Genotyping by Sequencing Analysis Pipeline. *PLoS ONE* 9(2): e90346. https://doi.org/10.1371/journal.pone.0090346

- Pereira GS, Garcia AAF, Margarido GRA. (2018) A fully automated pipeline for quantitative genotype calling from next generation sequencing data in autopolyploids. *BMC Bioinformatics* 19:398. https://doi.org/10.1186/s12859-018-2433-6.
