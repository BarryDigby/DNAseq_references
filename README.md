# DNAseq References for GRCh37, 38
## Create reference data for [Somatic_n-of-1 NextFlow pipeline](https://github.com/brucemoran/somatic_n-of-1)
### How to Setup
#### Dependencies:
[NextFlow](https://www.nextflow.io/index.html#GetStarted), [Singularity](https://sylabs.io/guides/3.0/user-guide/installation.html#)
#### Reference Generation:
```
nextflow run brucemoran/DNAseq_references
```
```
Optional arguments:
  --version   STRING    GRCh37 or GRCh38 (default)
  --outDir    STRING    output directory path; NB ${params.version} dir is created therein
  --exometag    STRING    naming for exome outputs when supplied; tag is then used in somatic_n-of-1 and batch_somatic pipelines to select relevant exome data
  and either
  --exomebedurl     STRING      URL to exome bed file for intervals; NB assumes GRCh37
  or
  --exomebedfile     STRING      locally downloaded exome bed file for intervals; NB assumes GRCh37
```
