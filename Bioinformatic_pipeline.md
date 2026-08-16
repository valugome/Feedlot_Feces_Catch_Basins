# BIOINFORMATIC PIPELINE 
## Trimmomatic 
Trimmomatic was run within the AMR++ pipeline. Installation followed the tutorial: https://github.com/Microbial-Ecology-Group/AMRplusplus/blob/master/docs/installation.md

The ‘trim_qc’ subworkflow from the ‘main_AMR++.nf’ pipeline was run:
```
nextflow run main_AMR++.nf --pipeline trim_qc --output AMR++_results
```
