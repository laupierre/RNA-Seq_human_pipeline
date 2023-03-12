### This is the working place for developing a RNA-Seq pipeline based on singularity

### 2- Human RNA-Seq: PE Illumina sequencing, in reverse strand orientation

This is the development version deposited for the production team: the current version is v0.0.1_dev  
There are two methods available for quantification: STAR and salmon  

The singularity images are available for:  
BBMap version 39.01  
FastQC version 0.11.9  
featureCounts version 2.0.4  
multiqc version 1.14  
rseqc version 5.0.1  
salmon version 1.10.0  
sambamba version 0.8.2  
STAR version 2.7.10b  
R version 4.2.2 Patched  

The Singularity containers are in /projects/ncrrbt_share_la/dev_pipe and are the same used for the mouse RNA-Seq pipeline.  
The human indexes are located in: /projects/ncrrbt_share_la/dev_pipe2/  
See human_indexes for details on how the STAR and salmon indexes are made  

This code was tested to produce STAR or salmon counts. The differential gene expression analysis was not implemented.  

