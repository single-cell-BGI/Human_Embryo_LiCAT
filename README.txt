Source code of the manuscript “An integrated chromatin accessibility and transcriptome landscape of human pre-implantation embryos”Requirements
These scripts have been tested on Mac OSX, but should be supported by Linux or Windows systems as well. Users should have R version 3.4.0 or higher, and several packages set up from CRAN and Bioconductor, as indicated in the scripts. Installing R version 3.4.2 on Ubuntu 16.04the latest version of R can be installed by adding the latest repository to apt:sudo echo "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" | sudo tee -a /etc/apt/sources.list gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9 gpg -a --export E084DAB9 | sudo apt-key add -sudo apt-get updatesudo apt-get install r-base r-base-devwhich should install in about 20 seconds.Descriptions of the scripts:  
Number_of_peak.R - Count the number of peaks at each developmental stage, using MACS2 called narrow peaks as input. 
Peak_annotation.R - Annotate the peaks in each developmental stage based on their distribution along human genome. Including the proportion of peaks in promoters, exons, introns et al. 
CA_data_matrix.R - Generate chromatin accessibility count matrix of human embryos and perform data normalization using DESeq2.  
GE_data_matrix.R - Generate gene experession count matrix of human embryos and perform data normalization using DESeq2.  
PCA_CA.R - PCA analysis of human embryos using normalized chromatin accessibility datasets as input.
PCA_GE -   PCA analysis of human embryos using normalized gene expression datasets as input.
Numer_of_gene.R - Count the number of genes using the normalized data matrix as input.
CA_sample_correlation_human_embryos.R - Calculate correlations between chromatin accessibility profiles of each replicate of human embryos.
GE_sample_correlation_human_embryos.R - Calculate the correlations between gene expression profiles of each replicate of human embryos.
 Gene_density_genomic_windows.R - Calculate the correlations between gene density and chromatin accessibility in each developmental stage.
 CA_association_with_GC_ratio_human_embryos.R - Calculate the chromatin accessibility level at regions with different GC ratios in human embryos.
 CA_association_with_GC_ratio_hESC_differentiated_Cells.R - Calculate the chromatin accessibility at regions with different GC ratios in hESCs/differentiated cells.
 CA_sample_correlation_hESC_differentiated_Cells.R - Calculate the correlations between chromatin accessibility profiles of each replicate of hESC/differentiated cells.
GE_sample_correlation_hESC_differentiated_Cells.R - Calculate the correlations between gene expression profiles of each replicate of hESC/differentiated cells.
 CA_and_NANOG_OCT4_binding_site.R - Calculate the overlaps between chromatin accessibility peaks and NANOG/OCT4 binding sites in hESCs and differentiated cells.
 GE_NANOG_OCT4_expression.R - Gene expression level of NANOG/OCT4 in hESCs and differentiated cells.
 CA_sample_correlation_mouse_embryos.R - Calculate the correlations between chromatin accessibility profiles of each replicate of mouse embryos.
 GE_sample_correlation_mouse_embryos.R - Calculate the correlations between gene expression profiles of each replicate of mouse embryos.
 CA_GE_dynamics.R - Analyze the overall dynamics of EGA gene expression by mFuzzy, following by chromatin accessibility dynamics of EGA genes and ploting the expression and accessibility dynamics for example genes (including different group of genes, and DNA methyltransferase and demethylase.
  Motif_enrichment_each_stage.R - Plot the enrichment of motifs (p-values calculated by Homer) at each developmental stage.
 DUX4_target_genes_expression.R -  Plot The expression level of DUX4 targeted genes.
 CA_dynamics_ERVs.R - The overall chromatin accessibility dynamics of ERVs.
 GE_dynamics_ERVs.R - The overall gene expression dynamics of ERVs Data availability
All raw sequencing data will be made available upon request at the peer review stage, and accession codes will be available before publication. The input data for the scripts, along with the parsed data can be found through this link: https://pan.genomics.cn/ucdisk/s/VRnaym 