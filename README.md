# Single cell transcriptomics reveals cerebrospinal fluid immune dysregulation during healthy brain aging and cognitive impairment
[Gate Lab at Northwestern University](https://sites.northwestern.edu/gatelab/)

## Abstract
Cerebrospinal fluid (CSF) contains a tightly regulated immune system. Yet, little is known about how CSF immunity is altered with aging or neurodegenerative disease. Here, we performed single cell RNA sequencing on CSF from 45 cognitively normal subjects ranging from 54-82 years old. We reveal upregulation of lipid transport genes in monocytes with age. We then compared this cohort to 14 cognitively impaired subjects. In cognitively impaired subjects, downregulation of lipid transport genes in monocytes occurred concomitantly with altered cytokine signaling to CD8 T cells. Clonal CD8 T effector memory cells upregulated _C-X-C Motif Chemokine Receptor 6 (CXCR6)_ in cognitively impaired subjects. The CXCR6 ligand, _C-X-C Motif Chemokine Ligand 16 (CXCL16)_, was elevated in CSF of cognitively impaired subjects, suggesting CXCL16-CXCR6 signaling as a mechanism for antigen-specific T cell entry into the brain. Cumulatively, these results reveal cerebrospinal fluid immune dysregulation during healthy brain aging and cognitive impairment. 

![merged_about_fig](https://user-images.githubusercontent.com/91904251/175093470-eb5fec04-98d8-46d8-b05e-e89b477e4b4c.png)

## About
This repository contains code used to process and analyze scRNA-TCRseq data from the **Single cell transcriptomics reveals cerebrospinal fluid immune dysregulation during healthy brain aging and cognitive impairment** study. 

All ```R``` dependencies are listed in the ```renv.lock``` file. Upon opening the ```csf_aging.Rproj``` file in ```Rstudio```, the ```renv``` package should automatically download and activate. Afterwards, the user can run ```renv::restore()``` to download all necessary packages for the study. 

The dataset can be viewed and analyzed interactively in our [modified ShinyCell app](https://gatelabnu.shinyapps.io/csf_aging/).

Raw .fastq files and gene expression matrices can be downloaded from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200164), post-publicatoin

For any comments/questions about the study please refer to the [Gate Lab website](https://sites.northwestern.edu/gatelab/) for contact information.

#
![Feinberg-linear-RGB](https://user-images.githubusercontent.com/91904251/164067720-937687c0-874b-4aaa-afd4-76f887e07025.png)
