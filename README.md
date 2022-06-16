# Activated monocytes recruit CD8 T cells to the cerebrospinal fluid during cognitive impairment
[Gate Lab at Northwestern University](https://sites.northwestern.edu/gatelab/)

## Abstract
Cerebrospinal fluid (CSF) contains a tightly regulated, specialized immune system. Yet, little is known about how CSF immunity is altered with aging or neurodegenerative disease. Here, we performed single cell RNA sequencing (scRNAseq) on CSF collected from 45 cognitively normal subjects ranging from 54-82 years old. We then assessed age-related transcriptomic changes using bioinformatic approaches, including linear and local polynomial regression. We reveal increased expression of lipid transport genes _Apolipoprotein E (APOE)_, _Apolipoprotein C1 (APOC1)_ and _Phospholipid transfer protein (PLTP)_ in activated monocytes with age. We then compared CSF immune systems from cognitively normal subjects to 14 subjects with mild cognitive impairment or Alzheimer’s disease. We detected upregulation of _C-X-C Motif Chemokine Receptor 6 (CXCR6)_ in clonally expanded T cells of cognitively impaired subjects. The CXCR6 ligand, _C-X-C Motif Chemokine Ligand 16 (CXCL16)_, was elevated in CSF and was associated with levels of neuroaxonal damage in cognitively impaired subjects. These results highlight dysregulation of lipid transport in CSF monocytes with age and identify the CXCR6-CXCL16 signaling axis as a potential route for T cell entry into brains with neurodegeneration. 

![merged_about_fig](https://user-images.githubusercontent.com/91904251/164067655-7c415284-46c4-42c9-8b06-af50763686fe.png)

## About
This repository contains code used to process and analyze scRNA-TCRseq data from the **Distinct cerebrospinal fluid immune perturbations in healthy brain aging and cognitive impairment** study. 

All ```R``` dependencies are listed in the ```renv.lock``` file. Upon opening the ```csf_aging.Rproj``` file in ```Rstudio```, the ```renv``` package should automatically download and activate. Afterwards, the user can run ```renv::restore()``` to download all necessary packages for the study. 

The dataset can be viewed and analyzed interactively in our [modified ShinyCell app](https://gatelabnu.shinyapps.io/csf_aging/).

Raw .fastq files and gene expression matrices can be downloaded from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200164), post-publicatoin

For any comments/questions about the study please refer to the [Gate Lab website](https://sites.northwestern.edu/gatelab/) for contact information.

#
![Feinberg-linear-RGB](https://user-images.githubusercontent.com/91904251/164067720-937687c0-874b-4aaa-afd4-76f887e07025.png)
