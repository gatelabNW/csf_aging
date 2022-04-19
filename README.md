# Distinct cerebrospinal fluid immune perturbations in healthy brain aging and cognitive impairment
[Gate Lab at Northwestern University](https://sites.northwestern.edu/gatelab/)

## Abstract
Cerebrospinal fluid (CSF) contains a tightly regulated, specialized immune system. Yet, little is known about how aging influences CSF immunity in cognitively typical versus cognitively impaired individuals. Here, we performed single cell RNA sequencing (scRNAseq) on CSF collected from 45 cognitively typical subjects ranging from 54-82 years old. We then assessed age-related transcriptomic changes using bioinformatic approaches, including linear and local polynomial regression. We reveal pronounced changes to several CSF immune cell types, underscored by increased expression of lipid processing genes in activated monocytes with age. We then compared CSF immune systems from cognitively typical subjects to 14 subjects with mild cognitive impairment or Alzheimer’s disease. Our results indicate disparate age-related CSF immune system perturbations, including the upregulation of C-X-C Motif Chemokine Receptor 6 in clonally expanded T cells of cognitively impaired subjects. These results highlight the utility of CSF immune system changes to identify therapeutic targets of neurodegenerative disease-associated neuroinflammation. 

![merged_about_fig](https://user-images.githubusercontent.com/91904251/164067655-7c415284-46c4-42c9-8b06-af50763686fe.png)

## About
This repository contains code used to process and analyze scRNA-TCRseq data from the **Distinct cerebrospinal fluid immune perturbations in healthy brain aging and cognitive impairment** study. 

All ```R``` dependencies are listed in the ```renv.lock``` file. Upon opening the ```csf_aging.Rproj``` file in ```Rstudio```, the ```renv``` package should automatically download and activate. Afterwards, the user can run ```renv::restore()``` from the project root to download all necessary packages for the study. 

The dataset can be viewed and analyzed interactively in our [modified ShinyCell app](https://gatelabnu.shinyapps.io/csf_aging/).

Raw .fastq files and gene expression matrices can be downloaded from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200164) post-publicatoin

For any comments/questions about the study please refer to the [Gate Lab website](https://sites.northwestern.edu/gatelab/) for contact information.

#
![Feinberg-linear-RGB](https://user-images.githubusercontent.com/91904251/164067720-937687c0-874b-4aaa-afd4-76f887e07025.png)
