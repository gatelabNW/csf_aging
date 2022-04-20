# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----        Distinct cerebrospinal fluid immune perturbations           -----
# -----          in healthy brain aging and cognitive impairment           -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 10-07-2021
# Written by: Natalie Piehl
# Summary: Run cellranger on scRNA-TCRseq data
#
#-------------------------------------------------------------------------------
#
# - cellranger v6.0.0 was used for this study
# - Each RNA sample was processed with "cellranger count" and
#   each TCR sample was processed with "cellranger vdj"
# - GRCh38 reference can be downloaded from the 10x Genomics website:
#   https://support.10xgenomics.com/single-cell-vdj/software/downloads/latest
# - Raw .fastqs can be downloaded from GEO with following the link:
#   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200164
#
# The "cellranger count" command was executed as follows:

cellranger count \
--id <sample_id> \
--fastqs <path_to_sample_fastq> \
--transcriptome <path_to_GRCh38_ref> \
--sample <sample_fastq_filename_prefix> \
--expect-cells <expected_cellnumber> \
--localcores 8

# The "cellranger vdj"" command was executed as follows:

cellranger vdj \
--id <sample_id> \
--fastqs <path_to_sample_fastq> \
--reference <path_to_GRCh38_ref> \
--sample <sample_fastq_filename_prefix> \
--localcores 8