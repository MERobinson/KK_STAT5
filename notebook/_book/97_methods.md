# Methods

## Clinical data
Gene expression RNA-seq and mutation WGS and WES data from The Cancer Genome Atlas LAML project were downloaded via the Broad Firehose repository (Run date: 2016-01-28) @cancer2013genomic. Data was processed in R to 173 patient samples with matched BCL6 expression data (gene-level normalised RSEM v2) & RAS mutational status. RAS mutational status was determined based on presence of absence of missense mutations in any RAS pathway gene, defined as: NRAS, KRAS, BRAF, FLT3, NF1 and PTPN11. Survival analysis was performed in R using the survival package (v3.2.7; @survival-book). P-values reported indicate Log-Rank test signficances, with BH multiple testr correction where appropriate.
