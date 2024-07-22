# Bioinformatics: <em> Arabidopsis thaliana </em> Project

## Abstract
We report upregulation in genes related to both strengthening pathogen defenses and innate
immune responses, as well as downregulation in genes related to growth development in
Arabidopsis thaliana infected with Pseudomonas syringae pv tomato DC3000. Analyzing
patterns in gene expression following infection of the model plant A. thaliana by P. syringae
pv tomato DC3000 can provide a basis for similar plant-pathogen interactions, as well as
provide insight on the potential effects that a P. syringae outbreak would have on local
agriculture. RNA-Seq data from Howard et al. (2013) was utilized for this study, alongside
the A. thaliana TAIR10.1 reference genome and associated annotations from the NCBI
database. MOCK and VIR samples at 1 hr, 6 hrs, and 12 hrs post-infection were analyzed.
STAR was used to index the reference genome and align paired-reads to the reference
genome in Bash. The resultant gene count data was uploaded into R Studio for count
normalization using DeSeq2 and GO term enrichment using topGO. Analysis of the
over-represented GO terms in A. thaliana following infection by P. syringae allowed for the
elucidation of affected cellular mechanisms.