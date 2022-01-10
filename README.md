# UMB-study

This repository contains data and code for the UTI microbiome (UMB) study.

*umb_analysis.r* - Main analysis and figure generation.

## Data

*umb_participants.tsv* - table of all participants in the study

*umb_samples.tsv* - table of all samples collected during the study

*umb_stool.tsv* - table of all stool samples and corresponding questionnaire responses at time of collection

*umb_utis.tsv* - table of UTI events during the study

## Outputs

*ecoli.scc.midpointroot.nwk* - SCC phylogenetic tree of E. coli database

*ecoli.scc.umb_only.nwk* - SCC phylogenetic tree of E. coli references identified in this study

*ecoli_strain_clades.txt* - reference to phylogroup mapping

*strainge_table.tsv* - StrainGE E. coli results for all samples

*umb_metaphlan2.tsv* - MetaPhlan2 results for all stool samples

## Code

*umb_prep.r* - Data preparation

*umb_functions.r* - R functions
