# non-canonical_activities_of_plant_miRNAs
This repository is for hosting the scripts and files to reproduce the analyses described in the listed manuscript "Expanding miRNA Regulatory Activities in Arabidopsis thaliana: What to Believe, and What to Leave?"\
Pranav Dawar<sup>1,2</sup>, Bhumika Jayee<sup>3</sup>, Isaac R. Eason<sup>4</sup>, Md. Fakhrul Azad<sup>1</sup>, Christopher D. Rock<sup>2</sup>*\
<sup>1</sup> Department of Biological Sciences, Texas Tech University, Lubbock, TX 79409, USA\
<sup>2</sup> Current address: Environmental Molecular Sciences Laboratory, Pacific Northwest National Laboratory, Richland, WA 99354, USA\
<sup>3</sup> Physical Sciences and Computational Division, Pacific Northwest National Laboratory, Richland, WA 99354, USA\
<sup>4</sup> Department of Chemistry and Biochemistry, Texas Tech University, Lubbock, TX 79409, USA\
*Correspondence. Email: chris.rock@ttu.edu\
Running head: Expanding miRNA Regulatory Activities in A. thaliana\
Key words: plant miRNA, canonical activity, non-canonical activity, Argonaute10, Post Transcriptional Gene Silencing, degradome analysis\

| File Name | Description |
| --- | --- |
| `PRJNA336058_abundance.txt` | kallisto output for fastq files used from project PRJNA336058 |
| `PRJNA388207_abundance.txt` | kallisto output for fastq files used from project PRJNA388207 |
| `PRJNA545832_abundance.txt` | kallisto output for fastq files used from project PRJNA545832 |
| `PRJNA788534_single_end_abundance.txt` | kallisto output for fastq files used from project PRJNA788534 |
| `RMSD_RNA.py` | python script to calculate RMSD for RNA duplexes within each conformation. The analysis focused exclusively on the RNA chains B (miRNA) and C (target mRNA) |
| `align_per_domain_plddt.py` | python script to calculate RMSD for protein (AGO10) domains using a custom Python workflow built on the Bio.PDB module from Biopython. The experimental cryo-EM structure (PDB ID: 7SWF) was used as the reference  |
| `com_distance_domain_rna.py` | python script to assess spatial relationship between individual protein domains and the RNA component |
| `run_align_per_domain_plddt.txt` | To run align_per_domain_plddt.py |
| `run_com_distance_domain_rna.txt` | To run com_distance_domain_rna.py |
| `structural_results_data_analysis.r` | R script used to generate Figures 4B to 4D in the main text |
| `cleaveland_output_data_analysis_Figures1Ato1F.r` | R script used to generate Figures 1A to 1F in the main text |
| `RNAseq_analysis_Figures3Ato3C.r` | R script used to generate Figures 3A to 3C in the main text |
| `degradome_results_all.zip` | CleaveLand output for degradome sequencing analysis reporting all slicing events with MFE ration >= 0.65 |
| `T_plots_MFE_0.6_0.65.txt` | CleaveLand output for degradome sequencing analysis reporting all slicing events with MFE ration >= 0.6 to < 0.65 |
