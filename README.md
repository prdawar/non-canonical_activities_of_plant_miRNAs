# non-canonical_activities_of_plant_miRNAs
This repository is for hosting the scripts and files to reproduce the analyses described in the listed manuscript "Expanding miRNA Regulatory Activities in Arabidopsis thaliana: What to Believe, and What to Leave?"\
Pranav Dawar<sup>1,2</sup>, Bhumika Jayee<sup>3</sup>, Isaac R. Eason<sup>4</sup>, Md. Fakhrul Azad<sup>1</sup>, Christopher D. Rock<sup>2</sup>*\
<sup>1</sup> Department of Biological Sciences, Texas Tech University, Lubbock, TX 79409, USA\
<sup>2</sup> Current address: Environmental Molecular Sciences Laboratory, Pacific Northwest National Laboratory, Richland, WA 99354, USA

<sup>3</sup> Physical Sciences and Computational Division, Pacific Northwest National Laboratory, Richland, WA 99354, USA

<sup>4</sup> Department of Chemistry and Biochemistry, Texas Tech University, Lubbock, TX 79409, USA

*Correspondence. Email: chris.rock@ttu.edu

Running head: Expanding miRNA Regulatory Activities in A. thaliana

Key words: plant miRNA, canonical activity, non-canonical activity, Argonaute10, Post Transcriptional Gene Silencing, degradome analysis

A description of each item in the repository:
limmaFit.R Wrapper function for LIMMA linear modeling. It preprocesses the data (e.g., drops rows/columns with all NAs, checks model specifications), fits linear models, and applies empirical Bayes moderation to assess differential expression across features. Implemented in "1_DEA.R".
limmaDEA.R Wrapper for performing differential analysis using the limma package, specifically through moderated t-tests or F-tests on coefficients derived from an MArrayLM object. It allows for flexible testing across multiple contrasts or coefficients, with options for global or individual adjustment of p-values to control the false discovery rate. Implemented in "1_DEA.R".
0a_loading_intensity_data.R Converts the output of TopPIC searches to an MSnSet object that stores label-free quantification data. Final data is on a relative log2-scale.
0b_loading_spectralcount_data.R Converts the output of TopPIC searches to an MSnSet object that stores spectral count data.
1_DEA.R Takes output of script "0a" and performs differential abundance analysis using the functions "limmaFit.R" and "limmaDEA.R".
Figure1_overview_plots.R Takes output of script "0a" and generates Figure 1, including quantifying the number of distinct and average proteoforms/genes and the completeness of those identifications. Also generates Figure S1. Histogram of observed proteoform monoisotopic masses.
Figure2_INS_GCG_rectangle_plots.R Takes output of script "0b" and generates Figure 2: Rectangle plots of most frequently observed proteoforms of insulin and glucagon.
Figure3_volcano_plots_rawp.R Takes output of script "1_DEA" and generates Figure 3: Volcano plots of all and curated proteoforms.
Figure4_CXCL.R Takes output of script "0a" and generates Figure 4: Heatmap and table summary of observed CXCL proteoforms. Performs hypergeometric probability distribution test to determine which unique observations are statistically significant.
FigureSI_HMGN.R Takes output of script "0a" and generates Figure S7: Heatmap and table summary of observed HMGN proteoforms.
FigureSI_RSD.R Takes output of script "0a" and generates Figure S2: Histogram of relative standard deviations based on label-free quantification data. Facets data into control and cytokine treated samples.
FigureSI_rectangle_plots.R Takes output of script "0b" and generates Figure S##-S##: Rectangle plots of most frequently observed proteoforms of chromogranin-A, secretogranin-1 a.k.a. chromogranin-B, secretogranin-2 a.k.a. chromogranin-C somatostatin, pancreatic polypeptide prohormone, islet amyloid polypeptide, and VGF.
Figure_1A_tdislets.png Graphical abstract for project. This file is not needed for analysis.
workflow_github_v2.png Workflow figure for running scripts in repository. This file is not needed for analysis.
