
[![DOI](https://zenodo.org/badge/99698135.svg)](https://zenodo.org/badge/latestdoi/99698135)  

## A list of R environment based tools for 16S rRNA gene data exploration, statistical analysis and visualization  
  
As a beginner, the entire process from sample collection to analysis for sequencing data is a daunting task. More specifically, the downstream processing of raw reads is the most time consuming and mentally draining stage. It is vital to understand the basic concepts in microbial ecology and then to use various tools at disposal to address specific research questions. Thankfully, several young researchers supported by their experienced principal investigators/supervisors are working on creating various tools for analysis and interpretation of microbial community data. A major achievement of the scientific community is the open science initiative which has led to sharing of knowledge worldwide. For microbial community analysis, several tools have been created in R, a free to use (GNU General Public License) programming language(Team, 2000). The power of R lies in its ease of working with individuals lacking programming skills and easy sharing of analysis scripts codes and packages aiding reproducibility. Using tools such as QIIME (the newer QIIME2) (Caporaso, Kuczynski, Stombaugh et al., 2010), Mothur (Schloss, Westcott, Ryabin et al., 2009), DADA2 (Callahan, McMurdie, Rosen et al., 2016) one can get from raw reads to species × samples table (OTU or ASVs amplicon sequence variants as suggested recently (Callahan, McMurdie & Holmes, 2017)). In this post, numerous resources that can be helpful for analysis of microbiome data are listed. This list may not have all the packages as this tool development space is ever growing. Feel free to add those packages or links to web tutorials related to microbiome data, there is a [google docs excel sheet at this link for a list of tools](https://docs.google.com/spreadsheets/d/1am-UyDVBGDOgm6jVQ5FDXxmg24iriHqeBeul14HRb1g/edit?usp=sharing) which can be edited to include more tools. These are mostly for improving statistical analysis and visualisation. These tools provide convenient options for data analysis and include several steps where the user has to make decisions. The work by [McMurdie PJ, Holmes S](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003531), [Weiss S](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0237-y) and [Tsilimigras M.C. and Fodor A.A](http://www.sciencedirect.com/science/article/pii/S1047279716300722) are useful resources to understand the data common to microbiome census. It can be tricky and frustrating in the beginning but patience and perseverance will be fruitful at the end (personal experience).

-------------------------------------------------------   

### Tools:
1.	Ampvis2	[Tools for visualising amplicon sequencing data](https://madsalbertsen.github.io/ampvis2/)  
2.	CCREPE	[Compositionality Corrected by PErmutation and REnormalization](http://bioconductor.org/packages/release/bioc/html/ccrepe.html)  
3.	DADA2	[Divisive Amplicon Denoising Algorithm](https://www.nature.com/nmeth/journal/v13/n7/full/nmeth.3869.html)  
4.	DESeq2	[Differential expression analysis for sequence count data](https://www.bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-mixture-models.html)  
5.	edgeR	[empirical analysis of DGE in R](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2796818/)  
6.	mare	[Microbiota Analysis in R Easily](https://github.com/katrikorpela/mare)  
7.	Metacoder	[An R package for visualization and manipulation of community taxonomic diversity data](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005404)  
8.	metagenomeSeq	[Differential abundance analysis for microbial marker-gene surveys](http://www.nature.com/nmeth/journal/v10/n12/full/nmeth.2658.html)  
9.	microbiome R package	[Tools for microbiome analysis in R](https://github.com/microbiome/microbiome)  
10.	MINT	[Multivariate INTegrative method](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1553-8)  
11.	mixDIABLO	[Data Integration Analysis for Biomarker discovery using Latent variable approaches for ‘Omics studies](http://mixomics.org/mixdiablo/)  
12.	mixMC	[Multivariate Statistical Framework to Gain Insight into Microbial Communities](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0160169)  
13.	MMinte	[Methodology for the large-scale assessment of microbial metabolic interactions (MMinte) from 16S rDNA data](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1230-3)  
14.	pathostat	[Statistical Microbiome Analysis on metagenomics results from sequencing data samples](https://bioconductor.org/packages/release/bioc/html/PathoStat.html)  
15.	phylofactor	[Phylogenetic factorization of compositional data](https://peerj.com/articles/2969/)  
16.	phylogeo	[Geographic analysis and visualization of microbiome data](https://www.ncbi.nlm.nih.gov/pubmed/25913208)  
17.	Phyloseq	[Import, share, and analyze microbiome census data using R](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0061217)  
18.	qiimer	[R tools compliment qiime](https://github.com/kylebittinger/qiimer)  
19.	RAM	[R for Amplicon-Sequencing-Based Microbial-Ecology](https://cran.r-project.org/web/packages/RAM/RAM.pdf)  
20.	ShinyPhyloseq	[Web-tool with user interface for Phyloseq](http://joey711.github.io/shiny-phyloseq/)  
21.	SigTree	[Identify and Visualize Significantly Responsive Branches in a Phylogenetic Tree](http://www.sciencedirect.com/science/article/pii/S2001037017300132)  
22.	SPIEC-EASI	[Sparse and Compositionally Robust Inference of Microbial Ecological Networks](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004226)  
23.	structSSI	[Simultaneous and Selective Inference for Grouped or Hierarchically Structured Data](https://www.jstatsoft.org/article/view/v059i13)  
24.	Tax4Fun	[Predicting functional profiles from metagenomic 16S rRNA gene data](https://www.ncbi.nlm.nih.gov/pubmed/25957349)  
25.	taxize	[Taxonomic Information from Around the Web](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3901538/)    
26. labdsv	[Ordination and Multivariate Analysis for Ecology](https://cran.r-project.org/web/packages/labdsv/labdsv.pdf)  
27. Vegan	[R package for community ecologists](https://github.com/vegandevs/vegan)  
28.	igraph	[Network Analysis and Visualization in R](http://igraph.org/r/)
29.	MicrobiomeHD	[A standardized database of human gut microbiome studies in health and disease *Case-Control*](http://www.biorxiv.org/content/early/2017/05/08/134031)   
30.	Rhea	[A pipeline with modular R scripts](https://peerj.com/articles/2836/)  
31. microbiomeutilities [Extending and supporting package based on microbiome and phyloseq R package](https://microsud.github.io/microbiomeutilities/)   
32. breakaway [Species Richness Estimation and Modeling](https://github.com/search?q=breakaway)   
33. corncob [Modeling microbial abundances and dysbiosis with beta-binomial regression](https://arxiv.org/abs/1902.02776)   
34. MicrobiomeR [MicrobiomeR: An R Package for Simplified and Standardized Microbiome Analysis Workflows](https://microbiomer.vallenderlab.org/)   
35. powmic [Power assessment in microbiome case-control studies](https://github.com/lichen-lab/powmic)  
36. yingtools2 [Tools and functions for working with clinical and microbiome data](https://github.com/ying14/yingtools2)  
37. animalcules [R shiny app for interactive microbiome analysis](https://compbiomed.github.io/animalcules-docs/)    
38. biome-shiny [GUI for microbiome visualization, based on the shiny package "microbiome"](https://github.com/Biodata-PT/Biome-Shiny)  
39. MelonnPan [Model-based Genomically Informed High-dimensional Predictor of Microbial Community Metabolic Profiles](https://github.com/biobakery/melonnpan)  
40. MaAsLin2 [MaAsLin2: Microbiome Multivariate Association with Linear Models](https://github.com/biobakery/Maaslin2)  
41. mbtools [Collection of helpers that we use to analyze microbiome data](https://gibbons-lab.github.io/mbtools)   
42. ANCOM [R scripts for Analysis of Composition of Microbiomes (ANCOM)](https://github.com/FrederickHuangLin/ANCOM)  
43. MetaLonDA [METAgenomic LONgitudinal Differential Abundance method](https://github.com/aametwally/MetaLonDA)  
44. dacomp [Testing for Differential Abundance in Compositional Counts Data, with Application to Microbiome Studies](https://github.com/barakbri/dacomp)  
45. BEEM [BEEM: Estimating Lotka-Volterra models from time-course microbiome sequencing data](https://github.com/lch14forever/BEEM)  
46. metamicrobiomeR [Analysis of Microbiome Relative Abundance Data using Zero Inflated Beta GAMLSS and Meta-Analysis Across Studies using Random Effects Model](https://github.com/cran/metamicrobiomeR)   
47. GLMMMiRKAT [A distance-based kernel association test based on the generalized linear mixed model](https://github.com/hk1785/GLMM-MiRKAT)  
48. MDPbiome [MDPbiome: microbiome engineering through prescriptive perturbations](https://academic.oup.com/bioinformatics/article/34/17/i838/5093255)  
49. bootLong [The Block Bootstrap Method for Longitudinal Microbiome Data](https://pratheepaj.github.io/bootLong/)   
50. OMiSA [Optimal Microbiome-based Survival Analysis (OMiSA)](https://github.com/hk1785/OMiSA)    
51. DECIPHER [Using DECIPHER v2.0 to Analyze Big Biological Sequence Data in R](http://www2.decipher.codes/index.html)  
52. DECIPHER/IIDTAXA [IDTAXA: a novel approach for accurate taxonomic classification of microbiome sequences](doi:10.1186/s40168-018-0521-5)  
53. curatedMetagenomicData [Accessible, curated metagenomic data through ExperimentHub](https://waldronlab.io/curatedMetagenomicData/)  
54. themetagenomics [Exploring Thematic Structure and Predicted Functionality of 16S rRNA Amplicon Data](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0219235)   
55. MDiNE [MDiNE: a model to estimate differential co-occurrence networks in microbiome studies](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btz824/5614428)  
56. StructFDR [False Discovery Rate Control Incorporating Phylogenetic Tree Increases Detection Power in Microbiome-Wide Multiple Testing](https://academic.oup.com/bioinformatics/article/33/18/2873/3824757)   
57. metamicrobiomeR [metamicrobiomeR: An R Package for Analysis of Microbiome Relative Abundance Data Using Zero-Inflated Beta GAMLSS and Meta-Analysis Across Studies Using Random Effects Models](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2744-2)  
58. Pldist [Pldist: Ecological Dissimilarities for Paired and Longitudinal Microbiome Association Analysis](https://academic.oup.com/bioinformatics/article-abstract/35/19/3567/5341424?redirectedFrom=fulltext)    
59. BDMMA [Batch Effects Correction for Microbiome Data With Dirichlet-multinomial Regression](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/bty729)  
60. RCM [A unified framework for unconstrained and constrained ordination of microbiome read count data](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0205474)   
61. decontam [Simple Statistical Identification and Removal of Contaminant Sequences in Marker-Gene and Metagenomics Data](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0605-2)  
62. ZIBBSeqDiscovery [A Zero-inflated Beta-binomial Model for Microbiome Data Analysis](https://onlinelibrary.wiley.com/doi/abs/10.1002/sta4.185)  
63. massMap [A Two-Stage Microbial Association Mapping Framework With Advanced FDR Control](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0517-1)  
64. SplinectomeR [SplinectomeR Enables Group Comparisons in Longitudinal Microbiome Studies](https://www.frontiersin.org/articles/10.3389/fmicb.2018.00785/full)  
65. DMBC [A Dirichlet-Multinomial Bayes Classifier for Disease Diagnosis With Microbial Compositions](https://msphere.asm.org/content/2/6/e00536-17)  
66. MicrobiomeDDA [An Omnibus Test for Differential Distribution Analysis of Microbiome Sequencing Data](https://academic.oup.com/bioinformatics/article/34/4/643/4470360)  
67. NMIT [Microbial Interdependence Association Test--a Non-parametric Microbial Interdependence Test](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5696116/#!po=75.0000)  
68. SparseMCMM [Estimating and testing the microbial causal mediation effect with the high-dimensional and compositional microbiome data (SparseMCMM)](https://sites.google.com/site/huilinli09/software)  
69. MTA [Microbial trend analysis (MTA) for common dynamic trend, group comparison and classification in longitudinal microbiome study](https://github.com/chanw0/MTA)  
70. miLineage [A General Framework for Association Analysis of Microbial Communities on a Taxonomic Tree](https://medschool.vanderbilt.edu/tang-lab/software/miLineage)   
71. zeroSum [Reference Point Insensitive Molecular Data Analysis](https://academic.oup.com/bioinformatics/article/33/2/219/2928229)  
72. MedTest [A Distance-Based Approach for Testing the Mediation Effect of the Human Microbiome](https://academic.oup.com/bioinformatics/article/34/11/1875/4810437)     
73. qgraph [Graph Plotting Methods, Psychometric Data Visualization and
Graphical Model Estimation](https://cran.r-project.org/web/packages/qgraph/qgraph.pdf)  
74. Adaptive gPCA [A method for structured dimensionality reduction](https://arxiv.org/abs/1702.00501)    

-------------------------------------------------------   

### Other tools    
1. ggplot2 [An implementation of the Grammar of Graphics in R](https://ggplot2.tidyverse.org/)  
    * Widely used package for data visualization  
2. ggvegan [ggplot-based versions of the plots produced by the vegan package](https://github.com/gavinsimpson/ggvegan)  
    * Convert base plots of vegan to ggplot.     
3. ggord [A simple package for creating ordination plots with ggplot2](https://fawda123.github.io/ggord/)  
    * Alternative to ggvegan    
4. cowplot [cowplot: Streamlined Plot Theme and Plot Annotations for ggplot2](https://wilkelab.org/cowplot/)  
    * Widely used package for combining multiple plots     
4. ggridges [Ridgeline plots in ggplot2](https://wilkelab.org/ggridges)   
5. ggtext [Improved text rendering support for ggplot2](https://wilkelab.org/ggtext/)   
    * More power in controlling annotations in plots (e.g. italicize taxa names in plots)  
6. patchwork [The Composer of ggplots](https://patchwork.data-imaginist.com)  
    * Combining multiple plots made easy   
7. ggpubr [Extension of ggplot2 based data visualization](http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/)  
    * Publication ready plots   
8. ggraph [Grammar of Graph Graphics](https://ggraph.data-imaginist.com)   
    * Network graphs using ggplot2    
9. gganimate [A Grammar of Animated Graphics](https://gganimate.com)  
    * Animate ggplot2 (Useful for presenting time-series dynamics of microbial communities)  
10. ggforce [Accelerating ggplot2](https://ggforce.data-imaginist.com)  
    * Zoom specific regions of the plots   
11. factoextra [Extract and Visualize the Results of Multivariate Data Analyses](https://rpkgs.datanovia.com/factoextra/index.html)  
    * Powerful package for multivvariate data analysis  
12. ggcorrplot [Visualization of a correlation matrix using ggplot2](https://rpkgs.datanovia.com/ggcorrplot/)   
13. tidyverse [R packages for data science](https://www.tidyverse.org/)  
    * Universe of several useful R packages for data handling, analysis and vidualization    
14. Extensions of ggplot [Gallary of numerous data visualistion R pacakges](https://www.ggplot2-exts.org/gallery/)  

-------------------------------------------------------   

### Proteomics resources<sup>*<sup>  

1. RforProteomics [Using R for proteomics data analysis](https://bioconductor.org/packages/release/data/experiment/vignettes/RforProteomics/inst/doc/RforProteomics.html)  
2. RforProteomics [Visualisation of proteomics data using R and Bioconductor](https://onlinelibrary.wiley.com/doi/full/10.1002/pmic.201400392)   
2. proteomics [proteomics: Mass spectrometry and proteomics data analysis](http://master.bioconductor.org/packages/release/workflows/vignettes/proteomics/inst/doc/proteomics.html)  

-------------------------------------------------------   

### RNAseq resources<sup>*<sup>    

1. RNA-seq analysis in R [Workflow by Shulin Cao](https://rstudio-pubs-static.s3.amazonaws.com/462299_a9bc385f89b94b0aa95de0f3b7040b04.html)  
2. RNA-seq workflow [RNA-seq workflow: gene-level exploratory analysis and differential expression](https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html)   

*Note: These are not focused towards microbiome data. These are listed as a reference point for beginners. If you have or know of workflows tools specific for microbiome data please let us know and we can add them here!  

-------------------------------------------------------   

### Useful resources are provided by:  
1. [Ben J. Callahan and Colleagues: Bioconductor Workflow for Microbiome Data Analysis: from raw reads to community analyses](https://f1000research.com/articles/5-1492/v2).   
2. [Comeau AM and Colleagues: Microbiome Helper: a Custom and Streamlined Workflow for Microbiome Research](http://msystems.asm.org/content/2/1/e00127-16)  
3. [Shetty SA, Lahti L., et al: Tutorial from microbiome data analysis spring school 2018, Wageningen University and Research](https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/introduction.html)     
4. [Holmes S, Huber W.: Modern statistics for modern biology. Cambridge University Press; 2018 Nov 30.](http://web.stanford.edu/class/bios221/book/)   

Note:  
A good practise is to use Rmarkdown for documenting your results and sharing with your collaborators and supervisors. For more information click here [RStudio](https://www.youtube.com/watch?v=cWJzjHh_3kk&t=337s) and  
[RStudio Overview](https://www.youtube.com/watch?v=n3uue28FD0w)  

-------------------------------------------------------   

[View this webiste repository on GitHub](https://github.com/microsud/Tools-Microbiome-Anlaysis)  
[Follow me on Twitter](https://twitter.com/gutmicrobe)  
[googlescholar](https://scholar.google.nl/citations?hl=en&user=Vahc6LUAAAAJ&view_op=list_works&sortby=pubdate)  
[ORCID ID: 0000-0001-7280-9915](http://orcid.org/0000-0001-7280-9915)   

-------------------------------------------------------   

### References:
1. Callahan, B. J., McMurdie, P. J. & Holmes, S. P. (2017). Exact sequence variants should replace operational taxonomic units in marker gene data analysis. bioRxiv, 113597.  
2. Callahan, B. J., McMurdie, P. J., Rosen, M. J., Han, A. W., Johnson, A. J. A. & Holmes, S. P. (2016). DADA2: high-resolution sample inference from Illumina amplicon data. Nature methods 13, 581-583.  
3. Caporaso, J. G., Kuczynski, J., Stombaugh, J., Bittinger, K., Bushman, F. D., Costello, E. K., Fierer, N., Peña, A. G., Goodrich, J. K. & Gordon, J. I. (2010). QIIME allows analysis of high-throughput community sequencing data. Nature methods 7, 335-336.  
4. Schloss, P. D., Westcott, S. L., Ryabin, T., Hall, J. R., Hartmann, M., Hollister, E. B., Lesniewski, R. A., Oakley, B. B., Parks, D. H. & Robinson, C. J. (2009). Introducing mothur: open-source, platform-independent, community-supported software for describing and comparing microbial communities. Applied and environmental microbiology 75, 7537-7541.  
5. Team, R. C. (2000). R language definition. Vienna, Austria: R foundation for statistical computing.  
-------------------------------------------------------   

### Was this website/resource useful for you? Then please share it with others too!  
You can cite this resource as:  
Shetty, Sudarshan A., and Leo Lahti. [Microbiome data science. Journal of biosciences 44, no. 5 (2019): 115](https://link.springer.com/article/10.1007%2Fs12038-019-9930-2).    
[*Pre-print*](https://openresearchlabs.github.io/publications/papers/2019-Shetty-MDS.pdf)   
Zendo: [![DOI](https://zenodo.org/badge/99698135.svg)](https://zenodo.org/badge/latestdoi/99698135)   

-------------------------------------------------------

### TODO 
Any help is welcome  
* Structure the list according to categories  
  * General purpose    
  * Visualization  
  * Snapshot/cross-sectional stats  
  * Time series/Longitudinal stats  
  * Integrative -Omics  
* Include metagenomics/metabolomics  
* Include more general microbiology oriented R packages/tools  
* List of 'good' research paper reproducible repositories   
* and so on .....
-------------------------------------------------------   
Google doc [link](https://docs.google.com/spreadsheets/d/1am-UyDVBGDOgm6jVQ5FDXxmg24iriHqeBeul14HRb1g/edit?usp=sharing)  

-------------------------------------------------------   
