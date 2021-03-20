
[![DOI](https://zenodo.org/badge/99698135.svg)](https://zenodo.org/badge/latestdoi/99698135)        

## A list of R environment based tools for microbiome data exploration, statistical analysis and visualization  
  
As a beginner, the entire process from sample collection to analysis for sequencing data is a daunting task. More specifically, the downstream processing of raw reads is the most time consuming and mentally draining stage. It is vital to understand the basic concepts in microbial ecology and then to use various tools at disposal to address specific research questions. Thankfully, several young researchers supported by their experienced principal investigators/supervisors are working on creating various tools for analysis and interpretation of microbial community data. A major achievement of the scientific community is the open science initiative which has led to sharing of knowledge worldwide. For microbial community analysis, several tools have been created in R, a free to use (GNU General Public License) programming language(Team, 2000). The power of R lies in its ease of working with individuals lacking programming skills and easy sharing of analysis scripts codes and packages aiding reproducibility. Using tools such as QIIME (the newer QIIME2) (Caporaso, Kuczynski, Stombaugh et al., 2010), Mothur (Schloss, Westcott, Ryabin et al., 2009), DADA2 (Callahan, McMurdie, Rosen et al., 2016) one can get from raw reads to species × samples table (OTU or ASVs amplicon sequence variants as suggested recently (Callahan, McMurdie & Holmes, 2017)). In this post, numerous resources that can be helpful for analysis of microbiome data are listed. This list may not have all the packages as this tool development space is ever growing. Feel free to add those packages or links to web tutorials related to microbiome data. These are mostly for improving statistical analysis and visualization. These tools provide convenient options for data analysis and include several steps where the user has to make decisions. The work by [McMurdie PJ, Holmes S](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003531), [Weiss S](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0237-y) and [Tsilimigras M.C. and Fodor A.A](http://www.sciencedirect.com/science/article/pii/S1047279716300722) are useful resources to understand the data common to microbiome census. It can be tricky and frustrating in the beginning but patience and perseverance will be fruitful at the end (personal experience).

-------------------------------------------------------   

### Tools     

> A   

Adaptive gPCA [A method for structured dimensionality reduction](https://arxiv.org/abs/1702.00501)     
Ampvis2 [Tools for visualising amplicon sequencing data](https://madsalbertsen.github.io/ampvis2/)   
ANCOM [R scripts for Analysis of Composition of Microbiomes (ANCOM)](https://github.com/FrederickHuangLin/ANCOM)   
animalcules [R shiny app for interactive microbiome analysis](https://compbiomed.github.io/animalcules-docs/)   

> B    

BDMMA [Batch Effects Correction for Microbiome Data With Dirichlet-multinomial Regression](https://academic.oup.com/bioinformatics/article-abstract/35/5/807/5078473)   
BEEM [BEEM: Estimating Lotka-Volterra models from time-course microbiome sequencing data](https://github.com/lch14forever/BEEM)   
biome-shiny [GUI for microbiome visualization, based on the shiny package "microbiome"](https://github.com/Biodata-PT/Biome-Shiny)   
bootLong [The Block Bootstrap Method for Longitudinal Microbiome Data](https://pratheepaj.github.io/bootLong/)   
breakaway [Species Richness Estimation and Modeling](https://github.com/search?q=breakaway)   

> C   

CCREPE [Compositionality Corrected by PErmutation and REnormalization](http://bioconductor.org/packages/release/bioc/html/ccrepe.html)   
corncob [Modeling microbial abundances and dysbiosis with beta-binomial regression](https://arxiv.org/abs/1902.02776)  
curatedMetagenomicData [Accessible, curated metagenomic data through ExperimentHub](https://waldronlab.io/curatedMetagenomicData/)    

> D     

dacomp [Testing for Differential Abundance in Compositional Counts Data, with Application to Microbiome Studies](https://github.com/barakbri/dacomp)   
DADA2 [Divisive Amplicon Denoising Algorithm](https://www.nature.com/nmeth/journal/v13/n7/full/nmeth.3869.html)    
DECIPHER [Using DECIPHER v2.0 to Analyze Big Biological Sequence Data in R](http://www2.decipher.codes/index.html)   
DECIPHER/IIDTAXA [IDTAXA: a novel approach for accurate taxonomic classification of microbiome sequences](doi:10.1186/s40168-018-0521-5)   
decontam [Simple Statistical Identification and Removal of Contaminant Sequences in Marker-Gene and Metagenomics Data](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0605-2)   
DESeq2 [Differential expression analysis for sequence count data](https://www.bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-mixture-models.html)   
DMBC [A Dirichlet-Multinomial Bayes Classifier for Disease Diagnosis With Microbial Compositions](https://msphere.asm.org/content/2/6/e00536-17)    

> E   

edgeR [empirical analysis of DGE in R](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2796818/)    

> G     

GLMMMiRKAT [A distance-based kernel association test based on the generalized linear mixed model](https://github.com/hk1785/GLMM-MiRKAT)   

> I   

igraph [Network Analysis and Visualization in R](http://igraph.org/r/)    

> L    

labdsv [Ordination and Multivariate Analysis for Ecology](https://cran.r-project.org/web/packages/labdsv/labdsv.pdf)  

> M   

MaAsLin2 [MaAsLin2: Microbiome Multivariate Association with Linear Models](https://github.com/biobakery/Maaslin2)   
mare [Microbiota Analysis in R Easily](https://github.com/katrikorpela/mare)   
massMap [A Two-Stage Microbial Association Mapping Framework With Advanced FDR Control](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0517-1)  
mbtools [Collection of R tools to analyze microbiome data](https://gibbons-lab.github.io/mbtools)   
MDiNE [MDiNE: a model to estimate differential co-occurrence networks in microbiome studies](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btz824/5614428)   
MDPbiome [MDPbiome: microbiome engineering through prescriptive perturbations](https://academic.oup.com/bioinformatics/article/34/17/i838/5093255)    
microDecon [An R package for removing contamination from metabarcoding (e.g., microbiome) datasets post-sequencing](https://github.com/donaldtmcknight/microDecon)  
microbiotaPair [An R Package for Simplified Microbiome Analysis Workflows on Intervention Study](https://github.com/rusher321/microbiotaPair)    
MedTest [A Distance-Based Approach for Testing the Mediation Effect of the Human Microbiome](https://academic.oup.com/bioinformatics/article/34/11/1875/4810437)    
MegaR [An Interactive R Package for Metagenomic Sample Classification and Disease Prediction using Microbiota and Machine Learning](https://github.com/BioHPC/MegaR)     
MelonnPan [Model-based Genomically Informed High-dimensional Predictor of Microbial Community Metabolic Profiles](https://github.com/biobakery/melonnpan)  
Metacoder [An R package for visualization and manipulation of community taxonomic diversity data](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005404)   
metagenomeSeq [Differential abundance analysis for microbial marker-gene surveys](http://www.nature.com/nmeth/journal/v10/n12/full/nmeth.2658.html)   
MetaLonDA [METAgenomic LONgitudinal Differential Abundance method](https://github.com/aametwally/MetaLonDA)   
metamicrobiomeR [Analysis of Microbiome Relative Abundance Data using Zero Inflated Beta GAMLSS and Meta-Analysis Across Studies using Random Effects Model](https://github.com/cran/metamicrobiomeR)   
microbiome R package  [Tools for microbiome analysis in R](https://github.com/microbiome/microbiome)    
MicrobiomeDDA [An Omnibus Test for Differential Distribution Analysis of Microbiome Sequencing Data](https://academic.oup.com/bioinformatics/article/34/4/643/4470360)   
MicrobiomeHD [A standardized database of human gut microbiome studies in health and disease *Case-Control*](http://www.biorxiv.org/content/early/2017/05/08/134031)   
microbiomeMarker [R package for microbiome biomarker discovery](https://github.com/yiluheihei/microbiomeMarker)   
MicrobiomeR [MicrobiomeR: An R Package for Simplified and Standardized Microbiome Analysis Workflows](https://microbiomer.vallenderlab.org/)   
microbiomeutilities [Extending and supporting package based on microbiome and phyloseq R package](https://microsud.github.io/microbiomeutilities/)   
MicrobiotaProcess [MicrobiotaProcess: an R package for analysis, visualization and biomarker discovery of microbiome](https://github.com/YuLab-SMU/MicrobiotaProcess)   
miLineage [A General Framework for Association Analysis of Microbial Communities on a Taxonomic Tree](https://medschool.vanderbilt.edu/tang-lab/software/miLineage)   
MINT [Multivariate INTegrative method](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1553-8)   
mixDIABLO [Data Integration Analysis for Biomarker discovery using Latent variable approaches for ‘Omics studies](http://mixomics.org/mixdiablo/)   
mixMC [Multivariate Statistical Framework to Gain Insight into Microbial Communities](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0160169)   
MMinte [Methodology for the large-scale assessment of microbial metabolic interactions (MMinte) from 16S rDNA data](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1230-3)   
MTA [Microbial trend analysis (MTA) for common dynamic trend, group comparison and classification in longitudinal microbiome study](https://github.com/chanw0/MTA)    
MicrobiomeExplor [An R package for the analysis and visualization of microbial communities](https://doi.org/10.1093/bioinformatics/btaa838)   
MicroBVS [Dirichlet-tree multinomial regression models with Bayesian variable selection - an R package](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03640-0)    
MicroNiche [An R package for assessing microbial niche breadth and overlap from amplicon sequencing data](https://academic.oup.com/femsec/article-abstract/96/8/fiaa131/5863182?redirectedFrom=fulltext)   

> N    
 
NMIT [Microbial Interdependence Association Test--a Non-parametric Microbial Interdependence Test](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5696116/#!po=75.0000)    
NetCoMi [Network construction and comparison for microbiome data in R](https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bbaa290/6017455)   
NBZIMM [Negative binomial and zero-inflated mixed models, with application to microbiome/metagenomics data analysis](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03803-z)    

> O    

OMiSA [Optimal Microbiome-based Survival Analysis (OMiSA)](https://github.com/hk1785/OMiSA)   

> P   

pathostat [Statistical Microbiome Analysis on metagenomics results from sequencing data samples](https://bioconductor.org/packages/release/bioc/html/PathoStat.html)   
phylofactor [Phylogenetic factorization of compositional data](https://peerj.com/articles/2969/)   
phylogeo [Geographic analysis and visualization of microbiome data](https://www.ncbi.nlm.nih.gov/pubmed/25913208)   
Phyloseq [Import, share, and analyze microbiome census data using R](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0061217)   
Pldist [Pldist: Ecological Dissimilarities for Paired and Longitudinal Microbiome Association Analysis](https://academic.oup.com/bioinformatics/article-abstract/35/19/3567/5341424?redirectedFrom=fulltext)   
powmic [Power assessment in microbiome case-control studies](https://github.com/lichen-lab/powmic)   

> Q   

qgraph [Graph Plotting Methods, Psychometric Data Visualization and Graphical Model Estimation](https://cran.r-project.org/web/packages/qgraph/qgraph.pdf)   
qiime2R [Importing QIIME2 artifacts and associated data into R sessions](https://github.com/jbisanz/qiime2R)   
qiimer [R tools compliment qiime](https://github.com/kylebittinger/qiimer)    

> R    

RAM [R for Amplicon-Sequencing-Based Microbial-Ecology](https://cran.r-project.org/web/packages/RAM/RAM.pdf)   
RCM [A unified framework for unconstrained and constrained ordination of microbiome read count data](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0205474)   
RevEcoR [Reverse Ecology analysis in R](https://github.com/yiluheihei/RevEcoR)   
Rhea [A pipeline with modular R scripts](https://peerj.com/articles/2836/)   

> S   

selbal [Balances: a New Perspective for Microbiome Analysis](https://msystems.asm.org/content/3/4/e00053-18)  
ShinyPhyloseq [Web-tool with user interface for Phyloseq](http://joey711.github.io/shiny-phyloseq/)     
SIAMCAT [Statistical Inference of Associations between Microbial Communities And host phenoTypes](https://bioconductor.org/packages/release/bioc/html/SIAMCAT.html)   
SigTree [Identify and Visualize Significantly Responsive Branches in a Phylogenetic Tree](http://www.sciencedirect.com/science/article/pii/S2001037017300132)    
SparseMCMM [Estimating and testing the microbial causal mediation effect with the high-dimensional and compositional microbiome data (SparseMCMM)](https://sites.google.com/site/huilinli09/software)   
SPIEC-EASI [Sparse and Compositionally Robust Inference of Microbial Ecological Networks](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004226)   
SplinectomeR [SplinectomeR Enables Group Comparisons in Longitudinal Microbiome Studies](https://www.frontiersin.org/articles/10.3389/fmicb.2018.00785/full)   
StructFDR [False Discovery Rate Control Incorporating Phylogenetic Tree Increases Detection Power in Microbiome-Wide Multiple Testing](https://academic.oup.com/bioinformatics/article/33/18/2873/3824757)   
structSSI [Simultaneous and Selective Inference for Grouped or Hierarchically Structured Data](https://www.jstatsoft.org/article/view/v059i13)   

> T   

Tax4Fun [Predicting functional profiles from metagenomic 16S rRNA gene data](https://www.ncbi.nlm.nih.gov/pubmed/25957349)    
taxize [Taxonomic Information from Around the Web](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3901538/)    
TreeSummarizedExperiment [An extension class of SummarizedExperiment to allow hierarchical structure on the row or column dimension](https://github.com/fionarhuang/TreeSummarizedExperiment)   
themetagenomics [Exploring Thematic Structure and Predicted Functionality of 16S rRNA Amplicon Data](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0219235)    
tidyMicro [A pipeline for microbiome data analysis and visualization using the tidyverse in R](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-03967-2)  

> V   

vegan  [R package for community ecologists](https://github.com/vegandevs/vegan)   

> Y   

yingtools2 [Tools and functions for working with clinical and microbiome data](https://github.com/ying14/yingtools2)   

> Z    

zeroSum [Reference Point Insensitive Molecular Data Analysis](https://academic.oup.com/bioinformatics/article/33/2/219/2928229)   
ZIBBSeqDiscovery [A Zero-inflated Beta-binomial Model for Microbiome Data Analysis](https://onlinelibrary.wiley.com/doi/abs/10.1002/sta4.185)   

-------------------------------------------------------    

### Microbiome data sets  

1. HMP16SData [Efficient Access to the Human Microbiome Project through Bioconductor](https://bioconductor.org/packages/release/data/experiment/html/HMP16SData.html)   
2. HMP2Data [16s rRNA sequencing data from the Human Microbiome Project 2](https://bioconductor.org/packages/release/data/experiment/html/HMP2Data.html)  
3. microbiomeDataSets [Experiment Hub based microbiome datasets](https://bioconductor.org/packages/devel/data/experiment/html/microbiomeDataSets.html)   
4. curatedMetagenomicData [Accessible, curated metagenomic data through ExperimentHub](https://bioconductor.org/packages/release/data/experiment/html/curatedMetagenomicData.html)   

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
6. ggpubr [Extension of ggplot2 based data visualization](http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/)  
    * Publication ready plots   
7. ggraph [Grammar of Graph Graphics](https://ggraph.data-imaginist.com)   
    * Network graphs using ggplot2    
8. gganimate [A Grammar of Animated Graphics](https://gganimate.com)  
    * Animate ggplot2 (Useful for presenting time-series dynamics of microbial communities)  
9. ggforce [Accelerating ggplot2](https://ggforce.data-imaginist.com)  
    * Zoom specific regions of the plots   
10. factoextra [Extract and Visualize the Results of Multivariate Data Analyses](https://rpkgs.datanovia.com/factoextra/index.html)  
    * Powerful package for multivvariate data analysis  
11. ggcorrplot [Visualization of a correlation matrix using ggplot2](https://rpkgs.datanovia.com/ggcorrplot/)   
12. tidyverse [R packages for data science](https://www.tidyverse.org/)  
    * Universe of several useful R packages for data handling, analysis and visualization    
13. Extensions of ggplot [Gallary of numerous data visualistion R pacakges](https://exts.ggplot2.tidyverse.org/gallery/)   
14. ggtree [Visualization and annotation of phylogenetic trees (in R) ](https://github.com/YuLab-SMU/ggtree)   
15. patchwork [The Composer of ggplots](https://patchwork.data-imaginist.com)  
    * Combining multiple plots made easy    
16. pheatmap [Pretty Heatmaps](https://cran.r-project.org/web/packages/pheatmap/pheatmap.pdf)    

-------------------------------------------------------   

### Proteomics resources    

1. <sup>*</sup>RforProteomics [Using R for proteomics data analysis](https://bioconductor.org/packages/release/data/experiment/vignettes/RforProteomics/inst/doc/RforProteomics.html)  
2. <sup>*</sup>RforProteomics [Visualisation of proteomics data using R and Bioconductor](https://onlinelibrary.wiley.com/doi/full/10.1002/pmic.201400392)   
3. <sup>*</sup>proteomics [proteomics: Mass spectrometry and proteomics data analysis](http://master.bioconductor.org/packages/release/workflows/vignettes/proteomics/inst/doc/proteomics.html)    
4. [Introduction to analysing microbial proteomics data in R](https://microsud.github.io/Bacterial-Proteomics-in-R/tutorial_main_doc.html)   

-------------------------------------------------------   

### RNAseq resources<sup>*<sup>    

1. RNA-seq analysis in R [Workflow by Shulin Cao](https://rstudio-pubs-static.s3.amazonaws.com/462299_a9bc385f89b94b0aa95de0f3b7040b04.html)  
2. RNA-seq workflow [RNA-seq workflow: gene-level exploratory analysis and differential expression](https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html)   

*Note: These are not focused towards microbiome data. These are listed as a reference point for beginners. If you have or know of workflows tools specific for microbiome data please let us know and we can add them here!  

-------------------------------------------------------   

### Useful resources are provided by:  
1. [Ben J. Callahan and Colleagues: Bioconductor Workflow for Microbiome Data Analysis: from raw reads to community analyses](https://f1000research.com/articles/5-1492/v2).   
2. [Comeau AM and Colleagues: Microbiome Helper: a Custom and Streamlined Workflow for Microbiome Research](http://msystems.asm.org/content/2/1/e00127-16)   
3. Schloss, P. D: [The Riffomonas Reproducible Research Tutorial Series](http://www.riffomonas.org/reproducible_research/)  
4. [Shetty SA, Lahti L., et al: Tutorial from microbiome data analysis spring school 2018, Wageningen University and Research](https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/introduction.html)     
5. [Holmes S, Huber W.: Modern statistics for modern biology. Cambridge University Press; 2018 Nov 30.](http://web.stanford.edu/class/bios221/book/)    
6. [Xu S, Yu G.: Workshop of microbiome dataset analysis using MicrobiotaProcess](https://yulab-smu.top/MicrobiotaProcessWorkshop/articles/MicrobiotaProcessWorkshop.html)   
7. [Antoni Susin, Yiwen Wang, Kim-Anh Lê Cao, M.Luz Calle.: Variable selection in microbiome compositional data analysis: tutorial](https://malucalle.github.io/Microbiome-Variable-Selection/)    

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
6. Salter SJ, Cox MJ, Turek EM, Calus ST, Cookson WO, Moffatt MF, Turner P, Parkhill J, Loman NJ, Walker AW. Reagent and laboratory contamination can critically impact sequence-based microbiome analyses. BMC biology. 2014 Dec 1;12(1):87.    
7. Karstens L, Asquith M, Davin S, Fair D, Gregory WT, Wolfe AJ, Braun J, McWeeney S. Controlling for contaminants in low-biomass 16S rRNA gene sequencing experiments. MSystems. 2019 Aug 27;4(4).   

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

* and so on .....  

-------------------------------------------------------   

### Was this website/resource useful for you? Then please share it with others too!  
You can cite this resource as:  
Shetty, Sudarshan A., and Leo Lahti. [Microbiome data science. Journal of biosciences 44, no. 5 (2019): 115](https://link.springer.com/article/10.1007%2Fs12038-019-9930-2).    
[*Pre-print*](https://openresearchlabs.github.io/publications/papers/2019-Shetty-MDS.pdf)   
Zendo: [![DOI](https://zenodo.org/badge/99698135.svg)](https://zenodo.org/badge/latestdoi/99698135)    


-------------------------------------------------------    

### Do you want to join us and contribute to the open source microbiome community?   
We are working on an open source project developing R/Bioc methods, benchmarking data, and educational material for microbiome research based on the `SummarizedExperiment` class and its derivatives.  
More information please check our website [R/Bioc microbiome ecosystem with SummarizedExperiment](https://microbiome.github.io/)   

-------------------------------------------------------    

### Acknowledgments        

* The Academy of Finland (Grants 295741; 307127)     
* The SIAM Gravitation Grant 024.002.002   
* UNLOCK project of the Netherlands Organization for Scientific Research (NWO)   
* University of Turku, Department of Mathematics and Statistics   
* Molecular Ecology group, Laboratory of Microbiology, Wageningen University, The Netherlands   
* Medische Microbiologie en Infectiepreventie, University Medical Center Groningen, The Netherlands   

-------------------------------------------------------     

### Updates   

* Updates that are added [See Changes log](CHANGES.html)   
  

-------------------------------------------------------   

