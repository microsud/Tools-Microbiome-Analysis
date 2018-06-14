
[![DOI](https://zenodo.org/badge/99698135.svg)](https://zenodo.org/badge/latestdoi/99698135)  

## A list of R environment based tools for 16S rRNA gene data exploration, statistical analysis and visualization  
  
As a beginner, the entire process from sample collection to analysis for sequencing data is a daunting task. More specifically, the downstream processing of raw reads is the most time consuming and mentally draining stage. It is vital to understand the basic concepts in microbial ecology and then to use various tools at disposal to address specific research questions. Thankfully, several young researchers supported by their experienced principal investigators/supervisors are working on creating various tools for analysis and interpretation of microbial community data. A major achievement of the scientific community is the open science initiative which has led to sharing of knowledge worldwide. For microbial community analysis, several tools have been created in R, a free to use (GNU General Public License) programming language(Team, 2000). The power of R lies in its ease of working with individuals lacking programming skills and easy sharing of analysis scripts codes and packages aiding reproducibility. Using tools such as QIIME (the newer QIIME2) (Caporaso, Kuczynski, Stombaugh et al., 2010), Mothur (Schloss, Westcott, Ryabin et al., 2009), DADA2 (Callahan, McMurdie, Rosen et al., 2016) one can get from raw reads to species × samples table (OTU or ASVs amplicon sequence variants as suggested recently (Callahan, McMurdie & Holmes, 2017)). In this post, numerous resources that can be helpful for analysis of microbiome data are listed. This list may not have all the packages as this tool development space is ever growing. Feel free to add those packages or links to web tutorials related to microbiome data, there is a [google docs excel sheet at this link for a list of tools](https://docs.google.com/spreadsheets/d/1am-UyDVBGDOgm6jVQ5FDXxmg24iriHqeBeul14HRb1g/edit?usp=sharing) which can be edited to include more tools. These are mostly for improving statistical analysis and visualisation. These tools provide convenient options for data analysis and include several steps where the user has to make decisions. The work by [McMurdie PJ, Holmes S](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003531), [Weiss S](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0237-y) and [Tsilimigras M.C. and Fodor A.A](http://www.sciencedirect.com/science/article/pii/S1047279716300722) are useful resources to understand the data common to microbiome census. It can be tricky and frustrating in the beginning but patience and perseverance will be fruitful at the end (personal experience).

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

Google doc [link](https://docs.google.com/spreadsheets/d/1am-UyDVBGDOgm6jVQ5FDXxmg24iriHqeBeul14HRb1g/edit?usp=sharing)

### Useful resources are provided by:  
1. [Ben J. Callahan and Colleagues: Bioconductor Workflow for Microbiome Data Analysis: from raw reads to community analyses](https://f1000research.com/articles/5-1492/v2).   
2. [Comeau AM and Colleagues: Microbiome Helper: a Custom and Streamlined Workflow for Microbiome Research](http://msystems.asm.org/content/2/1/e00127-16)  
3. [Shetty SA, Lahti L., et al: Tutorial from microbiome data analysis spring school 2018, Wageningen University and Research](https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/introduction.html)  

Note: 
A good practise is to use Rmarkdown for documenting your results and sharing with your collaborators and supervisors. For an introduction to [RStudio](https://www.youtube.com/watch?v=cWJzjHh_3kk&t=337s) and an 
[RStudio Overview](https://www.youtube.com/watch?v=n3uue28FD0w)  


[View this webiste repository on GitHub](https://github.com/microsud/Tools-Microbiome-Anlaysis)  
[Follow me on Twitter](https://twitter.com/gutmicrobe)  
[googlescholar](https://scholar.google.nl/citations?hl=en&user=Vahc6LUAAAAJ&view_op=list_works&sortby=pubdate)  
[ORCID ID: 0000-0001-7280-9915](http://orcid.org/0000-0001-7280-9915)   

### References:
1. Callahan, B. J., McMurdie, P. J. & Holmes, S. P. (2017). Exact sequence variants should replace operational taxonomic units in marker gene data analysis. bioRxiv, 113597.  
2. Callahan, B. J., McMurdie, P. J., Rosen, M. J., Han, A. W., Johnson, A. J. A. & Holmes, S. P. (2016). DADA2: high-resolution sample inference from Illumina amplicon data. Nature methods 13, 581-583.  
3. Caporaso, J. G., Kuczynski, J., Stombaugh, J., Bittinger, K., Bushman, F. D., Costello, E. K., Fierer, N., Peña, A. G., Goodrich, J. K. & Gordon, J. I. (2010). QIIME allows analysis of high-throughput community sequencing data. Nature methods 7, 335-336.  
4. Schloss, P. D., Westcott, S. L., Ryabin, T., Hall, J. R., Hartmann, M., Hollister, E. B., Lesniewski, R. A., Oakley, B. B., Parks, D. H. & Robinson, C. J. (2009). Introducing mothur: open-source, platform-independent, community-supported software for describing and comparing microbial communities. Applied and environmental microbiology 75, 7537-7541.  
5. Team, R. C. (2000). R language definition. Vienna, Austria: R foundation for statistical computing.  

### Was this website/resource useful for you? Then please share it with others too!  
You can cite this resource as:  
Shetty SA and Lahti L (2018). A list of R environment based tools for 16S rRNA gene data exploration, statistical analysis and visualization. [![DOI](https://zenodo.org/badge/99698135.svg)](https://zenodo.org/badge/latestdoi/99698135)  



