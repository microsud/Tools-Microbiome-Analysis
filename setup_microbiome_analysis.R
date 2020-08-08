
############################################################################################################################################
#																																	                                                                    	   #
#						                            START OF CODE								                                                										   #
#																														                                                                    				   #
############################################################################################################################################

# These are some R pacakges I commonly use for analysis. Everytime I update R to  latest version, I run this installation to in setup my environment.  

# Copy from here

setup_microbiome_analysis <- function(){
  
  .packages = c("ape", 
                "gridExtra", 
                "picante", 
                "data.table", 
                "RColorBrewer", 
                "DT", 
                "reshape", 
                "reshape2", 
                "magrittr", 
                "markdown",
                "ggpubr", 
                "tibble", 
                "pheatmap", 
                "dplyr", 
                "viridis", 
                "devtools", 
                "rmdformats",
                "intergraph",
                "network",
                "igraph",
                "ggplot2", 
                "gridExtra", 
                "knitr", 
                "vegan", 
                "plyr", 
                "dplyr",
                "ggrepel", 
                "ggnetwork", 
                "ade4", 
                "rmarkdown",
                "formatR",
                "caTools",
			"GGally")
  
  .bioc_packages <- c("phyloseq",
                      "microbiome", 
                      "phangorn", 
                      "genefilter")
  
  # Install CRAN packages (if not already installed)
  .inst <- .packages %in% installed.packages()
  if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])

  .inst <- .bioc_packages %in% installed.packages()
  if(any(!.inst)) {
    source("http://bioconductor.org/biocLite.R")
    biocLite(.bioc_packages[!.inst], ask = F)
  }
  
  
  if (!("microbiomeutilities"  %in% installed.packages())) {
    devtools::install_github("microsud/microbiomeutilities")
  }
  
  if (!("SpiecEasi"  %in% installed.packages())) {
    devtools::install_github("zdk123/SpiecEasi")
  }
  
  if (!("ggnet"  %in% installed.packages())) {
    devtools::install_github("briatte/ggnet")
  }
  
  message("If there was no error then you are ready to do microbiome data analysis")
  
}

setup_microbiome_analysis()

# Copy until here previous line!

####################################################END OF CODE###########################################################################
