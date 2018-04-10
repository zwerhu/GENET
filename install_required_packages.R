cran.packages <- c("glmnet", "methods", "dplyr", "reshape2")

gv_install_packages <- function(cran.packages, bioc.packages, verbose = TRUE) {
  # Ensures that all the CRAN and biocLite packages provided are installed and loaded
  # Args:
  #	cran.packages: vector of CRAN packages
  #	bioc.packages: vector of Bioconductor packages
  if (length(cran.packages)>=1) {
    message("Checking CRAN dependencies")
    new.packages <- cran.packages[!(cran.packages %in% installed.packages()[,"Package"])] # identify missing packages
    if(length(new.packages)) {
      install.packages(new.packages, repos = "http://cran.us.r-project.org", dependencies=TRUE) # install missing packages
      message("Necessary packages installed")
    }
    else message ("All necessary packages already installed")
    invisible(lapply(cran.packages, require, character.only=TRUE)) # to load the list of required packages
    message("All necessary CRAN packages loaded")
  }
  
  if (length(bioc.packages)>=1) {
    message("Checking Bioconductor dependencies")
    source("https://bioconductor.org/biocLite.R")
    biocLite(ask=FALSE) # updates all packages
    biocLite(bioc.packages, ask=FALSE)
    message("Necessary packages installed")
    invisible(lapply(bioc.packages, require, character.only=TRUE)) # to load the list of required bioconductor packages
    message("All necessary Bioconductor packages loaded")
  }
}


gv_install_packages(cran.packages)
