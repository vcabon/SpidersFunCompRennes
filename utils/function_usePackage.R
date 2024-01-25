# Function install and/or load package

usePackage <- function(p) {
    if(!is.element(p, installed.packages()[,1])) {
        install.packages(p, dep = TRUE)
        library(p, character.only = TRUE)
    }
    else {
        library(p, character.only = TRUE)
    }
}
