# Load packages quietly
sshhh <- function(a.package){
  suppressWarnings(suppressPackageStartupMessages(
    library(a.package, character.only=TRUE)))
}


# Load R essentials
load_R_essentials <- function() {
  sshhh("readr")
  sshhh("tidyr")
  sshhh("dplyr")
  sshhh("stringr")
  sshhh("magrittr")
  sshhh("cowplot")
  sshhh("scales")
}