pkgs <- c("MASS")
options(warn = -1)
for (i in pkgs){
  if (!require(i, quietly = TRUE, character.only = TRUE)){
    install.packages(i)
  }
}
devtools::install_github("david-borchers/mt5751a",build = TRUE)