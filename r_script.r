if(!(require(rmarkdown))){ 
install.packages("rmarkdown")}
library(rmarkdown)
args <- commandArgs(TRUE)
f <- as.character(args[1])
d <- as.character(args[2])
rmarkdown::render(f,output_dir=d)
