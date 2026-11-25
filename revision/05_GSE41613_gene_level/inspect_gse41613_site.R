suppressPackageStartupMessages(library(data.table))
p <- as.data.table(readRDS(paste0(Sys.getenv("OSCC_DATA"), "/GSE41613_phenotype_local.rds")))
cat('columns\n'); print(names(p))
ch <- grep('characteristics',names(p),value=TRUE)
for(n in ch){cat('\n',n,'\n',sep='');print(unique(p[[n]]))}
