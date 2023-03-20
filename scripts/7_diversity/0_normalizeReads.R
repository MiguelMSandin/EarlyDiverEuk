#----
#---- loading packages ----

library(data.table)
library(tidyr)
library(dplyr)

#----
setwd("~/Desktop/Uppsala/1_ecoEvo/repository/resources/")
#----

data = fread("otu_reads.tsv")
colnames(data) = c("otus", "reads")

data$db = ifelse(grepl("^PacBio_", data$otus), "PacBio", "PR2")

(max = data %>% group_by(db) %>% summarise(max=max(reads)))
# The OTU with the maximum numer of reads of the PR2 database (which is the database with the smaller maximum number of reads) is 1418
# So we will normalize the reads at 1000
data$norm = ifelse(grepl("^PacBio_", data$otus), data$reads/max$max[1]*1000, data$reads/max$max[2]*1000)
data$norm = ceiling(data$norm)

dataOut = data %>% select(otus, norm)
write.table(dataOut, "tipNames_OTUreads_norm.tsv", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

#----
