## Do this first "qrsh -l leek -q leek.q@compute-085"

library(dplyr)
library(stringr)
library(GenomicFeatures)


jx_sql = src_sqlite(path="/dcl01/leek/data/cwilks/combined_junctions")
jx = tbl(jx_sql,"intron")


### Try a really simple thing

chrominfo <- getChromInfoFromUCSC(genome="hg38")

## Try counting commas for the different files

njx = jx %>% group_by(chromosome) %>% summarize(count = n()) %>% collect()

chrs = paste0("chr",c(1:25,"X","Y"))

for(i in 1:length(chrs)){

}



