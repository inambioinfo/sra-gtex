
##########
### Doing metadata prediction
### using GTEX as a training set
##########

### Load libraries
library(splines)
library(limma)
library(dplyr)
library(readr)

### Load in information on the selected regions
load("/home/bst/faculty/jleek/projects/2016/metapredict/regions_to_subset.rda")
regions = regions_to_subset
rm(regions_to_subset)


#### Figure out which are the unique ones
#### and subset them out 

load("/dcl01/leek/data/gtex_work/runs/sra/DER_analysis/coverageMatrix/ers_gtex/regioninfo.rda")
regioninfo = regioninfo %>% 
  group_by(index,chr) %>% 
  mutate(newtype=paste(unique(type),collapse=",")) %>%
  filter(row_number()==1)

### Load in SRA metadata
load('/dcl01/leek/data/gtex_work/runs/recount2/metadata/metadata_sra.Rdata')
sra_meta = metadata
rm(metadata)

### Load gtex metadata

load("/dcl01/leek/data/gtex_work/runs/sra/DER_analysis/coverageMatrix/ers_gtex/gtexmetadata.rda")
gtex_meta = gtexmetadata
gtex_meta = cbind(gtex_meta,usegtex)
rm(gtexmetadata,usegtex)

### Load SRA coverage data

load('/dcl01/leek/data/gtex_work/runs/sra/DER_analysis/coverageMatrix/ers_gtex/coverageMatrix-cut0.5.Rdata')
sra_cov = coverageMatrix
rm(coverageMatrix)

### Load GTEX coverage data 
load("/dcl01/leek/data/gtex_work/runs/sra/DER_analysis/coverageMatrix/ers_gtex/gtex-coverageMatrix-cut0.5.Rdata")
gtex_cov = coverageMatrixGtex
rm(coverageMatrixGtex)

### Keep only the useful GTEX 

gtex_cov = gtex_cov[,gtex_meta$usegtex]
gtex_meta = gtex_meta[gtex_meta$usegtex,]


### Keep only the good SRA
mm = match(colnames(sra_cov),sra_meta$run)
sra_meta = sra_meta[mm,]
pd = read_csv("https://raw.githubusercontent.com/nellore/runs/master/sra/v2/hg38/SraRunInfo.csv")
sra_meta = left_join(sra_meta,pd,by=c("run"="Run","sample"="Sample"))

####
## Build a predictor of sex in gtex
####

#####
### Build dumb predictor
#####

gtex_sex = gtex_cov[regioninfo$chr %in% c("chrY"),]/(gtex_meta$SumCoverage/4e9)

design = model.matrix(~ Sex,
                      data=gtex_meta)
fit = lmFit(log2(gtex_sex + 1),design)
eb = eBayes(fit)
topSEX = topTable(eb,n=20)

sex_pred = t(gtex_sex[as.numeric(rownames(topSEX)),])
colnames(sex_pred) = paste0("sex_pred",1:20)
sex_pred = as.data.frame(log2(sex_pred+1))
lm1 = lm((gtex_meta$Sex=="male") ~ .,data=sex_pred)


######
### Apply to sra
######


sra_sex = sra_cov[regioninfo$chr %in% c("chrY"),]/(sra_meta$auc/4e9)
sex_pred_sra = log2(t(sra_sex[as.numeric(rownames(topSEX)),])+1)
colnames(sex_pred_sra) = paste0("sex_pred",1:20)

sexhat = predict(lm1,newdata=as.data.frame(sex_pred_sra))


sra_sex_var = tolower(sra_meta$Sex)
usesex = sra_sex_var %in% c("male","female")

table(sexhat[usesex & sra_meta$auc > 339579983] > 0.5, sra_sex_var[usesex & sra_meta$auc > 339579983])
