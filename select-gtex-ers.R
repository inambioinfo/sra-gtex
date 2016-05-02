
##########
### Get GTEX expressed regions associated with 
### sex, age, BMI, and tissue
##########


### Load necessary libraries

library(splines)
library(sva)
library(Hmisc)
library(genefilter)
library(limma)
library(GenomicRanges)


### Load GTEX loading function
source("/dcs01/ajaffe/GTEX/Leek/coverageMatrix/simpleLoad/gtexLoad.R")



### Loop over chromosomes and load data, fitting the models

chrs = c("X","Y",1:22)

index_list = vector(mode="list",length=24)

for(i in 1:24){
  chrname=paste0("chr",chrs[i])
  dat = gtexLoad(chrname)
  pheno = dat$pheno
  cm = dat$coverageMatrix
  cm = log2(cm + 1)
  wid = width(dat$regions)
  
  ### These samples seem to have missing values or to be inappropriate for analysis
  
  use = (pheno$SAMPLE_USE != "Seq_RNA_WTSS" & pheno$SUBJID !="K-562")
  rownames(cm) = as.character(1:dim(cm)[1])
  cm = cm[wid > 20,use]
  pheno = pheno[use,]


  #### Some of these variables have too many levels
  design = model.matrix(~ns(AGE,df=10) + Sex + ns(BMI,10) 
                       + droplevels(Body_Site) + droplevels(SMNABTCHT) + droplevels(SMGEBTCH),data=pheno)
  fit = lmFit(cm,design)
  eb = eBayes(fit)

  #### Take the top 40 of each variable
  topAGE = topTable(eb,2:11,n=40)
  topSEX = topTable(eb,12,n=40)
  topBMI = topTable(eb,13:22,n=40)
  topSite = topTable(eb,23:74,n=40)

  #### Also take the 250 most variable. 
  varSites = rownames(cm)[which(rank(-rowSds(cm)) <= 250)]

  index = c(as.numeric(rownames(topAGE)),as.numeric(rownames(topSEX)),
            as.numeric(rownames(topBMI)),as.numeric(rownames(topSite)),
            as.numeric(varSites))

  type=c(rep(c("age","sex","bmi","tissue"),each=40),
    rep("variable",250))

  out = data.frame(type=type, 
                  index=index)
  save(out,file=paste0(chrname,".rda"))
  index_list[[i]] = out
  cat(i)
}


##########
### Get the unique regions
### from the GTEX selection proces
##########

regions_to_subset = GRanges()
for(i in 1:24){
  chrname=paste0("chr",chrs[i])
  dat = gtexLoad(chrname)
  load(paste0(chrname,".rda"))
  regions_to_subset <- append(regions_to_subset, dat$regions[unique(out$index)])
  rm(out,dat)
  cat(i)
}

save(regions_to_subset,file="regions_to_subset.rda")



#####
### Save GTEX coverage
#####

### Load GTEX coverage

chrs = c("X","Y",1:22)
covmat = matrix(NA,nrow=1,ncol=9662)
setwd("/home/bst/faculty/jleek/projects/2016/metapredict/")
for(i in 1:24){
  chrname=paste0("chr",chrs[i])
  dat = gtexLoad(chrname)
  load(paste0(chrname,".rda"))
  covmat = rbind(covmat,dat$coverageMatrix[unique(out$index),])
  cat(i)
}
covmat = covmat[-1,]

coverageMatrixGtex = covmat

save(coverageMatrixGtex,file="/dcl01/leek/data/gtex_work/runs/sra/DER_analysis/coverageMatrix/ers_gtex/gtex-coverageMatrix-cut0.5.Rdata")


####
## Save GTEX metadata
####

source("/dcs01/ajaffe/GTEX/Leek/coverageMatrix/simpleLoad/gtexLoad.R")
dat = gtexLoad("chrY")
pheno = dat$pheno
gtexmetadata = dat$pheno
usegtex = (pheno$SAMPLE_USE != "Seq_RNA_WTSS" & pheno$SUBJID !="K-562")
save(usegtex,gtexmetadata,file="/dcl01/leek/data/gtex_work/runs/sra/DER_analysis/coverageMatrix/ers_gtex/gtexmetadata.rda")


#### 
## Save region information
####

setwd("/home/bst/faculty/jleek/projects/2016/metapredict/")
load("regions_to_subset.rda")
chrs = c("X","Y",1:22)
regioninfo = data.frame(type=NA,index=NA,chr=NA)
for(i in 1:24){
  chrname=paste0("chr",chrs[i])
  load(paste0(chrname,".rda"))
  out$chr = chrname
  regioninfo = rbind(regioninfo,out)
  cat(i)
}

regioninfo=regioninfo[-1,]
save(regioninfo,file="/dcl01/leek/data/gtex_work/runs/sra/DER_analysis/coverageMatrix/ers_gtex/regioninfo.rda")

