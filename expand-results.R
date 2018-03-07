#!/usr/bin/env Rscript

## libraries
library(randomFunctions, quietly=TRUE) # devtools::install_github("chr1swallace/random-functions")
library(annotSnpStats, quietly=TRUE)# devtools::install_github("chr1swallace/annotSnpStats")
library(snpStats) # bioconductor package
library(magrittr) # CRAN
library(methods) # just in case
library(data.table)
library(flashClust)
library(GUESSFM) # devtools::install_github("chr1swallace/GUESSFM",ref="groups")
library(ggplot2); theme_set(theme_bw())
library(gridExtra)
library(parallel)
NCORES <- 1 # or more?
options(mc.cores=NCORES) 
    library(reshape)
    library(plyr)
    library(speedglm)


## hard coded variables
#DIR="/rds/user/cew54/hpc-work/eqtlfm-mikhail"
## source("~/DIRS.txt")
## DIR=MIKHAILDIR
## EXPRFILE=file.path(DIR,"quantile_normalised_peer_residuals_LCLs.txt")
## ASSOCFILE=file.path(DIR,"all_tested_affinity_expression_associations.txt")
## GENOFILE=function(i) { file.path(DIR,paste0("jGeno-",i,".vcf.gz")) }
## WINDOW=5e+4 # bp either side of index SNP to consider
## MINP.THRESHOLD=1e-5 # don't fine map if no univariate p < MINP.THRESHOLD
## R2WINDOW.THRESHOLD=150 # maximum number of SNPs for r2 calculation to estimate approximate LD blocks. NB - allows blocks to be larger than this, should be ok.
## TAG.R2=0.98
## NSWEEP=50000
## NSAVE=5000
## NCHAINS=5
NEXP=2
CPPTHR=0.99

## PLINK="/scratch/wallace/local/bin/plink" # plink binary
## BCFTOOLS="/scratch/wallace/local/bin/bcftools" # bcftools binary

## ## EQTLFILE="/home/cew54/rds/hpc-work/eqtlfm/Eqtl.22"
## ## OUTDIR="/home/cew54/rds/hpc-work/eqtlfm/working"
## ## data row
## ## todo <- fread(paste0("awk 'NR==1 || NR==",args$row+1,"' ",ASSOCFILE))
## todo <- fread(ASSOCFILE)
## todo <- todo[sig==TRUE,]
## todo2 <- todo[,.(minpos=min(snp_pos),maxpos=max(snp_pos)),by=c("snp_chr","gene_id")]
## todo3 <- todo2[maxpos-minpos>2e+4,]
## todo2 <- todo2[maxpos-minpos<2e+4,] # for now
## todo2[,myid:=gene_id]
## todo[,myid:=paste(gene_id,1:.N,sep="."),by="gene_id"]
## todo[,minpos:=snp_pos]
## todo[,maxpos:=snp_pos]
## todo2 <- rbind(todo2,todo[gene_id %in% todo3$gene_id,colnames(todo2),with=FALSE])


args <- getArgs(list(d="/scratch/cew54/eqtlfm-mikhail/results/ENSG00000204642.2"))
args <- getArgs(list(d="/scratch/cew54/eqtlfm-mikhail/results/ENSG00000164587"))

    message(args$d)
    SKIPFILE <- file.path(args$d,"skip")
    if(file.exists(SKIPFILE))
        next

    files <- list.files(args$d,full=TRUE,pattern="_features")
    if(!length(files)) {
        message("No output files found")
        next
    }

(load(file.path(args$d,"all-data.RData")))

    ## trace plot
    ess <- try(read.ess(files))
    if(class("ess")=="try-error") {
         message("failed to read guess output for ",args$line)
         print(attr(ess,"condition")$message)
         next
    }
png(file.path(args$d,"ess.png"),height=6,width=6,units="in",res=300)
    plot(ess)
dev.off()

    ## load guess 
    dd <- read.snpmod(files)

    ## diffusion plot
    #plot_diffuse(dd)
    bm <- best.models(dd)
    
    pp <- pp.nsnp(dd,expected=2)
    qc.pp <- qc(pp)
    if(qc.pp$flag) {
        system(paste("touch",file.path(args$d,"qcflag")))
        warning("pp.nsnp profile is suspect")
    }

    ## stop if null model has > 50% posterior
    dropmodels <- function(x,minpp=0.5) {
        if(!any(x$str==""))
            return(FALSE)
        w <- which(x$str=="")
        x[w,"PP"]
    }
    nullpp <- dropmodels(bm) 
    message("null model >50% of posterior?")
    print(drop <- nullpp>0.5)
    
    ## these next two lines discard the least interesting models - these
    ## are an approximation to make things run faster - delete if you
    ## want greater accuracy at the cost of lower speed
    dd@models <- dd@models[1:nrow(bm),] # trim
    dd@model.snps <- dd@model.snps[1:nrow(bm)] # trim
    
    ## expand tags
    (load(file.path(args$d,"tags.RData")))
    dx <- expand.tags(dd,tags=tags)
    bm <- best.models(dx,cpp.thr=CPPTHR)
bm <- bm[!is.na(bm$rank),,drop=FALSE]

    ## STOP HERE IF NULL MODEL HAS > 50% POSTERIOR ACROSS ALL TRAITS
    nullpp <- dropmodels(bm) 
    message("null model posterior")
    nullpp
    message("null model >50% of posterior?")
    print(drop <- nullpp>0.5)
    if(drop) {
        system(paste0("touch ",file.path(d,"skip")))
        next
    }

    
    ## refits
    message(ncol(DATA)," SNPs in region")

    nrow(bm)

    ## fits <- structure(vector("list",length(files)), names=names(files))

## SDATA <- as(DATA,"SnpMatrix")
    snp.data0 <- DATA[ , tags@.Data ]
    ## rn <- rownames(snp.data0)
    ## rs <- row.summary(snp.data0)
    ## skip qc
    ## use <- rs[,"Call.rate"]==1 & complete.cases(covars) & complete.cases(cc.ph)
todo <- unique(c(colnames(DATA),bm$str))
fits <- abf.calc(y=y,x=DATA,models=todo,family="gaussian",
                     snp.data=snp.data0)

    ## make snpmods
    SM <- abf2snpmod(fits,expected=NEXP,overdispersion=1,nsnps=length(tags@.Data))

    ## SM contains all information you need - save it somewhere
    ## you can look at the posterior distribution of the number of causal snps
    pp.nsnp(SM)
    plot(pp.nsnp(SM))

    ## most likely models
    best.models(SM,pp.thr=0.01)

    ## most likely snps
    best.snps(SM,pp.thr=0.01)

    ## marginal posterior probability of inclusion for each SNPlocs
    mppi <- best.snps(SM,pp.thr=0)

    of <-file.path(args$d,paste0(basename(args$d),"-snpmod-",
                                      sub("0.","",as.character(CPPTHR)),
                                      ".RData"))
    message("saving SM to ",of)
save(SM,file=of)

of  %<>%  sub("snpmod.*","mppi.csv",.)
z <- best.snps(SM,pp.thr=0)
write.table(z[,c("var","Marg_Prob_Incl")],file=of,row.names = FALSE,sep="\t",quote=FALSE)
of  %<>%  sub("mppi.*","nsnp.csv",.)
z <- as.data.frame(pp.nsnp(SM)$trait)
z$NCausalVars <- as.numeric(rownames(z))
colnames(z)[1] <- "PostProb"
write.table(z[,2:1],file=of,row.names = FALSE,sep="\t",quote=FALSE)
