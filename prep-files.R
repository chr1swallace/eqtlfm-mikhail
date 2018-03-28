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
library(randomFunctions) # devtools::install_github("chr1swallace/random-functions")

## hard coded variables
#DIR="/rds/user/cew54/hpc-work/eqtlfm-mikhail"
source("~/DIRS.txt")
DIR=MIKHAILDIR
EXPRFILE=file.path(DIR,"standard_normal_transformed_peer_residuals_LCLs.txt") #quantile_normalised_peer_residuals_LCLs.txt")
ASSOCFILE=file.path(DIR,"all_tested_affinity_expression_associations_v2.txt")
GENOFILE=function(i) { file.path(DIR,paste0("jGeno-",i,".vcf.gz")) }
WINDOW=2e+5 # bp either side of index SNP to consider
MINP.THRESHOLD=1e-5 # don't fine map if no univariate p < MINP.THRESHOLD
R2WINDOW.THRESHOLD=150 # maximum number of SNPs for r2 calculation to estimate approximate LD blocks. NB - allows blocks to be larger than this, should be ok.
TAG.R2=0.98
NSWEEP=500000
NSAVE=5000
NCHAINS=5
NEXP=3

PLINK="/scratch/wallace/local/bin/plink" # plink binary
BCFTOOLS="/home/cew54/locald/bin/bcftools" # bcftools binary

## EQTLFILE="/home/cew54/rds/hpc-work/eqtlfm/Eqtl.22"
## OUTDIR="/home/cew54/rds/hpc-work/eqtlfm/working"
## data row
## todo <- fread(paste0("awk 'NR==1 || NR==",args$row+1,"' ",ASSOCFILE))
todo <- fread(ASSOCFILE)
todo <- todo[sig==TRUE,]
todo2 <- todo[,.(minpos=min(snp_pos),maxpos=max(snp_pos)),by=c("snp_chr","gene_id")]
todo3 <- todo2[maxpos-minpos>1e+5,]
todo2 <- todo2[maxpos-minpos<=1e+5,] # for now
todo2[,myid:=gene_id]
todo[,myid:=paste(gene_id,1:.N,sep="."),by="gene_id"]
todo[,minpos:=snp_pos]
todo[,maxpos:=snp_pos]
todo2 <- rbind(todo2,todo[gene_id %in% todo3$gene_id,colnames(todo2),with=FALSE])

fwrite(todo2,file="assoc-processed.tab")
# i=42
for(i in sample(1:nrow(todo2))) {
    message("\n!!!",todo2$myid[i],"\n")
    outd <- file.path(DIR,"results",todo2$myid[i])
    SKIPFILE <- file.path(outd,"skip")
    COMFILE <- file.path(outd,"runme.sh")
    if(file.exists(COMFILE) || file.exists(SKIPFILE))
        next
    if(!file.exists(GENOFILE(todo2$snp_chr[i])))
        next
    
    ## read in expression data
    expr <- fread(paste0("awk 'NR==1 || $464==\"",todo2$gene_id[i],"\"' ",EXPRFILE),sep="\t")[,-c(1,464)]
    if(!(nrow(expr)==1)) { stop("gene ",todo2$gene_id[i]," not uniquely found") } ## check
    
    ## read in genetic data
    tmp <- tempfile()
    
    ## /home/cew54/local/bin/bcftools view  -f \"VTYPE==SNP\" /scratch/cew54/eqtlfm-mikhail/jGeno-16.vcf.gz --regions 16:30994897-31094897 |grep rs112607901
comm <- paste0(BCFTOOLS," view  -v snps ",GENOFILE(todo2$snp_chr[i]),
                  " --regions ",todo2$snp_chr[i],":",todo2$minpos[i]-WINDOW,"-",todo2$maxpos[i]+WINDOW," -Ob | ",
                   PLINK," --bcf /dev/stdin --keep-allele-order --vcf-idspace-to _ --const-fid --allow-extra-chr 0 --split-x b37 no-fail --make-bed --out ", tmp)
    system(comm)
    ## x <- read.plink(tmp)
    X <- annot.read.plink(tmp)
    int <- intersect(colnames(expr),rownames(X))
    y <- t(expr[,int,with=FALSE])[,1]
    y <- y-mean(y) # center
    X <- X[int,]
    X@samples$y <- y
    phenotype(X) <- "y"

    message("\tsnps found: ",ncol(X))
    ## what is min p?  Only go ahead if min p < 10-5
    ss <- snp.rhs.tests(X@samples$y ~ 1, snp.data=sm(X), family="Gaussian")
    p <- p.value(ss)
    names(p) <- colnames(X)
    if(!file.exists(outd)) {
        dir.create(outd,recursive = TRUE)
    }
    if(min(p) > MINP.THRESHOLD) {
        system(paste("touch",SKIPFILE))
        next()
    }

    sinfo <- samples(X)
    
    ## runguess
    message("Running GUESS")
    options(scipen=1000000000)
    
    ## tags
    f.tags <- file.path(outd,"tags.RData")
    if(TAG.R2<1 && !file.exists(f.tags)) {
        message("Preparing to tag.")
        DATA <- as(X,"SnpMatrix")
        message("Input matrix has ",ncol(DATA)," SNPs.")

        ## commented, but here is where QC normally applied
        cs0 <- col.summary(DATA)
        wh <- which(is.na(cs0[,"z.HWE"]) |
                    cs0[,"MAF"]<0.001 |
                    cs0[,"Call.rate"]<0.99 |
                    cs0[,"Certain.calls"]<0.75 | 
                    abs(cs0[,"z.HWE"])>4)
        if(length(wh)) {
            message("Dropping ",length(wh)," SNPs with |z.HWE|>5, MAF < 0.001 in controls or call rate <0.99")
            DATA <- DATA[,-wh]
        }
        tags <- tag(DATA, method="single", tag.threshold=TAG.R2)
        message("after tagging at ",TAG.R2,", matrix now has ",length(unique(tags(tags)))," SNPs.")
        save(tags, file=f.tags)
    }

    save(DATA,y,file=file.path(outd,"all-data.RData"))

    f.data <- file.path(outd,"data.RData")
    f.par <- file.path(outd, "par.xml")
    if(file.exists(f.par) ) {
        message("output file already exists, skipping: ",f.par)
    }
    DATA <- DATA[, unique(tags(tags))]

  ## PCs
  ## covars <- sinfo[wh.DATA,c("PC1","PC2","PC3","PC4")]
  ## fix <- which(sapply(covars,class)=="matrix")
  ## if(length(fix)) {
  ##   for(j in fix) 
  ##     covars[,j] <- as.vector(covars[,j])
  ## }
  
    ## save data
    y <- samples(X)$y
    save(DATA,y,file=f.data)
 
  message("Samples: ", length(y))
  message("SNPs: ",ncol(DATA))
  #message("Tags: ",length(unique(tags(tags))))

  ## command=paste("/home/chrisw/local/bin/GUESS")
  coms <- run.bvs(X=DATA,Y=y, #covars=covars,
                  gdir=outd,nsweep=NSWEEP, family="gaussian",tag.r2=NA,
                  nsave=NSAVE,nchains=NCHAINS,nexp=NEXP,run=FALSE)

    cat(coms,file=COMFILE,sep="\n")
}
