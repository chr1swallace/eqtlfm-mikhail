#!/usr/bin/env Rscript
library(data.table)
library(mvtnorm)
library(randomFunctions)
library(magrittr)

args <- getArgs(defaults=list(g1="ENSG00000204531.2",
                              g2="ENSG00000234745"))
args <- getArgs(defaults=list(g1="ENSG00000165985",
                              g2="ENSG00000165983"))
source("~/DIRS.txt")
print(args)

d1 <- file.path(MIKHAILDIR,"results",args$g1)

library(snpStats)
(load(file.path(d1,"all-data.RData")))
y1 <- y
N <- as(DATA,"numeric")

## add in y2
DIR=MIKHAILDIR
EXPRFILE=file.path(DIR,"standard_normal_transformed_peer_residuals_LCLs.txt") #quantile_normalised_peer_residuals_LCLs.txt")
 expr <- fread(paste0("awk 'NR==1 || $464==\"",args$g2,"\"' ",EXPRFILE),sep="\t") #[,-c(1,464)] ## nb warning about 1 too few column names is header is expected
if(!(nrow(expr)==1)) { message("gene ",args$g2," not uniquely found"); stop()} ## check
y2 <- t(expr[,names(y1),with=FALSE])[,1]
y2 <- y2 - mean(y2) # center, as for y1

## N <- N[,1:100]
## DATA <- DATA[,1:100]
## dim(N)

## step one, conditional testing
library(GUESSFM)
cond <- function(X,Y) {
    best <- NULL
    ret <- cond.best(X,Y)
    while (length(newbest <- cond.best(X, Y, best, stepwise.p.thr=1e-6, family="gaussian"))) {
        best <- c(best, newbest)
    }
    best
}

y <- list(y1,y2)
names(y) <- c("g1","g2")
condres <- lapply(y, function(Y) { cond(DATA,Y) })

cs <- col.summary(DATA)
## LD <- snpStats::ld(DATA,depth=200,stat="R",symmetric=TRUE)
devtools::load_all("~/RP/coloc")
LD <- cor2(N)

getmods <- function(y,f) {
    Ndf <- as.data.frame(N)
    Ndf$y <- y
    mods <- data.table(snp=colnames(N))
    Ndf$yr <- lm(paste("y ~ ",f), data=Ndf)  %>%  residuals()
    for(i in 1:ncol(N)) {
        fi <- paste("yr ~",colnames(N)[i])
        m <- lm(fi, data=Ndf)
        mods$B[i] <- coefficients(m)[2]
        mods$V[i] <- vcov(m)[2,2]
    }
    mods
}
    
message("models to fit:")
cmods <- lapply(1:2, function(i) {
    csnps <- condres[[i]]
    if(length(csnps)==1)
        return("1")
    sapply(1:length(csnps), function(i) {
        ## paste(c("1",csnps)[1:i],collapse="+")
        paste(csnps[-i],collapse="+")
    })
})
print(cmods)

## summary stats for each model
mods <- lapply(1:2, function(i) {
    lapply(cmods[[i]], function(f) getmods(y[[i]], f))
})

do.coloc <- function(mod1,mod2) {
    coloc.abf(dataset1=list(snp=mod1$snp,
                            beta=structure(mod1$B,names=mod1$snp),
                            varbeta=structure(mod1$V,names=mod1$snp),
                            N=length(y1),sdY=sd(y1),type="quant"),
              dataset2=list(snp=mod2$snp,
                            beta=structure(mod2$B,names=mod2$snp),
                            varbeta=structure(mod2$V,names=mod2$snp),
                            N=length(y2),sdY=sd(y2),type="quant"),
              MAF=structure(cs[,"RAF"],names=colnames(DATA)))
}
library(lazyeval)
myouter <- function(L1,L2, FUN) {
    todo <- expand.grid(i=seq_along(L1),j=seq_along(L2))
    ret <- lapply(1:nrow(todo), function(r) {
        FUN(L1[[todo$i[r]]], L2[[ todo$j[r] ]])$summary
    })  %>% do.call("rbind",.)  %>% as.data.table()
    ret[,cond1:=cmods[[1]][todo$i]]
    ret[,cond2:=cmods[[2]][todo$j]]
    ## ret[,method:=lazyeval::expr_text(FUN)]
    ret
}


ret <- myouter(mods[[1]], mods[[2]], do.coloc)
ret$g1 <- args$g1
ret$g2 <- args$g2
ret

write.table(ret,row.names=FALSE,quote=FALSE,file=file.path(MIKHAILDIR,"coloc-v5",paste0(args$g1,"-",args$g2,".csv")))
print(ret)
