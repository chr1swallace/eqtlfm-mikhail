## setwd('/home/ng414/nfg/coloccc')

library(mvtnorm)
library(MASS)
library(simGWAS)
library(parallel)
library(coloc)
library(randomFunctions)
library(snpStats)
                                        # devtools::load_all("/home/cew54/RP/coloc")
library(data.table)

#univariate analysis for genomatrix g and trait y
maker2 <- function(g, y) {
    ss <- snp.rhs.estimates(y ~ 1, family = "gaussian", snp.data = g)
    b <- sapply(ss, "[[", "beta")
    names(b) <- colnames(g)
    list(snp = colnames(g),
         beta = b,
         varbeta = sapply(ss,"[[","Var.beta"),
         N = nrow(g),
         sdY = sd(y),
         type = "quant"
      )
}

## source("functions.R")
## source("/rds/user/ng414/hpc-work/coloccc/qq_cew54.R")

#simulation parameters
nmax <- 1000
N <- 201

## outputdir <- sprintf("~/scratch/qqsim/qq_sym_N%s-nmax%s-%s", N, nmax, sig)
outputdir <- sprintf("/home/cew54/scratch/qqsim")
if(!dir.exists(outputdir)) dir.create(outputdir,recursive = TRUE)


  args <- getArgs(defaults=list(N = 201), numeric = "N")
  ## TODO - set temp file prefix - then won't need this next line
  file.ldd="/home/cew54/share/Data/reference/lddetect/EUR/fourier_ls-chr22.bed"
  file.vcf="/home/cew54/share/Data/reference/UK10K/BCF/chr22.bcf.gz"

  ## ldblocks
  ldd <- fread(file.ldd)

  ## split bcf by ldblocks
  ldd[,blocknum:=1:.N]
  ldd[,dist:=stop-start]
  ldd[,comm:=paste0("/home/cew54/localc/bin/bcftools view ",file.vcf,
                    " --min-af 0.01:minor --max-alleles 2 --min-alleles 2 ",
                    " -r chr",22,":",start,"-",stop," -Ov ")] # -o ",tmp)]
  gethap <- function(i) {
      y=fread(ldd$comm[i])
      ha <- simGWAS:::vcf2haps(as.matrix(y[,-c(1:9)]))
      rownames(ha) <- paste0("pos",y$POS)
      t(ha)
  }

    # rho <- sample(c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), 1)
    block <- sample(nrow(ldd), 1)
    h <- gethap(block) # rows=samples, cols=snps
    use <- apply(h,2,var)>0 & colMeans(h) > 0.01 & colMeans(h)<0.99  ## no monomorphs
    h <- h[ ,use,drop=FALSE]

    if(ncol(h)>nmax)
    h <- h[,1:nmax] ## max 1000 snps, for speed in simulations
    maf <- colMeans(h)
    LD <- coloc:::cor2(h)

    ## LOR <- rnorm(200,sd=0.2)
    ## LOR <- sample( LOR[ abs(LOR) > 0.05 ], 2) # just exclude the values vv close to 0
    message("nsnps: ",ncol(h))
    ## message("scenario: effects: ", round(LOR[1], 3), " ", round(LOR[2], 3))

    ## generate genotypic data
    G <- new("SnpMatrix",
                h[sample(1:nrow(h),args$N,replace=TRUE),] +
                h[sample(1:nrow(h),args$N,replace=TRUE),] + 1)
    cs <- col.summary(G)
    G <- G[, cs$MAF != 0 & !is.na(cs$z.HWE)]

    #  #generating correlated errors, fixed rho, variable variance
    #  sig <- sample(c(0.1, 0.2, 0.3, 0.5, 0.8, 1), 1) #noise variance
    #  rho <- 0.5 #noise correlation
    ##generating correlated errors, varied rho, fixed variance

RESULT <- lapply(rep(c("h0", "h1", "h3", "h4"),100), function(hyp) {
    CVS <- sample(which(maf>0.1 & maf<0.9), 2) # pick common causal SNPs for power
    #simulate effect size
    LOR <- sample(c(0.2,0.4,0.6),2,replace=TRUE)
sig <- sample(c(0.1,0.2,0.4),1) #noise variance
    rho <- sample(c(0, 0.2,  0.5,  0.8), 1)
## rho <- 0
    sig.mat <- matrix(c(sig, rho * sig, rho * sig, sig), nrow = 2)
    error <- mvrnorm(n = nrow(G), mu = c(0, 0), Sigma = sig.mat)
  #generate first output
    if(hyp=="h0") {
          CVS[1] <- 0
          y1 <- error[, 1]
      } else {
          y1 <- LOR[1] * as.numeric(G[, CVS[1]]) + error[, 1]
    }

    #generate second output
    if(hyp=="h4") { # as above; h4=same source of variation
         CVS[2] <- CVS[1]
         y2 <- LOR[2] * as.numeric(G[, CVS[1]]) + error[, 2]
    }
    if(hyp=="h3") { # new cv; h3=two different sources of variation
         y2 <- LOR[2] * as.numeric(G[, CVS[2]]) + error[, 2]
    }
    if(hyp %in% c("h0","h1")) { # no signal
         CVS[2] <- 0
         names(CVS)[2] <- "--"
         y2 <- error[, 2]
    }

    ## create input data
    s1 <- maker2(G, y1)
    s2 <- maker2(G, y2)

    N <- args$N; p1 <- 1e-4; p2 <- 1e-4; p12 <- 1e-5

    message("running coloc")
    result <- as.data.table(as.list(coloc.abf(s1, s2, MAF = maf)$summary)) # classic analysis

    result[,hyp:=hyp]
    ## result[,test:=c("coloc","coloc-qq")]
    result[,N:=args$N]
    result[,block:=block]
    result[,cv:=colnames(h)[CVS[1]]]
    result[,dv:=colnames(h)[CVS[2]]]
    result[,beta1:=LOR[1]]
    result[,beta2:=LOR[2]]
    result[,minp1:=min(pchisq(s1$beta^2/s1$varbeta,1,lower=FALSE))]
    result[,minp2:=min(pchisq(s2$beta^2/s2$varbeta,1,lower=FALSE))]
    result[,r:=if(hyp=="h3") { LD[ CVS[1], CVS[2] ] } else { NA }]
    result[, sigsq:=sig]
    result[, rho:=rho]
    result
})

result <- rbindlist(RESULT)
of <- tempfile(tmpdir=outputdir,fileext=".rds")

message("saving to ", of)
saveRDS(result, file = of)


if(!interactive())
    q("no")

