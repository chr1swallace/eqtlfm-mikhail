library(data.table)
library(magrittr)
library(ggplot2)
library(cowplot)
library(ggpubr)

#assemble and save
ff <- list.files("~/scratch/qqsim",pattern="rds",full=TRUE)
length(ff)
res <- lapply(ff, readRDS)  %>% rbindlist()

res <- res[order(hyp,sigsq),]
res[,.(median(minp1),median(minp2)),by=c("hyp","sigsq")]
res <- res[sigsq!=0.3,]
res <- res[rho %in% c(0,0.2,0.5,0.8),]

table(res$rho,res$sigsq)

#plot all
pby <- function(x) {
   print(table(x$hyp))
   xs <- x[,lapply(.SD,mean),by=c("hyp","rho","sigsq"), .SDcols=2:6]
   xs <- melt(xs, c("hyp","rho","sigsq"))
   xs[,variable:=gsub("PP.|.abf","",variable)]
   h <- toupper(unique(xs$hyp))
   xs[,variable:=factor(variable,
                            levels=paste0("H",0:4))]
   xs[,simhyp:=variable==toupper(hyp)]
   print(xs[simhyp==TRUE,])
   xs[,sigsq:=paste0("s2 = ",sigsq)]
   xs[,hyp:=paste0("sim = ",toupper(hyp))]
   xs[,rho:=factor(rho)]
   ggplot(xs,aes(x=rho,y=value,fill=variable)) +
     geom_col(position="stack") +
     facet_grid(hyp ~ sigsq) +
     ## scale_x_continuous("rho",breaks=seq(0,0.9,by=0.1)) +
     scale_fill_manual("posterior",
                       values=c("H0"="white",
                                "H1"="grey30",
                                "H2"="grey70",
                                "H3"="blue",
                                "H4"="orange"))  +
     theme(axis.text.x = element_text(angle = 90)) +
     ylab("Average posterior probability")
}
pby(res) + theme_pubr() + theme(legend.position="right")
## pby(res) + labs_pubr() + theme(legend.position="right")
ggsave("~/qqsim.pdf",height=8,width=8)
##

#plot qq
pby <- function(x) {
   print(table(x$hyp))
   x[,rhocat:=cut(rho,c(0,0.05,0.33,0.66,1),include.lowest=TRUE)]
   m <- melt(x,c("hyp","sigsq","rhocat"),grep("abf",names(x)))
   m[,variable:=gsub("PP.|.abf","",variable)]
   m[,variable:=factor(variable,
                            levels=paste0("H",0:4))]
   m[,x:=factor(paste(sigsq,rhocat))]
   xs <- copy(m)
   xs <- m[,.(iq=seq(0.01,0.99,by=0.01),
              q=quantile(value,seq(0.01,0.99,by=0.01))),
           by=c("hyp","sigsq","rhocat","variable")]
   xs[,uncor:=rhocat=="[0,0.05]"]
   xs <- merge(xs[uncor==FALSE,],xs[uncor==TRUE,],
               by=c("hyp","sigsq","variable","iq"),
               suffixes=c(".un",".cor"))
   ## xs <- m[,.(m=mean(value),
   ##            l=quantile(value,0.25),
   ##            u=quantile(value,0.75)), by=c("hyp","sigsq","rhocat","variable")]
   ## xs <- m[,.(m=median(value),
   ##            l=quantile(value,0.25),
   ##            u=quantile(value,0.75)), by=c("hyp","sigsq","rhocat","variable")]
   ## ggplot(xs,aes(x=rhocat+sigsq,y=m,ymin=l,ymax=u,col=variable)) +
   ## ggplot(xs,aes(x=x,y=m,ymin=l,ymax=u,col=variable)) +
   ggplot(xs,aes(x=q.un,y=q.cor,col=rhocat.un)) +
     ## geom_col(position="stack") +
     geom_abline() + 
     geom_path() +
     facet_grid(hyp ~ variable + sigsq) +
     ## scale_x_continuous("rho",breaks=seq(0.1,0.9,by=0.1)) +
     scale_fill_manual(values=c("H0"="white",
                                "H1"="grey30",
                                "H2"="grey70",
                                "H3"="blue",
                                "H4"="orange"))  +
     theme(axis.text.x = element_text(angle = 90)) +
     ylab("Average posterior probability")
}
pby(res)
##

#ridges 
library(ggridges)
pby <- function(x) {
   print(table(x$hyp))
   x[,rhocat:=cut(rho,c(0,0.05,0.5,1),include.lowest=TRUE)]
   m <- melt(x,c("hyp","sigsq","rhocat"),grep("abf",names(x)))
   m[,variable:=gsub("PP.|.abf","",variable)]
   m[,variable:=factor(variable,
                            levels=paste0("H",0:4))]
   m[,x:=factor(paste(sigsq,rhocat))]
   xs <- copy(m)
   ## xs <- m[,.(m=mean(value),
   ##            l=quantile(value,0.25),
   ##            u=quantile(value,0.75)), by=c("hyp","sigsq","rhocat","variable")]
   ## xs <- m[,.(m=median(value),
   ##            l=quantile(value,0.25),
   ##            u=quantile(value,0.75)), by=c("hyp","sigsq","rhocat","variable")]
   ## ggplot(xs,aes(x=rhocat+sigsq,y=m,ymin=l,ymax=u,col=variable)) +
   ## ggplot(xs,aes(x=x,y=m,ymin=l,ymax=u,col=variable)) +
   ggplot(xs,aes(x=value,y=x,height=..density..)) +
     ## geom_col(position="stack") +
     geom_density_ridges(stat="density") +
     facet_grid(hyp ~ variable) +
     ## scale_x_continuous("rho",breaks=seq(0.1,0.9,by=0.1)) +
     theme(axis.text.x = element_text(angle = 90)) +
     ylab("Average posterior probability")
}
pby(res) 
##


#boxplots 
pby <- function(x) {
   print(table(x$hyp))
   x[,rhocat:=cut(rho,c(0,0.05,0.5,1),include.lowest=TRUE)]
   m <- melt(x,c("hyp","sigsq","rhocat"),grep("abf",names(x)))
   m[,variable:=gsub("PP.|.abf","",variable)]
   m[,variable:=factor(variable,
                            levels=paste0("H",0:4))]
   m[,x:=factor(paste(sigsq,rhocat))]
   xs <- copy(m)
   ## xs <- m[,.(m=mean(value),
   ##            l=quantile(value,0.25),
   ##            u=quantile(value,0.75)), by=c("hyp","sigsq","rhocat","variable")]
   ## xs <- m[,.(m=median(value),
   ##            l=quantile(value,0.25),
   ##            u=quantile(value,0.75)), by=c("hyp","sigsq","rhocat","variable")]
   ## ggplot(xs,aes(x=rhocat+sigsq,y=m,ymin=l,ymax=u,col=variable)) +
   ## ggplot(xs,aes(x=x,y=m,ymin=l,ymax=u,col=variable)) +
   ggplot(xs,aes(x=x,y=value,col=variable,fill=sigsq)) +
     ## geom_col(position="stack") +
     geom_boxplot(notch=TRUE) +
     facet_grid(hyp ~ variable) +
     ## scale_x_continuous("rho",breaks=seq(0.1,0.9,by=0.1)) +
     ## scale_fill_manual(values=c("H0"="white",
     ##                            "H1"="grey30",
     ##                            "H2"="grey70",
     ##                            "H3"="blue",
     ##                            "H4"="orange"))  +
     theme(axis.text.x = element_text(angle = 90)) +
     ylab("Average posterior probability")
}
pby(res)
##


#plot classic coloc
x <- zz[test == 'coloc', ]
x$rho <- as.factor(x$rho)
xs <- x[,lapply(.SD,mean),by=c("hyp","rho"), .SDcols=2:6]
xs <- melt(xs, c("hyp","rho"))
xs[,variable:=gsub("PP.|.abf","",variable)]
h <- toupper(unique(xs$hyp))
xs[,variable:=factor(variable,
                        levels=paste0("H",0:4))]
xs[,simhyp:=variable==toupper(hyp)]
print(xs[simhyp==TRUE,])
ggplot(xs,aes(x=rho,y=value,fill=variable)) + geom_col(position="stack") +
 facet_grid(paste("hyp ~ .")) +
 scale_fill_manual(values=c("H0"="white",
                            "H1"="grey30",
                            "H2"="grey70",
                            "H3"="blue",
                            "H4"="orange"))  +
 theme(axis.text.x = element_text(angle = 90)) +
 xlab("") + ylab("Average posterior probability")

x[, max.lor:= ifelse(abs(lor1) > abs(lor2), abs(lor1), abs(lor2))]
ggplot(dat = x, aes(x = max.lor, y = PP.H4.abf)) + geom_point() + facet_grid("hyp ~ rho")


##
