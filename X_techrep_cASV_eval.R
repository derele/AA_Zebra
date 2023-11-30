library(tidyverse)
library(MultiAmplicon)
library(Hmisc)
library(Matrix)
library(igraph)

redoBioInf <- TRUE

## Using my fork version of phyloseq for this (for now)
devtools::load_all("/home/ele/git_projects/phyloseq/")


if(redoBioInf){
    source("1_Zebra_general_MA.R")
} else {
    PS <- readRDS("intermediate_data/phyloseqRAW.RDS")
}


## This gives us MA7, MAsample, MA8 and PS, which is derived from MA8
## ... need to improve on this!


dupes <- sample_data(PS)$Animal.No[duplicated(sample_data(PS)$Animal.No)]
table(unique(sample_data(PS)$Animal.No)%in%dupes)
##  FALSE  TRUE 
##   106   131

## only 131 were duplicated and can be used here...
## and we also exclude the negative conrols

P <- subset_taxa(PS, taxa_sums(PS)>100)

P <- subset_samples(P,
                    Animal.No%in%dupes &
                    !grepl("neg|Neg|NEG", Animal.No))


### Susanna's contribution!!!
combineASVs <- function(phyloseq, taxon.level="genus", p.valt = 0.05) {
    taxa <- get_taxa_unique(phyloseq, taxon.level)
    for (i in 1:length(taxa)){
        kaza <- prune_taxa(tax_table(phyloseq)[,taxon.level]%in%taxa[i], phyloseq)
        kaza <- kaza@otu_table
        ## correlation matrix################
        otu.cor <- rcorr(as.matrix(kaza), type="spearman")
        otu.pval <- forceSymmetric(otu.cor$P)
        cor.p <- p.adjust(otu.pval, method="BH") # adjusting for multiple testing
        otu.pval@x<- cor.p
        p.yes <- otu.pval < p.valt # only significant p values
        r.val <- otu.cor$r # select all the correlation values
        p.yes.r <- r.val*p.yes # only select correlation values based on p-value criterion
        ## sanity check
        if(!all(rownames(p.yes.r)==colnames(kaza))){
            stop("Error rownames in correlation matrix do not match ASV names, please report this as a bug at https://github.com/derele/MultiAmplicon/issues")
        }
        ## network basded on the correlation adjancency matrix
        adjm <- as.matrix(p.yes.r)
        ##ignoring NAs
        adjm[is.na(adjm)] <- 0
        net.grph <- graph.adjacency(adjm,mode="undirected",weighted=TRUE,diag=FALSE)
        ## remove negative edges
        net <- delete.edges(net.grph, which(E(net.grph)$weight<0)) # here's my condition.
        ## # to plot this
        ## plot(net, vertex.label="")
        oc <- cluster_fast_greedy(net) # cluster
        ## and now we merge based on the clustered modules
        group <- list()
        for (i in 1:length(levels(as.factor(oc$membership)))){
            group[[i]] <- oc$names[which(oc$membership==i)]
            phyloseq <- merge_taxa(phyloseq, group[[i]],
                                   multi.fn=function(x) paste(x,
                                                              collapse="|"))
        }
    }
    phyloseq
}


getTechRepPaired <- function (ps){
    p <- psmelt(ps)
    pivot_wider(p,
                id_cols = c(OTU, Animal.No, amplicon),
                values_from = Abundance,
                names_from = pool,
                values_fn = ~(sum(.x)))
}


ASVsONETWO <- getTechRepPaired(P)
genusONETWO <- getTechRepPaired(
    tax_glom(P, "genus",
             multi.fn=function(x) paste(x, collapse="|")))
familyONETWO <- getTechRepPaired(
    tax_glom(P, "family",
             multi.fn=function(x) paste(x, collapse="|")))
cASVgONETWO <- getTechRepPaired(combineASVs(P, "genus"))
cASVfONETWO <- getTechRepPaired(combineASVs(P, "family"))
cASVpONETWO <- getTechRepPaired(combineASVs(P, "phylum"))

allPairs <- rbind(cbind(ASVsONETWO, meth="ASVs"),      
                  cbind(genusONETWO, meth="GLOg"),
                  cbind(familyONETWO, meth="GLOf"),
                  cbind(cASVgONETWO, meth="COMg"),
                  cbind(cASVfONETWO, meth="COMf"),
                  cbind(cASVpONETWO, meth="COMp"))

## Hi!
by(allPairs, allPairs$meth, function(x){
    cor(x[, "P_II"], x[, "P_I"], use="pairwise.complete.obs")
})

repDat <- by(allPairs, list(allPairs$OTU, allPairs$meth), function (x) {
    cbind(Cor=cor(x[, "P_I"], x[, "P_II"], use="pairwise.complete.obs"),
          Sum=sum(x[, "P_I"], x[, "P_II"]),
          Prev=sum(x[, "P_I"]>0, x[, "P_II"]>0),
          amplicon=unique(x[, "amplicon"]),
          meth=unique(x[,"meth"]))
})



repDat <- data.frame(do.call("rbind", as.list(repDat)))
rownames(repDat) <- 1:nrow(repDat)

repDat$Namplicon <- sapply(strsplit(repDat$amplicon, "\\|"), length)
repDat <- repDat[!is.na(repDat[, "Cor"]), ]
repDat$Sum <- as.numeric(repDat$Sum)
repDat$Cor <- as.numeric(repDat$Cor)
repDat$Prev <- as.numeric(repDat$Prev)

tapply(repDat$Cor, list(cut(repDat$Sum, 10^(1:6)), repDat$meth), mean, na.rm=TRUE)

tapply(repDat$Cor, list(cut(repDat$Prev, 1:10*20), repDat$meth), mean, na.rm=TRUE)

tapply(repDat$Cor, list(cut(repDat$Namplicon, c(0, 1, 2, 5, 10, 20, 40, 60, Inf)) , repDat$meth),
       mean, na.rm=TRUE)

## plot
pdf("figures/Eval_tech_cor.pdf", width=20, height=8)
ggplot(repDat, aes(Cor, Prev, size=Sum,
                   color=cut(repDat$Namplicon, c(0, 1, 2, 5, 10, 20, 40, 60, Inf)))) +
    geom_point() +
    facet_wrap(~meth)
dev.off()

### super interesting, now tabulate proportion of reads and ASVSs
### "lost" when excluding based on tech-rep correlations

## t
## cASVsGENUS <- combineASVs(P, "genus")

