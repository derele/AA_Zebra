library(tidyverse)
library(MultiAmplicon)
library(Hmisc)
library(Matrix)
library(igraph)

redoBioInf <- FALSE

if(!exists("PS", mode="S4")|redoBioInf){
    source("1_Zebra_general_MA.R")
}

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

Pm <- psmelt(P)
    
ASVsONETWO <- pivot_wider(Pm,
                          id_cols = c(OTU, Animal.No, Date, Sex, Age, Repro,
                                      Season, dens, hab, veg, StrongyEPG,
                                      CryptoEPG, AscaridEPG, AnoploEPG, FGM, 
                                      rawash, ndf, adf, adl, cellul, hemic, energy, N, 
                                      Prot, Ca, Cu, Fe, K, Mg, Mn, Mo, Na, P, S, Zn,
                                      superkingdom, phylum, class, order, family,
                                      genus, species),
            values_from = Abundance,
            names_from = pool,
            values_fn = ~(sum(.x)))
}


cor(ASVsONETWO$P_II, ASVsONETWO$P_I, use="pairwise.complete.obs")

repDat <- by(ASVsONETWO, ASVsONETWO$OTU, function (x) {
    c(Cor=cor(x[, "P_I"], x[, "P_II"]),
      lCor=cor(log10(x[, "P_I"] +1 ), log10(x[, "P_II"] +1)),
      Sum=sum(x[, "P_I"], x[, "P_II"]))
})
    
repDat <- data.frame(do.call("rbind", as.list(repDat)))
rownames(repDat) <- 1:nrow(repDat)

tapply(repDat$Cor,cut(repDat$Sum, 10^(1:6)), mean, na.rm=TRUE)

### interesting to look at!
### But now: does correlation improve when collapsing genera?

PSG <- tax_glom(P, "genus")
PSGm <- psmelt(PSG)
    
genusONETWO  <- pivot_wider(PSGm,
                            id_cols = c(Animal.No, Date, Sex, Age, Repro, Season,
                                        dens, hab, veg, StrongyEPG, CryptoEPG, AscaridEPG,
                                        AnoploEPG, FGM, rawash, ndf, adf, adl, cellul,
                                        hemic, energy, N, Prot, Ca, Cu, Fe, K, Mg, Mn,
                                        Mo, Na, P, S, Zn, superkingdom, phylum, class,
                                        order, family, genus),
                    values_from = Abundance,
                    names_from = pool,
                    values_fn = ~(sum(.x)))

cor(genusONETWO$P_II, genusONETWO$P_I, use="pairwise.complete.obs")


PSF <- tax_glom(P, "family")


PSFm <- psmelt(PSF)
    
bar  <- pivot_wider(PSFm,
                    id_cols = c(Animal.No, Date, Sex, Age, Repro, Season, dens, hab,
                                veg, StrongyEPG, CryptoEPG, AscaridEPG, AnoploEPG, FGM, 
                                rawash, ndf, adf, adl, cellul, hemic, energy, N, 
                                Prot, Ca, Cu, Fe, K, Mg, Mn, Mo, Na, P, S, Zn,
                                superkingdom, phylum, class, order, family),
                    values_from = Abundance,
                    names_from = pool,
                    values_fn = ~(sum(.x)))

cor(bar$P_II, bar$P_I, use="pairwise.complete.obs")



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
            phyloseq <- merge_taxa(phyloseq, group[[i]])
        }
    }
    phyloseq
}

PcASV <- combineASVs(P)

PcASVm <- psmelt(PcASV)
    
baz  <- pivot_wider(PcASVm,
                    id_cols = c(Animal.No, Date, Sex, Age, Repro, Season, dens, hab,
                                veg, StrongyEPG, CryptoEPG, AscaridEPG, AnoploEPG, FGM, 
                                rawash, ndf, adf, adl, cellul, hemic, energy, N, 
                                Prot, Ca, Cu, Fe, K, Mg, Mn, Mo, Na, P, S, Zn,
                                superkingdom, phylum, class, order, family, OTU),
                    values_from = Abundance,
                    names_from = pool,
                    values_fn = ~(sum(.x)))

cor(baz$P_II, baz$P_I, use="pairwise.complete.obs")



PcASVF <- combineASVs(P, "family")
PcASVFm <- psmelt(PcASVF)
    
boom  <- pivot_wider(PcASVFm,
                    id_cols = c(Animal.No, Date, Sex, Age, Repro, Season, dens, hab,
                                veg, StrongyEPG, CryptoEPG, AscaridEPG, AnoploEPG, FGM, 
                                rawash, ndf, adf, adl, cellul, hemic, energy, N, 
                                Prot, Ca, Cu, Fe, K, Mg, Mn, Mo, Na, P, S, Zn,
                                superkingdom, phylum, class, order, family, OTU),
                    values_from = Abundance,
                    names_from = pool,
                    values_fn = ~(sum(.x)))

cor(boom$P_II, boom$P_I, use="pairwise.complete.obs")



### going extreme... relying only on correlations...
PcASVULTRA <- combineASVs(P, "phylum")
PcASVULTRAm <- psmelt(PcASVULTRA)
    
bazooka  <- pivot_wider(PcASVULTRAm,
                    id_cols = c(Animal.No, Date, Sex, Age, Repro, Season, dens, hab,
                                veg, StrongyEPG, CryptoEPG, AscaridEPG, AnoploEPG, FGM, 
                                rawash, ndf, adf, adl, cellul, hemic, energy, N, 
                                Prot, Ca, Cu, Fe, K, Mg, Mn, Mo, Na, P, S, Zn,
                                superkingdom, phylum, OTU),
                    values_from = Abundance,
                    names_from = pool,
                    values_fn = ~(sum(.x)))

cor(bazooka$P_II, bazooka$P_I, use="pairwise.complete.obs")
