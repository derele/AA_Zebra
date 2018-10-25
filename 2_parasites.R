library(pheatmap)
library(DESeq2)
library(ggplot2)


if(!exists("PM")){
    source("1_general_MA.R")
}


################################## 18S and 28S  #########################
PM.para <- subset_taxa(PM, phylum%in%c("Nematoda",
                                       "Apicomplexa",
                                       "Platyhelminthes"))

PG.para <- tax_glom(PM.para, "genus", NArm = TRUE)

table(tax_table(PG.para)[,6])

mat <- as.matrix(t(otu_table(PG.para)))
rownames(mat) <- make.names(tax_table(PG.para)[, 6])

colanot <- data.frame(lapply(sample_data(PG.para)[, c("veg", "Sex", "Age")],
                             factor))
ann_colors <- list(
    veg = c("orange", "green", "darkgreen"),
    Sex = c("white", "red", "blue"),
    Age = c("blue", "red", "green", "orange")
)


pdf("figures/parasite_genera_heat.pdf", width=8, height=14)
pheatmap(log10(mat+1),  labels_col=rep("", times=ncol(mat)),
         annotation_col = colanot)
dev.off()


## high vervalence parasites
hPP <- c("Lamanema", "Strongylus", "Necator", "Cylicostephanus",
         "Gregarina", "Strongyloides", "Dictyocaulus", "Haemonchus", "Oxyuris")

pdf("figures/parasite_Abu_genera_heat.pdf", width=8, height=4)
pheatmap(log10(mat[hPP,]+1),  labels_col=rep("", times=ncol(mat[hPP,])),
         annotation_col = colanot)
dev.off()


## By ASV for prevalent genera
PM.hPP <- subset_taxa(PM, genus%in%hPP)

PM.hPP <- subset_taxa(PM.hPP,
                      ## more than 100 as count
                      colSums(otu_table(PM.hPP))>100 &
                      ## and in more than 5 samples
                      colSums(otu_table(PM.hPP)>0)>5)
PM.hPP.Ages <- subset_samples(PM.hPP,
                              Age%in%c("fl", "juv", "mat", "sa") &
                              Season%in%c("d", "w") &
                              hab%in%c("lg", "sgp"))


deseqTestNULL <- function(P, full, red, alpha=0.1){
    diagdds <- phyloseq_to_deseq2(P, full)
    diagdds <- estimateSizeFactors(diagdds, type="poscounts") 
    diagdds.red <- DESeq(diagdds, test="LRT", reduced = red, fitType="parametric")
    res.red <- results(diagdds.red, cooksCutoff = FALSE)
    sigtab.red <- res.red[which(res.red$padj < alpha), ]
    if (nrow(sigtab.red)>0){
        sigtab.red <- cbind(as(sigtab.red, "data.frame"),
                            as(tax_table(P)[rownames(sigtab.red), ], "matrix"))
        rownames(sigtab.red) <- NULL
        sigtab.red
    } else NULL
}


deseqTestNULL(PM.hPP.Ages, ~ Age*Season, ~Age+Season)
deseqTestNULL(PM.hPP.Ages, ~ Age*Season, ~Age)
deseqTestNULL(PM.hPP.Ages, ~ Age*Season, ~Season)
deseqTestNULL(PM.hPP.Ages, ~ Age*Season, ~1)
deseqTestNULL(PM.hPP.Ages, ~ Age, ~1)
deseqTestNULL(PM.hPP.Ages, ~ Season, ~1)
## absolutely nothing!!!


## By Genus for prevalent genera
PG.hPP <- subset_taxa(PG.para, genus%in%hPP)

PG.hPP.Ages <- subset_samples(PG.hPP,
                              Age%in%c("fl", "juv", "mat", "sa") &
                              Season%in%c("d", "w") &
                              hab%in%c("lg", "sgp"))

deseqTestNULL(PG.hPP.Ages, ~ Age*Season, ~Age+Season)
deseqTestNULL(PG.hPP.Ages, ~ Age*Season, ~Age)
deseqTestNULL(PG.hPP.Ages, ~ Age*Season, ~Season)
deseqTestNULL(PG.hPP.Ages, ~ Age*Season, ~1)
deseqTestNULL(PG.hPP.Ages, ~ Age, ~1)
deseqTestNULL(PG.hPP.Ages, ~ Season, ~1)
## Age vs the null model is somewhat significant for gGregarina!!!
## This might credible and interesting Gregarines more in foals,
## juveniles and subadults??

Gregarina <- unname(otu_table(subset_taxa(PG.hPP.Ages, genus%in%"Gregarina")))
## need to sum up because of alternative taxonomic paths
Gregarina <- rowSums(Gregarina)

Gregarina <- cbind(sample_data(PG.hPP.Ages), Gregarina)

sort(tapply(Gregarina$Gregarina, as.factor(Gregarina$Season):as.factor(Gregarina$Age), median))
sort(tapply(Gregarina$Gregarina, as.factor(Gregarina$Season):as.factor(Gregarina$Age), function(x) mean(log10(x+1))))

sort(tapply(Gregarina$Gregarina, as.factor(Gregarina$Age), median))
sort(tapply(Gregarina$Gregarina, as.factor(Gregarina$Age), function(x) mean(log10(x+1))))
## Probably NOT! Very dodgy...


## Testing for reproductive status and season, ONLY for adult females!
PM.hPP.AdultF <- subset_samples(PM.hPP,
                                   Age%in%"mat" &
                                   Repro%in%c( "lact", "preg", "u") &
                                   Sex%in%c("f") &
                                   Season%in%c("d", "w") &
                                   hab%in%c("lg", "sgp")
                                )

deseqTestNULL(PM.hPP.AdultF, ~ Repro*Season, ~Repro+Season)
deseqTestNULL(PM.hPP.AdultF, ~ Repro*Season, ~Repro)
deseqTestNULL(PM.hPP.AdultF, ~ Repro*Season, ~Season)
deseqTestNULL(PM.hPP.AdultF, ~ Repro*Season, ~1)
deseqTestNULL(PM.hPP.AdultF, ~ Repro, ~1)
deseqTestNULL(PM.hPP.AdultF, ~ Season, ~1)
## Absolutely nothing!

## Same for genera
PG.hPP.AdultF <- subset_samples(PG.hPP,
                                Age%in%"mat" &
                                Repro%in%c("lact", "preg", "u") &
                                Sex%in%c("f") &
                                Season%in%c("d", "w") &
                                hab%in%c("lg", "sgp")
                                )

deseqTestNULL(PG.hPP.AdultF, ~ Repro*Season, ~Repro+Season)
deseqTestNULL(PG.hPP.AdultF, ~ Repro*Season, ~Repro)
deseqTestNULL(PG.hPP.AdultF, ~ Repro*Season, ~Season)
deseqTestNULL(PG.hPP.AdultF, ~ Repro*Season, ~1)
deseqTestNULL(PG.hPP.AdultF, ~ Repro, ~1)
deseqTestNULL(PG.hPP.AdultF, ~ Season, ~1)
## Absolutely nothing!


## Testing for reproductive status and season, ONLY for adult females!
PM.hPP.AdultM <- subset_samples(PM.hPP,
                                   Age%in%"mat" &
                                   Repro%in%c("b", "h") &
                                   Sex%in%c("m") &
                                   Season%in%c("d", "w") &
                                   hab%in%c("lg", "sgp")
                                )

deseqTestNULL(PM.hPP.AdultM, ~ Repro*Season, ~Repro+Season)
deseqTestNULL(PM.hPP.AdultM, ~ Repro*Season, ~Repro)
deseqTestNULL(PM.hPP.AdultM, ~ Repro*Season, ~Season)
deseqTestNULL(PM.hPP.AdultM, ~ Repro*Season, ~1)
deseqTestNULL(PM.hPP.AdultM, ~ Repro, ~1)
deseqTestNULL(PM.hPP.AdultM, ~ Season, ~1)
## Absolutely nothing!

## Same for genera
PG.hPP.AdultM <- subset_samples(PG.hPP,
                                Age%in%"mat" &
                                Repro%in%c("b", "h") &
                                Sex%in%c("m") &
                                Season%in%c("d", "w") &
                                hab%in%c("lg", "sgp")
                                )

deseqTestNULL(PG.hPP.AdultM, ~ Repro*Season, ~Repro+Season)
deseqTestNULL(PG.hPP.AdultM, ~ Repro*Season, ~Repro)
deseqTestNULL(PG.hPP.AdultM, ~ Repro*Season, ~Season)
deseqTestNULL(PG.hPP.AdultM, ~ Repro*Season, ~1)
deseqTestNULL(PG.hPP.AdultM, ~ Repro, ~1)
deseqTestNULL(PG.hPP.AdultM, ~ Season, ~1)

## more Haemonchus in harem stallions?
Haemonchus <- unname(otu_table(subset_taxa(PG.hPP.AdultM, genus%in%"Haemonchus")))
## need to sum up because of alternative taxonomic paths
Haemonchus <- rowSums(Haemonchus)

Haemonchus <- cbind(sample_data(PG.hPP.AdultM), Haemonchus)

sort(tapply(Haemonchus$Haemonchus, as.factor(Haemonchus$Season):as.factor(Haemonchus$Repro), median))
sort(tapply(Haemonchus$Haemonchus, as.factor(Haemonchus$Season):as.factor(Haemonchus$Repro), function(x) mean(log10(x+1))))

sort(tapply(Haemonchus$Haemonchus, as.factor(Haemonchus$Repro), median))
sort(tapply(Haemonchus$Haemonchus, as.factor(Haemonchus$Repro), function(x) mean(log10(x+1))))

## Bachelors have more!!!

## ggplot(Haemonchus, aes(Repro, Haemonchus)) +
##     geom_boxplot() +
##     scale_y_log10()
       

### Extract and revisit all Necator ASVs
Necator <- unname(otu_table(subset_taxa(PM.hPP.Ages, genus%in%"Necator")))

Necator.seq <- colnames(otu_table(subset_taxa(PM.hPP.Ages, genus%in%"Necator")))
Necator.seq <- unlist(strsplit(Necator.seq, "NNNNNNNNNN"))

writeFasta(Necator.seq, "Necator.fasta")


### Extract and revisit all Cylicostephanus ASVs
Cylico <- unname(otu_table(subset_taxa(PM.hPP.Ages, genus%in%"Cylicostephanus")))

Cylico.seq <- colnames(otu_table(subset_taxa(PM.hPP.Ages, genus%in%"Cylicostephanus")))
Cylico.seq <- unlist(strsplit(Cylico.seq, "NNNNNNNNNN"))

writeFasta(Cylico.seq, "Cylico.fasta")



