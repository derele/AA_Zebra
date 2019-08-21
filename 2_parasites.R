library(pheatmap)
library(DESeq2)
library(ggplot2)


if(!exists("PM")){
    source("1_general_MA.R")
}


################################## 18S and 28S  #########################
PS.para <- subset_taxa(PS, phylum%in%c("Nematoda",
                                       "Apicomplexa",
                                       "Platyhelminthes"))

PM.para <- subset_taxa(PM, phylum%in%c("Nematoda",
                                       "Apicomplexa",
                                       "Platyhelminthes"))

ordi.l <- lapply(PM.l, function (P) { 
    GP <- prune_taxa(taxa_sums(P) > 0, P)
    GP <- prune_samples(sample_sums(GP) > 0, GP)
    plot_ordination(GP, ordinate(GP, "MDS"), color = "Season") + geom_point(size = 5)
})

names(ordi.l) <- names(PS.l)

ordi.l <- ordi.l[!unlist(lapply(ordi.l, is.null))]

lapply(seq_along(ordi.l), function (i) {
    ordi.l[[i]] + ggtitle(names(ordi.l[i]))
})

GP <- prune_taxa(taxa_sums(PM) > 0, PS)

plot_richness(GP, x="Repro", measures="Observed")  + geom_violin() + geom_jitter() 
dev.off()


foo <- tax_glom(subset_taxa(PM, genus%in%"Toxocara"), "genus")
bar <- as(otu_table(foo), "vector")
foobar <- cbind(bar, sample_data(foo))



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

## high prervalence Strongyles!!!
hPP <- c("Strongylus", "Cylicostephanus", "Cylicocyclus", "Caenorhabditis")

nameOtuByTax <- function(ps, taxon="genus"){
    otable <- otu_table(ps)
    colnames(otable) <- make.unique(as.character(tax_table(ps)[, "genus"]))
    otable
}

## A correlation analysis to figure out what these "Genera are

#### Approaches to correlation
foo <- lapply(PS.l[c(27, 28)], function (x) subset_taxa(x, genus%in%hPP))

r1 <- rownames(otu_table(foo[[1]]))

bar <- lapply(PM.l[c(27, 28)], function (x) subset_taxa(x, genus%in%hPP))

CCC <- cor(nameOtuByTax(bar[[1]]), nameOtuByTax(bar[[2]]))

no.na.row <- apply(CCC, 1, function (x) !any(is.na(x)))

### pheatmap(CCC, cluster_rows=F, cluster_cols=F)

## An alignment approach

ref18S <- readDNAStringSet("/SAN/db/RDP/annotated_18S_ena.fasta")
ref28 <- readDNAStringSet("/SAN/db/RDP/annotated_28S_ena.fasta")

is.taxon <- sapply(hPP, function (x){
    grepl(x, names(ref18S))
})

is.long <- width(ref18S)>1600

names(ref18S[rowSums(is.long&is.taxon)>0])

sel.seq <- ref18S[rowSums(is.long&is.taxon)>0][1:12]

names(sel.seq) <- make.unique(gsub(".*;(.*)$", "\\1", names(sel.seq)))

PS.l.18S <- PS.l[grepl("18S$", names(PS.l))]

## Just to test #### - > Why does this not work while
hPP <- "Cylicostephanus"

## This (or any other genus) seems to be less problematic????
## hPP <- "Cylicocyclus"

## hPP <- "Entamoeba"

PS.l.18S.hpp <- lapply(PS.l.18S, function(x) {
    if(nrow(tax_table(x))>0 &
       any(tax_table(x)[, "genus"]%in%hPP)) {
        subset_taxa(x, genus%in%hPP)
    }
})

PS.l.18S.hpp <- PS.l.18S.hpp[(unlist(lapply(PS.l.18S.hpp, function (x) !is.null(x))))]

hpp.seq.18S <- lapply(seq_along(PS.l.18S.hpp), function (i) {
    seq <- colnames(otu_table(PS.l.18S.hpp[[i]]))
    hpp.names <- colnames(nameOtuByTax(PS.l.18S.hpp[[i]]))
    preN <- paste0("A", i, hpp.names)
    ppreN<- gsub("phanus|Cylico|enorhabditis|rongylus|clus", "", preN)
    names(seq) <- gsub("\\.", "_", preN)
    seq
})



## foo <- subset_taxa(PS, genus%in%hPP)

## seq <- colnames(otu_table(foo))
## hpp.names <- colnames(nameOtuByTax(foo))
## preN <- paste0("A", hpp.names)
## ## ppreN<- gsub("phanus|Cylico|enorhabditis|rongylus|clus", "", preN)
## names(seq) <- gsub("\\.", "_", preN)
## seq


sp.hpp.seq.18S <- lapply(hpp.seq.18S, strsplit, "NNNNNNNNNN")

## sp.hpp.seq.18S <- strsplit(seq, "NNNNNNNNNN")

OMGseq <- unlist(sp.hpp.seq.18S)

OMGseq <- OMGseq[!duplicated(OMGseq)]

set.seed(129)
selseq <- OMGseq[sample(1:length(OMGseq), 20)]

writeFasta(c(sel.seq, selseq), "strongyle.fasta" )

writeFasta(sel.seq, "ref.fasta" )

writeFasta(selseq, "query.fasta" )

## high prevalence parasites
hPP_more <- c("Lamanema", "Strongylus", "Necator", "Cylicostephanus",
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
PGM.hPP <- subset_taxa(PGM, genus%in%hPP)

PGM.hPP.Ages <- subset_samples(PGM.hPP,
                              Age%in%c("fl", "juv", "mat", "sa") &
                              Season%in%c("d", "w") &
                              hab%in%c("lg", "sgp"))

deseqTestNULL(PGM.hPP.Ages, ~ Age*Season, ~Age+Season)
deseqTestNULL(PGM.hPP.Ages, ~ Age*Season, ~Age)
deseqTestNULL(PGM.hPP.Ages, ~ Age*Season, ~Season)
deseqTestNULL(PGM.hPP.Ages, ~ Age*Season, ~1)
deseqTestNULL(PGM.hPP.Ages, ~ Age, ~1)
deseqTestNULL(PGM.hPP.Ages, ~ Season, ~1)
## Age vs the null model is somewhat significant for gGregarina!!!
## This might credible and interesting Gregarines more in foals,
## juveniles and subadults??

Gregarina <- unname(otu_table(subset_taxa(PGM.hPP.Ages, genus%in%"Gregarina")))
## need to sum up because of alternative taxonomic paths
Gregarina <- rowSums(Gregarina)

Gregarina <- cbind(sample_data(PGM.hPP.Ages), Gregarina)

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



