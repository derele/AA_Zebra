if(!exists("PS")){
    source("1_general_MA.R")
}


################################## 18S and 28S  ##########################
## library(DESeq2)
PS.para <- subset_taxa(PS, phylum%in%c("Nematoda",
                                       "Apicomplexa",
                                       "Platyhelminthes"))

PG.para <- tax_glom(PS.para, "genus", NArm = TRUE)

table(tax_table(PG.para)[,6])

Zeb <- by(otu_table(PG.para), sample_data(PG.para)$Animal.No, colSums)
Zebra <- do.call(cbind, Zeb)
Zebra <- Zebra[, !colnames(Zebra)%in%"neg"]

PSM.para <- phyloseq(otu_table(Zebra, taxa_are_rows=TRUE),
                     sample_data(sample.data[colnames(Zebra), ]),
                     tax_table(tax_table(PG.para)))

PSM.para.Ages <- subset_samples(PSM.para,
                                Age%in%c("fl", "juv", "mat", "sa") &
                                Season%in%c("d", "w") &
                                hab%in%c("lg", "sgp"))

library(pheatmap)

mat <- as.matrix(otu_table(PSM.para.Ages))
rownames(mat) <- make.names(tax_table(PSM.para.Ages)[, 6])

colanot <- data.frame(lapply(sample_data(PSM.para.Ages)[, c("veg", "Sex", "Age")],
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




diagdds <- phyloseq_to_deseq2(PSM.para.Ages, ~ Age*Season)
diagdds <- estimateSizeFactors(diagdds, type="poscounts") # type="iterate"
diagdds.interact <- DESeq(diagdds, test="LRT",
                          reduced = ~ Age+Season, fitType="parametric")
diagdds.Season <- DESeq(diagdds, test="LRT", reduced = ~ Age, fitType="parametric")
diagdds.Age <- DESeq(diagdds, test="LRT", reduced = ~ Season, fitType="parametric")
diagdds.null <- DESeq(diagdds, test="LRT", reduced = ~ 1, fitType="parametric")


res.interact <- results(diagdds.interact, cooksCutoff = FALSE)
res.Season <- results(diagdds.Season, cooksCutoff = FALSE)
res.Age <- results(diagdds.Age, cooksCutoff = FALSE)
res.null <- results(diagdds.null, cooksCutoff = FALSE)

alpha <- 0.5

sigtab.interact <- res.interact[which(res.interact$padj < alpha), ]
sigtab.interact <- cbind(as(sigtab.interact, "data.frame"),
                as(tax_table(PG.para.Ages)[rownames(sigtab.interact), ], "matrix"))
rownames(sigtab.interact) <- NULL


sigtab.Age <- res.Age[which(res.Age$padj < alpha), ]
sigtab.Age <- cbind(as(sigtab.Age, "data.frame"),
                as(tax_table(PG.para.Ages)[rownames(sigtab.Age), ], "matrix"))
rownames(sigtab.Age) <- NULL

sigtab.Season <- res.Season[which(res.Season$padj < alpha), ]
sigtab.Season <- cbind(as(sigtab.Season, "data.frame"),
                as(tax_table(PG.para.Ages)[rownames(sigtab.Season), ], "matrix"))
rownames(sigtab.Season) <- NULL


sigtab.null <- res.null[which(res.null$padj < alpha), ]
sigtab.null <- cbind(as(sigtab.null, "data.frame"),
                as(tax_table(PG.para.Ages)[rownames(sigtab.null), ], "matrix"))
rownames(sigtab.null) <- NULL

### An interaction effect of Season and Age on Ostertagia and Cylicoccus
foo <- cbind(sample_data(PSM),
             t(unname(otu_table(subset_taxa(PSM, Genus%in%"Cylicocyclus")))))

sort(tapply(foo$sp1, as.factor(foo$Season):as.factor(foo$Age), median))

PSM.para.AdultF <- subset_samples(PSM,
                                  Age%in%"mat" &
                                  Repro%in%c("b", "h", "lact", "preg", "u") &
                                  Sex%in%c("f") &
                                  Season%in%c("d", "w") &
                                  hab%in%c("lg", "sgp")
                                  )


diagdds <- phyloseq_to_deseq2(PSM.para.AdultF, ~ Repro*Season)
diagdds <- estimateSizeFactors(diagdds, type="poscounts") # type="iterate"
diagdds.interact <- DESeq(diagdds, test="LRT",
                          reduced = ~ Repro+Season, fitType="parametric")

diagdds.Season <- DESeq(diagdds, test="LRT", reduced = ~ Repro, fitType="parametric")
diagdds.Repro <- DESeq(diagdds, test="LRT", reduced = ~ Season, fitType="parametric")
diagdds.null <- DESeq(diagdds, test="LRT", reduced = ~ 1, fitType="parametric")


res.interact <- results(diagdds.interact, cooksCutoff = FALSE)
res.Season <- results(diagdds.Season, cooksCutoff = FALSE)
res.Repro <- results(diagdds.Repro, cooksCutoff = FALSE)
res.null <- results(diagdds.null, cooksCutoff = FALSE)

alpha <- 0.5

sigtab.interact <- res.interact[which(res.interact$padj < alpha), ]
sigtab.Repro <- res.Repro[which(res.Repro$padj < alpha), ]
sigtab.Season <- res.Season[which(res.Season$padj < alpha), ]
sigtab.null <- res.null[which(res.null$padj < alpha), ]
## empty

PSM.para.AdultM <- subset_samples(PSM,
                                 Age%in%"mat" &
                                 Repro%in%c("b", "h", "lact", "preg", "u") &
                                 Sex%in%c("m") &
                                 Season%in%c("d", "w") &
                                 hab%in%c("lg", "sgp")
                                 )


diagdds <- phyloseq_to_deseq2(PSM.para.AdultM, ~ Repro*Season)
diagdds <- estimateSizeFactors(diagdds, type="poscounts") # type="iterate"
diagdds.interact <- DESeq(diagdds, test="LRT",
                          reduced = ~ Repro+Season, fitType="parametric")

diagdds.Season <- DESeq(diagdds, test="LRT", reduced = ~ Repro, fitType="parametric")
diagdds.Repro <- DESeq(diagdds, test="LRT", reduced = ~ Season, fitType="parametric")
diagdds.null <- DESeq(diagdds, test="LRT", reduced = ~ 1, fitType="parametric")


res.interact <- results(diagdds.interact, cooksCutoff = FALSE)
res.Season <- results(diagdds.Season, cooksCutoff = FALSE)
res.Repro <- results(diagdds.Repro, cooksCutoff = FALSE)
res.null <- results(diagdds.null, cooksCutoff = FALSE)

alpha <- 0.5

sigtab.interact <- res.interact[which(res.interact$padj < alpha), ]
sigtab.Repro <- res.Repro[which(res.Repro$padj < alpha), ]
sigtab.Season <- res.Season[which(res.Season$padj < alpha), ]
sigtab.null <- res.null[which(res.null$padj < alpha), ]

library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}
## Phylum order
x <- tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x <- sort(x, TRUE)
sigtab$Phylum <- factor(as.character(sigtab$Phylum), levels=names(x))

## Genus order
x <- tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x <- sort(x, TRUE)
sigtab$Genus <- factor(as.character(sigtab$Genus), levels=names(x))

ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) +
    geom_point(size=6) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))



