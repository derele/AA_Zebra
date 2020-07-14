library(vegan)

## need to see what's required and source the respective data...


getOrdinations <- function(ps){
    ps <- prune_taxa(taxa_sums(ps) > 0, ps)
    ## at least 100 otus in table otherwise ordination breaks
    if(ncol(otu_table(ps))<100){
        return(NULL)
    } else{
        psD <- sample_data(ps)
        psNA <- subset_samples(ps, !is.na(cellul)&
                                   !is.na(FGM))
        ps.dat <- sample_data(psNA)
        class(ps.dat) <- "data.frame"
        ## dbrda because it can do maringals  
        adonis <- dbrda(otu_table(psNA)~veg+Repro+dens,
                        data=ps.dat, distance="bray")
        marginals <- anova(adonis, by = 'margin')
        marginals
    }
}

allO <- getOrdinations(ps=PS)

allO <- getOrdinations(ps=PM)

allO <- getOrdinations(ps=PGM)

allO <- getOrdinations(ps=PMharsh)

allO <- getOrdinations(ps=PGMharsh)

allO <- getOrdinations(ps=PMGharsh)

EukO <- getOrdinations(subset_taxa(PMGharsh,
                                   superkingdom%in%"Eukaryota"))
EukO[[1]]


BacO <- getOrdinations(subset_taxa(prune_taxa(taxa_sums(PMGharsh) > 5, PMGharsh),
                       superkingdom%in%"Bacteria"))

FunO <- getOrdinations(subset_taxa(PGMharsh,
                                   phylum%in%c("Ascomycota",
                                               "Basidiomycota",
                                               "Blastocladiomycota",
                                               "Chytridiomycota",
                                               "Cryptomycota",
                                               "Mucoromycota",
                                               "Olpidiomycota",
                                               "Planctomycetes",
                                               "Zoopagomycota")))
FunO

ord.l <- lapply(PS.l, getOrdinations)
names(ord.l) <- names(PS.l)

ordM.l <- lapply(PM.l, getOrdinations)
names(ord.l) <- names(PS.l)


### ord.l[["ADM330.Klin0785_CR 16S"]] ## has significant density

## "NMDS", "bray" does not work ... PCoA works
getOrdPlot <- function (ps, color="dens", shape="veg"){
    psp <- prune_samples(sample_sums(ps)>=20, ps)
    ord.unsup <- ordinate(psp, "PCoA", "bray")
    plot.unsup <- plot_ordination(psp, ord.unsup,
                                  color=color, shape=shape) +
        geom_point(size=5) 
    psD <- sample_data(psp)
    psNA <- subset_samples(psp, !is.na(cellul)&
                                !is.na(Age)&
                                !is.na(dens))
    formula <- as.formula(paste0("psNA ~",  color,  "+",  shape))
    ord.sup <- ordinate(psNA,
                        "CCA", formula = formula)
    plot.sup <- plot_ordination(psNA, ord.sup,
                                color=color, shape=shape) +
        geom_point(size=5) 
    list(ord.unsup, plot.unsup, ord.sup, plot.sup)
}

ordPlot <- getOrdPlot(PM.l[["ACM_008.Klin0341_CR 16S"]], "Repro")

ordPlot <- getOrdPlot(PM.l[["Klin0341.Klin0785_CR 16S"]], "Repro")

ordPlot <- getOrdPlot(PM.l[["ADM330.Klin0785_CR 16S"]], "Repro")

ordPlot <- getOrdPlot(PS.l[["Proti15.Proti440R 18S"]], "Sex")

ordPlot <- getOrdPlot(PS.l[["D3A_5_3Mod.D3B_5_3Mod 28S"]], "Repro")

dev.off()

ordPlot <- getOrdPlot(PGMharsh)


## 18S: no comositional difference between sexes!?

pdf("figures/ordi_18S_raw_brayPCoA_vegcol.pdf")
plot_ordination(ps18S, X18S.ord, color="veg", shape="Season") +
        geom_point(size=5) 
dev.off()

## 18S: high comositional difference between vegetation!?

pdf("figures/ordi_18S_raw_brayPCoA_seascol.pdf")
plot_ordination(ps18S, X18S.ord, color="Season", shape="veg") +
    geom_point(size=5) 
dev.off()

## 18S: high compositional difference between seasons, maybe less than
## 16S!?

pdf("figures/ordi_18S_raw_brayPCoA_cellulcol.pdf")
plot_ordination(ps18S, X18S.ord, color="cellul", shape="Season") +
    geom_point(size=5) 
dev.off()

## 18S: composition varying with cellulose content in diet!?


ps18sdat.cellul <- ps18sdat[!is.na(ps18sdat$cellul), ]


############################################# 16 S ###############
ps16S <- phyloseq(otu_table(X16S, taxa_are_rows=FALSE),
                  sample_data(samples.final))


X16S.ord <- ordinate(ps16S, "PCoA", "bray")

pdf("figures/ordi_16S_raw_brayPCoA_sexcol.pdf")
plot_ordination(ps16S, X16S.ord, color="Sex", shape="Season")+
geom_point(size=5)
dev.off()

## 16S: no comositional difference between sexes!?

pdf("figures/ordi_16S_raw_brayPCoA_vegcol.pdf")
plot_ordination(ps16S, X16S.ord, color="veg", shape="Season") +
    geom_point(size=5) 
dev.off()

## 16S: high comositional difference between vegetation!?

ano.16S.veg <- anosim(otu_table(ps16S), sample_data(ps16S)$veg)

ps16sdat <- sample_data(ps16S)
class(ps16sdat) <- "data.frame"

pdf("figures/ordi_16S_raw_brayPCoA_seascol.pdf")
plot_ordination(ps16S, X16S.ord, color="Season", shape="veg") +
    geom_point(size=5) 
dev.off()

## 16S: high compositional difference between seasons!?

pdf("figures/ordi_16S_raw_brayPCoA_cellulcol.pdf")
plot_ordination(ps16S, X16S.ord, color="cellul", shape="Season") +
    geom_point(size=5) 
dev.off()

## 16S: composition varying with cellulose content in diet!?

ps16S.cellul <- subset_samples(ps16S, !is.na(sample_data(ps16S)$cellul))

ps16sdat.cellul <- ps16sdat[!is.na(ps16sdat$cellul), ]

###  Overall
##https://stats.stackexchange.com/questions/312302/adonis-in-vegan-order-of-variables-non-nested-with-one-degree-of-freedom-for

## dbrda because it can do maringals  
ado.16S.all <- dbrda(otu_table(ps16S.cellul)~cellul+Season+veg+
                     Sex+Age+Repro+dens,
                     data=ps16sdat.cellul)

marg.ado.16S.all <- anova(ado.16S.all, by = 'margin')

########################### 28S ##############################
ps28S <- phyloseq(otu_table(X28S, taxa_are_rows=FALSE),
                  sample_data(sample.data[rownames(X28S), ]),
                  tax_table(as.matrix(X28S.tax)))

ps28S.genus.own <- tax_glom(ps28S, "Genus_own")

Gen.tab.28S <- otu_table(subset_taxa(ps28S.genus.own))

genus.names28 <- tax_table(ps28S.genus.own)[, "Genus_own"]

gen.c.l <- by(t(Gen.tab.28S), genus.names28, colSums) 

Gen.tab.28S <- t(do.call(rbind, gen.c.l))

round(colSums(Gen.tab.28S)[order(colSums(Gen.tab.28S))])

Gen.tab.28S <- as.data.frame(merge(Gen.tab.28S, sample.data, by=0))

get.best.cors <- function(column, cutoff) {
    foo <- cor(Gen.tab.28S[, unique(genus.names28)], Gen.tab.28S[, column],
               use="na.or.complete", method="spearman")
    rows <- apply(foo, 1, function (x) any(abs(x)>cutoff))
    foo[which(rows), ]
}


get.best.cors("StrongyEPG", 0.1)[order(get.best.cors("StrongyEPG", 0.1))]

pdf("figures/Strongy28SCor.pdf")
ggplot(Gen.tab.28S, aes(Petrovinema + Cylicocyclus, StrongyEPG)) +
    geom_point() +
    scale_y_log10() + scale_x_log10() +
    geom_smooth(method="lm", se=FALSE) +
    stat_smooth_func(geom="text", method="lm", hjust=0, parse=TRUE) +
    theme_bw()
dev.off()


X28S.ord <- ordinate(ps28S, "PCoA", "bray")

pdf("figures/ordi_28S_raw_brayPCoA_sexcol.pdf")
plot_ordination(ps28S, X28S.ord, color="Sex", shape="Season")+
geom_point(size=5)
dev.off()

## 28S: no comositional difference between sexes!?

pdf("figures/ordi_28S_raw_brayPCoA_vegcol.pdf")
plot_ordination(ps28S, X28S.ord, color="veg", shape="Season") +
    geom_point(size=5) 
dev.off()

## 28S: comositional difference between vegetation (but less than
## 16S)!?

ano.28S.veg <- anosim(otu_table(ps28S), sample_data(ps28S)$veg)

ps28sdat <- sample_data(ps28S)
class(ps28sdat) <- "data.frame"

pdf("figures/ordi_28S_raw_brayPCoA_seascol.pdf")
plot_ordination(ps28S, X28S.ord, color="Season", shape="veg") +
    geom_point(size=5) 
dev.off()

## 28S: not a strong comositional difference between seasons?!

pdf("figures/ordi_28S_raw_brayPCoA_cellulcol.pdf")
plot_ordination(ps28S, X28S.ord, color="cellul", shape="Season") +
    geom_point(size=5) 
dev.off()

ps28S.cellul <- subset_samples(ps28S, !is.na(sample_data(ps28S)$cellul))

ano.28S.cell <- anosim(otu_table(ps28S.cellul),
                       sample_data(ps28S.cellul)$cellul)


ps28sdat.cellul <- ps28sdat[!is.na(ps28sdat$cellul), ]

###  Overall
##https://stats.stackexchange.com/questions/312302/adonis-in-vegan-order-of-variables-non-nested-with-one-degree-of-freedom-for

## dbrda because it can do maringals  
ado.28S.all <- dbrda(otu_table(ps28S.cellul)~cellul+Season+veg+
                     Sex+Age+Repro+dens, ###+FGM,
                     data=ps28sdat.cellul)

anova(ado.28S.all, by = 'margin')

## 28S: composition varying with cellulose content in diet!?


psALL <- phyloseq(otu_table(ALL, taxa_are_rows=FALSE),
                  sample_data(samples.final))

