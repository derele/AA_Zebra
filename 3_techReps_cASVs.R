
## This should acutally work with phyloseq's merge_samples function
## but doesn't as this messes up sample_data.
## CANDIDATE FOR INCLUSION IN PACKAGE...!

sumTecRep <- function (PS, by.sample, fun=sum){
    otab <- setDT(apply(otu_table(PS), 2, as.list))
    ## the columns giving numbers for sequences
    numcols <- colnames(otab)[nchar(colnames(otab))>10]
    sdat <- sample_data(PS, errorIfNULL = FALSE)
    otab[, (numcols):=lapply(.SD, as.numeric), .SDcols=numcols]
    otab[, sfac := as.factor(sdat[[by.sample]])]
    setkey(otab, sfac)
    otabN <- otab[, lapply(.SD, fun), by=sfac]
    setkey(otabN, sfac)
    OTN <- as.matrix(otabN, rownames=TRUE)
    ## now select the entries from colums that have the same values in
    ## the sample table...
    sdatN <- by(sdat, sdat[[by.sample]], function(x){
        sapply(x, function (y){
            uy <- unique(y)
            if(length(uy)==1) uy else paste(uy, collapse=";")
        })
    })
    sdatN <- as.data.frame(do.call(rbind, sdatN))
    phyloseq(otu_table(OTN, taxa_are_rows=FALSE),
             sample_data(sdatN),
             tax_table(PS))
}

PM <- sumTecRep(PS, by.sample="Animal.No")

OTT <- as.data.frame(otu_table(PM))
pdf("figures/primers_PMphyloseqTRmerged.pdf", width=45, height=15, onefile=FALSE)
pheatmap(log10(OTT+1))
dev.off()


### the individual amplicons to phyloseq
PS.l <- toPhyloseq(MA8, samples=colnames(MA8), multi2Single=FALSE)

PS.l <- lapply(PS.l, function (ps) prune_taxa(taxa_sums(ps) > 0, ps))

PM.l <- mclapply(PS.l, function(ps){
    p <- sumTecRep(ps, by.sample="Animal.No")
    prune_taxa(taxa_sums(p) > 0, p)
}, mc.cores=20)


############### Agglomerating by genus / OR NOT  ##########################
##### Raw non-merged technical replicates -> objects named "PSG": "PSG"
##### for merged amplicons "PSG.l" for the raw list per amplicosn
PSG <- tax_glom(PS, "genus")
PSG.l <- mclapply(PS.l, tax_glom, "genus", mc.cores=12)

PMGharsh <- sumTecRep(PSG, by.sample="Animal.No", 
                      fun=function(...) if (any(c(...) ==0)) {0} else {sum(...)})



##### Merged technical replicates -> objects named "PMG": "PMG" for
##### merged amplicons "PMG.l" for the raw list per amplicosn
## PGM <- tax_glom(PM, "genus")
## PGM.l <- mclapply(PM.l, tax_glom, "genus", mc.cores=12)


PGMharsh <- tax_glom(PMharsh, "genus")


## gimmeFoo <- function (Phy, tax){
##     sumSeqByTax <- function (Phy, tax) {
##         counts <- data.frame(cbind(asvCount=colSums(otu_table(Phy)), tax_table(Phy)))
##         counts$asvCount <- as.numeric(as.character(counts$asvCount))
##         tapply(counts$asvCount, counts[, tax], sum)
##     }

##     readNumByPhylum <- lapply(Phy, sumSeqByTax, tax)
##     names(readNumByPhylum) <- rownames(MA8)

##     mergeDf <- function(x, y) {
##         m <- merge(x, y,  by=0, all=TRUE)
##         rownames(m) <- m$Row.names
##         m$Row.names <- NULL
##         m
##     }

##     foo <- Reduce(mergeDf, readNumByPhylum)
##     colnames(foo) <- rownames(MA8)
##     foo[is.na(foo)] <- 0
##     foo
## }

## foo <- gimmeFoo(PS.l, "Phylum")

## foo <- gimmeFoo(PMharsh, "Phylum")

## ## Phyla detected by how many 16S primers?
## data.frame(
##     X16S=rowSums(foo[, grepl("16S", colnames(foo))]>0), 
##     ## detected by how many 28S primers
##     X28S=rowSums(foo[, grepl("28S", colnames(foo))]>0),
##     ## detected by how many 18S primers
##     X18S=rowSums(foo[, grepl("18S", colnames(foo))]>0)
## )

## bar <- gimmeFoo(PS.l, "genus")

## pdf("figures/how_many_per_primer.pdf")
## pheatmap(log10(bar+1))
## dev.off()

## pdf("figures/which_primer.pdf")
## pheatmap(apply(bar, 2, function (x) as.numeric(x>0)),
##          clustering_distance_rows="manhattan",
##          clustering_distance_cols="manhattan"
##          )
## dev.off()




## ## ## in how many amplicons are particular genera observed
## ## ## A function that could be included as quality check option in the PACKAGE 
## all.genera <- unique(unlist(lapply(readNumByGenus, names)))

## foo <- lapply(all.genera, function (x) {
##     lapply(readNumByGenus, function(y){
##         x%in%names(y)
##     })
## })

## num.prim <- unlist(lapply(foo, function(x) sum(unlist(x))))
## names(num.prim) <- all.genera

## tail(num.prim[order(num.prim)], n=300)

## ## ## in how many amplicons are genera observed at a high read proportion
## bar <- lapply(all.genera, function (x) {
##     lapply(readNumByGenus, function(y){
##         total <- sum(y)
##         prop <- y/total
##         x%in%(names(prop)[prop>0.00001])
##     })
## })

## num.H.prim <- unlist(lapply(bar, function(x) sum(unlist(x))))
## names(num.H.prim) <- all.genera

## tail(num.H.prim[order(num.H.prim)], n=100)


## ## ## in how many amplicons are genera observed at a high read numbers
## baz <- lapply(all.genera, function (x) {
##     lapply(readNumByGenus, function(y){
##         x%in%(names(y)[y>10])
##     })
## })

## num.H.prim <- unlist(lapply(baz, function(x) sum(unlist(x))))
## names(num.H.prim) <- all.genera

## tail(num.H.prim[order(num.H.prim)], n=100)




## This should acutally work with phyloseq's merge_samples function
## but doesn't as this messes up sample_data.
## CANDIDATE FOR INCLUSION IN PACKAGE...!

sumTecRep <- function (PS, by.sample, fun=sum){
    otab <- setDT(apply(otu_table(PS), 2, as.list))
    ## the columns giving numbers for sequences
    numcols <- colnames(otab)[nchar(colnames(otab))>10]
    sdat <- sample_data(PS, errorIfNULL = FALSE)
    otab[, (numcols):=lapply(.SD, as.numeric), .SDcols=numcols]
    otab[, sfac := as.factor(sdat[[by.sample]])]
    setkey(otab, sfac)
    otabN <- otab[, lapply(.SD, fun), by=sfac]
    setkey(otabN, sfac)
    OTN <- as.matrix(otabN, rownames=TRUE)
    ## now select the entries from colums that have the same values in
    ## the sample table...
    sdatN <- by(sdat, sdat[[by.sample]], function(x){
        sapply(x, function (y){
            uy <- unique(y)
            if(length(uy)==1) uy else paste(uy, collapse=";")
        })
    })
    sdatN <- as.data.frame(do.call(rbind, sdatN))
    phyloseq(otu_table(OTN, taxa_are_rows=FALSE),
             sample_data(sdatN),
             tax_table(PS))
}

PM <- sumTecRep(PS, by.sample="Animal.No")

OTT <- as.data.frame(otu_table(PM))
pdf("figures/primers_PMphyloseqTRmerged.pdf", width=45, height=15, onefile=FALSE)
pheatmap(log10(OTT+1))
dev.off()


### the individual amplicons to phyloseq
PS.l <- toPhyloseq(MA8, samples=colnames(MA8), multi2Single=FALSE)

PS.l <- lapply(PS.l, function (ps) prune_taxa(taxa_sums(ps) > 0, ps))

PM.l <- mclapply(PS.l, function(ps){
    p <- sumTecRep(ps, by.sample="Animal.No")
    prune_taxa(taxa_sums(p) > 0, p)
}, mc.cores=20)


############### Agglomerating by genus / OR NOT  ##########################
##### Raw non-merged technical replicates -> objects named "PSG": "PSG"
##### for merged amplicons "PSG.l" for the raw list per amplicosn
PSG <- tax_glom(PS, "genus")
PSG.l <- mclapply(PS.l, tax_glom, "genus", mc.cores=12)

PMGharsh <- sumTecRep(PSG, by.sample="Animal.No", 
                      fun=function(...) if (any(c(...) ==0)) {0} else {sum(...)})



##### Merged technical replicates -> objects named "PMG": "PMG" for
##### merged amplicons "PMG.l" for the raw list per amplicosn
## PGM <- tax_glom(PM, "genus")
## PGM.l <- mclapply(PM.l, tax_glom, "genus", mc.cores=12)


PGMharsh <- tax_glom(PMharsh, "genus")


## gimmeFoo <- function (Phy, tax){
##     sumSeqByTax <- function (Phy, tax) {
##         counts <- data.frame(cbind(asvCount=colSums(otu_table(Phy)), tax_table(Phy)))
##         counts$asvCount <- as.numeric(as.character(counts$asvCount))
##         tapply(counts$asvCount, counts[, tax], sum)
##     }

##     readNumByPhylum <- lapply(Phy, sumSeqByTax, tax)
##     names(readNumByPhylum) <- rownames(MA8)

##     mergeDf <- function(x, y) {
##         m <- merge(x, y,  by=0, all=TRUE)
##         rownames(m) <- m$Row.names
##         m$Row.names <- NULL
##         m
##     }

##     foo <- Reduce(mergeDf, readNumByPhylum)
##     colnames(foo) <- rownames(MA8)
##     foo[is.na(foo)] <- 0
##     foo
## }

## foo <- gimmeFoo(PS.l, "Phylum")

## foo <- gimmeFoo(PMharsh, "Phylum")

## ## Phyla detected by how many 16S primers?
## data.frame(
##     X16S=rowSums(foo[, grepl("16S", colnames(foo))]>0), 
##     ## detected by how many 28S primers
##     X28S=rowSums(foo[, grepl("28S", colnames(foo))]>0),
##     ## detected by how many 18S primers
##     X18S=rowSums(foo[, grepl("18S", colnames(foo))]>0)
## )

## bar <- gimmeFoo(PS.l, "genus")

## pdf("figures/how_many_per_primer.pdf")
## pheatmap(log10(bar+1))
## dev.off()

## pdf("figures/which_primer.pdf")
## pheatmap(apply(bar, 2, function (x) as.numeric(x>0)),
##          clustering_distance_rows="manhattan",
##          clustering_distance_cols="manhattan"
##          )
## dev.off()




## ## ## in how many amplicons are particular genera observed
## ## ## A function that could be included as quality check option in the PACKAGE 
## all.genera <- unique(unlist(lapply(readNumByGenus, names)))

## foo <- lapply(all.genera, function (x) {
##     lapply(readNumByGenus, function(y){
##         x%in%names(y)
##     })
## })

## num.prim <- unlist(lapply(foo, function(x) sum(unlist(x))))
## names(num.prim) <- all.genera

## tail(num.prim[order(num.prim)], n=300)

## ## ## in how many amplicons are genera observed at a high read proportion
## bar <- lapply(all.genera, function (x) {
##     lapply(readNumByGenus, function(y){
##         total <- sum(y)
##         prop <- y/total
##         x%in%(names(prop)[prop>0.00001])
##     })
## })

## num.H.prim <- unlist(lapply(bar, function(x) sum(unlist(x))))
## names(num.H.prim) <- all.genera

## tail(num.H.prim[order(num.H.prim)], n=100)


## ## ## in how many amplicons are genera observed at a high read numbers
## baz <- lapply(all.genera, function (x) {
##     lapply(readNumByGenus, function(y){
##         x%in%(names(y)[y>10])
##     })
## })

## num.H.prim <- unlist(lapply(baz, function(x) sum(unlist(x))))
## names(num.H.prim) <- all.genera

## tail(num.H.prim[order(num.H.prim)], n=100)

