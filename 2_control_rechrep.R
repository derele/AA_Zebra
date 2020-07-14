library(tidyverse)



plotcor <- function(ps, plot=TRUE){
    PSdup <- subset_samples(ps,
                            Animal.No%in%Animal.No[duplicated(Animal.No)])
    PSdupDat <- t(otu_table(PSdup))

    foo <- sapply(unique(sample_data(PSdup)$Animal.No),
                  function (x) PSdupDat[, sample_data(PSdup)$Animal.No%in%x])

    foo <- cbind(foo, repls=rep(c("rep1", "rep2"), each=nrow(foo)/2),
                 ASV=rep(1:(nrow(foo)/2), times=2))

    foo[,which(!colnames(foo)%in%c("repls", "ASV"))] <-
        sapply(foo[,which(!colnames(foo)%in%c("repls", "ASV"))],
               function (x) as.numeric(as.character(x)))
    
    bar <- gather(as.data.frame(foo), animal, count,
                  colnames(foo)[!colnames(foo)%in%c("repls", "ASV")])

    bar$count <- as.numeric(as.character(bar$count))

    baz <- spread(bar, repls,  count)

    corplot <- ggplot(baz, aes(rep1+1, rep2+1)) +
        geom_jitter(alpha=0.1, height=0.1, width=0.1) +
        scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                      labels = scales::trans_format("log10", scales::math_format(10^.x))) +
        scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                      labels = scales::trans_format("log10", scales::math_format(10^.x))) +
        annotation_logticks() +
        theme_bw()
    
    if (isTRUE(plot)) {
        corplot
    } else if (isFALSE(plot)) {
        baz
    } else {
        list(corplot, baz)
    }
}

corAll <- plotcor(PS, plot="all")

corHarsh <- plotcor(PMharsh, plot="all")

png("figures/tech_rep_overall.png")
plot(corAll[[1]])
dev.off()

png("figures/tech_rep_genus.png")
plotcor(tax_glom(PS, "genus"))
dev.off()

png("figures/tech_rep_family.png")
plotcor(tax_glom(PS, "family"))
dev.off()


