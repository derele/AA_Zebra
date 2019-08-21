## weird idea...

    blt <- blt[,.(bitdiff= max(bitsum) - bitsum, bitsum,
                  superkingdom, phylum, class, order, family, genus, species),
               by="ampProd"]
    get.opt.gen <- function(bdiff, what, Nsupport=1, eval=TRUE){
        get.unique.or.na <- function (x){
            ## unique taxa at that level excluding potential NA's 
            agnostic <- as.character(x)
            taxa <- agnostic[!is.na(agnostic)]
            ux <- unique(taxa)
            ## but return NA if they are not unique
            if(length(taxa)>=Nsupport && ## number of supporting annotations
               length(ux)==1){ ## has to be a unique answer
                return(ux)
            } else {as.character(NA)}
        }
        df <- blt[bitdiff <= bdiff, .(get.unique.or.na(eval(as.name((what))))),
                     by=ampProd]  
        setnames(df, "V1", what)
        if(!eval){ ## fro final result return the data
            df
        } else{ ## for optimization return the number of NA results to
                ## minimize
            nrow(df) - 
            nrow(df[!is.na(df[, eval(as.name(what))])])
        }
    }
    taxa <- c("superkingdom", "phylum", "order", "family", "genus")
    annot.l <- lapply(taxa, function (x) {
        ## Optim <- optimize(get.opt.gen, upper=10^5, lower=0,
        ##                   what=x)
        ## tax <- get.opt.gen(Optim$minimum, what=x, eval=FALSE)
        ## list(tax, Optim$minimum)
        tax <- get.opt.gen(2, what=x, eval=FALSE)
        list(tax, 2)

    })
    annot <- Reduce("merge", lapply(annot.l, "[[", 1))
    message("optimized at bitscore differeence:\n",
            lapply(annot.l, function (x) {
                paste(round(x[[2]], 2), "\n")}))
    ## now we have to break this up into an annotation list for each
    ## amplicon
    annot.l <- by(annot, gsub("_S_\\d+", "", annot$ampProd),
                     function (x) return(x))
    annot.i <- as.numeric(gsub("A_(\\d+)(_R_)?", "\\1", names(annot.l)))
    names(annot.l) <- rownames(MA)[annot.i]

