# ------------------------------------------------------------------- #
# a function that gets the first n principal components for plotting
getPCA = function (x, intgroup = "condition", ntop = 500, pcs=3, cdata=NULL) 
    {     
        library(genefilter)
        library(GenomicRanges)
        
        if(class(x) == 'SeqExpressionSet'){
            cdata = pData(x)
            x = counts(x)
            
        }
            
        if(class(x) %in% c('SummarizedExperiment','DESeqDataSet')){
            cdata = colData(x)
            x = assays(x)[[1]]
            
        }
        if(class(x) == 'matrix'){
            cdata = cdata
            x = x
            
        }
            
        
        rv <- rowVars(x)
        select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                           length(rv)))]
        library(pcaMethods)
        pca <- prcomp(t(x[select, ]))
        percentVar <- pca$sdev^2/sum(pca$sdev^2)
        if (!all(intgroup %in% names(cdata))){
            stop("the argument 'intgroup' should specify columns of cdata")
        }         
        intgroup.df <- as.data.frame(cdata[, intgroup, drop = FALSE])
        group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
        pc = data.frame(pca$x[,1:min(pcs, ncol(pca$x))])
        d <- data.frame(pc, group = group, 
                        intgroup.df, names = colnames(x))
        attr(d, "percentVar") <- percentVar[1:pcs]
        rownames(d) = colnames(x)
        return(d) 
    }


# ------------------------------------------------------------------- #
# a function that gets the first n principal components for plotting



# ------------------------------------------------------------------- #
# a function that gets the first n principal components for plotting