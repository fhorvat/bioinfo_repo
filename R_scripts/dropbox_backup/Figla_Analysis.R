load('GSE5558_Figla_KO.RData')
library(limma)
library(dplyr)
library(stringr)


# ------------------------------------------------------------------------------ #
# functions for making contrasts
makeBinaryContrasts = function(data,column='sample'){
  
  if(class(data) == 'data.frame')
    contrasts = expand.grid(sort(unique(data[[column]])), sort(unique(data[[column]])))
  
  if(class(data) == 'character')
    contrasts = expand.grid(sort(unique(data)), sort(unique(data)))
  
  if(class(data) == 'factor')
    contrasts = expand.grid(sort(levels(data)), sort(levels(data)))
  
  contrasts = subset(contrasts, Var1 != Var2)
  contrasts = contrasts[order(contrasts$Var1),]
  contrasts = contrasts[!duplicated(apply(contrasts, 1, function(x)paste(sort(x), collapse='-'))),]
  contrasts = contrasts[order(contrasts[,1]),]
  contrasts = with(contrasts, paste(Var1, Var2, sep='-'))
  return(contrasts)
}

# ------------------------------------------------------------------------------ #
makeContlist = function(
  Factor,
  column  = 'Factor',
  sort_by = 'up'
){
  
  suppressPackageStartupMessages(library(dplyr))
  if(class(Factor) == 'data.frame')
    Factor = Factor[[column]]
  
  Factor = as.character(unique(Factor))
  contrasts = expand.grid(Factor, Factor, stringsAsFactors=FALSE)
  contrasts = with(contrasts,
                   data.frame(
                     Var1 = pmin(Var1, Var2),
                     Var2 = pmax(Var1, Var2)
                   )) %>%
    distinct() %>%
    dplyr::filter(Var1 != Var2)
  
  if(sort_by == 'up'){
    contrasts =  contrasts %>%
      dplyr::rename(X1 = Var1, X2 = Var2)
  }
  if(sort_by == 'down'){
    contrasts =  contrasts %>%
      dplyr::rename(X1 = Var2, X2 = Var1)
  }
  contnames =  with(contrasts, (results = paste(X1, X2, sep='_')))
  
  contlist = split(contrasts, contnames)
  contlist = lapply(contlist,unlist)
  return(contlist)
}
# ---------------------------------------------------------------------------- #
get_limma_tab=function(expr, samps, lfc=1, padj=0.05, method='ls', covar=NULL){
  
  library(limma)
  lm = get_limma(expr, samps, method=method, covar=covar)
  cont = makeBinaryContrasts(unique(samps))
  res = getResults_limma(lm, cont, lfc=lfc, pval=padj)
  
  if(class(expr) == 'expressionSet'){
    
    dat = as(featureData(expr),'data.frame') %>%
      dplyr::select(1,9,10,11)
    colnames(dat) = str_replace(colnames(dat),' ','_')
    colnames(dat) = tolower(colnames(dat))
    means = getMeans(exprs(expr), samps, unique=FALSE)
    means$id = rownames(means)
    tab = merge(dat, res, by='id')
  }else{
    means = getMeans(expr, samps, unique=FALSE)
    means$id = rownames(means)
    tab = res
  }
  
  tab = merge(tab, means, by='id')
  return(tab)
}


# ---------------------------------------------------------------------------- #
get_limma = function(eset, samps, method='ls', covar=NULL){
  
  message('Contrasts... ')
  cont=makeBinaryContrasts(samps)
  contrast.matrix = makeContrasts(contrasts=cont,
                                  levels=unique(samps))
  
  message('Design... ')
  design = model.matrix(~0+samps)
  colnames(design) = str_replace(colnames(design),'samps','')
  design = design[,match(rownames(contrast.matrix), colnames(design))]
  if(!is.null(covar)){
    design = cbind(design,covar)
    contrast.matrix = rbind(contrast.matrix,matrix(0, nrow=ncol(covar), ncol=ncol(contrast.matrix)))
  }
  
  
  message('Fit... ')
  fit  = lmFit(eset, design, method=method)
  fit2 = contrasts.fit(fit, contrast.matrix)
  message('eBayes... ')
  fit2 = eBayes(fit2, robust=TRUE)
  return(fit2)
}


e = g[[1]]@assayData$exprs
d = g[[1]]@phenoData@data %>%
  dplyr::select(title) %>% 
  mutate(sample_names = rownames(.)) %>% 
  mutate(title = sub('ovary dyeswap','WT',title)) %>% 
  mutate(title = sub('ovary','KO',title)) %>% 
  mutate(title = sub('replicate ','rep',title)) %>% 
  mutate(title = gsub(' ','_',title)) %>% 
  mutate(sample = sub('_rep.','',title)) %>% 
  magrittr::set_rownames(.$sample_names)

pheno = AnnotatedDataFrame(d)
eset = ExpressionSet(assayData = e, phenoData = pheno)
l = get_limma(eset, pheno$sample)

  
