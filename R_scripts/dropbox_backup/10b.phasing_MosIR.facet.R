### INFO: 
### DATE: Tue Sep 11 09:14:41 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/phasing")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)
library(Biostrings)
library(ggiraphExtra)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS
# get coverage from bam file
phasingRegister <- function(bam_path, which_gr, read_length){
  
  # read bam
  bam_gr <- 
    GenomicAlignments::readGAlignmentsList(file = bam_path, 
                                           use.names = TRUE, 
                                           param = ScanBamParam(which = which_gr, 
                                                                flag = scanBamFlag(isMinusStrand = F))) %>% 
    unlist(.)
  
  # take only reads between 21-23nt
  bam_gr <- bam_gr[str_detect(cigar(bam_gr), str_c(read_length, "M"))]
  
  # get starting positions of reads, normalize to start of mosIR arm, get mod 21
  if(values(which_gr)$name == "MosIR_arm1"){
    
    # if read is on MosIR arm 1 just set start of arm 1 to 1
    values(bam_gr)$start_norm <- start(bam_gr) - start(which_gr) + 1
    values(bam_gr)$end_norm <- end(bam_gr) - start(which_gr) + 1
    
  }else{
    if(values(which_gr)$name == "MosIR_arm2"){
      
      # if read is on MosIR arm 2 then start is really an end of mosIR arm 2
      values(bam_gr)$start_norm <- end(which_gr) - end(bam_gr) + 1
      values(bam_gr)$end_norm <- end(which_gr) - start(bam_gr) + 1
      
    }else{
      
      # if non true then input is wrong
      stop("Something went wrong!")
      
    }
    
  }
  
  # get mod 21 of normalized read starts
  read_start_mod <- (values(bam_gr)$start_norm %% read_length) + 1
  
  # table for plot
  read_register_df <- 
    tibble(read_register = read_start_mod) %>% 
    dplyr::group_by(read_register) %>% 
    dplyr::summarize(count = n()) %>% 
    dplyr::mutate(percent = (count / sum(count)) * 100) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(-count) %>% 
    dplyr::right_join(tibble(read_register = 1:read_length), by = "read_register") %>% 
    dplyr::mutate(percent = replace(percent, is.na(percent), 0)) %>% 
    tidyr::spread(read_register, percent)
  
  # return 
  return(read_register_df)
  
}

ggRadar2 <- function(data, mapping = NULL, rescale = TRUE, legend.position = "top",
                     colour = "red", alpha = 0.3, size = 3, ylim = NULL, scales = "fixed",
                     use.label = FALSE, interactive = FALSE, ...){
  
  data = as.data.frame(data)
  (groupname = setdiff(names(mapping), c("x", "y")))
  groupname
  mapping
  length(groupname)
  if (length(groupname) == 0) {
    groupvar <- NULL
  }
  else {
    groupvar = getMapping(mapping, groupname)
  }
  groupvar
  facetname <- colorname <- NULL
  if ("facet" %in% names(mapping)) {
    facetname <- getMapping(mapping, "facet")
  }
  (colorname = setdiff(groupvar, facetname))
  if ((length(colorname) == 0) & !is.null(facetname))
    colorname <- facetname
  data = num2factorDf(data, groupvar)
  (select = sapply(data, is.numeric))
  if ("x" %in% names(mapping)) {
    xvars = getMapping(mapping, "x")
    xvars
    if (length(xvars) < 3)
      warning("At least three variables are required")
  }
  else {
    xvars = colnames(data)[select]
  }
  (xvars = setdiff(xvars, groupvar))
  if (rescale)
    data = rescale_df(data, groupvar)
  temp = sjlabelled::get_label(data)
  cols = ifelse(temp == "", colnames(data), temp)
  if (is.null(groupvar)) {
    id = newColName(data)
    data[[id]] = 1
    longdf = reshape2::melt(data, id.vars = id, measure.vars = xvars)
  }
  else {
    cols = setdiff(cols, groupvar)
    longdf = reshape2::melt(data, id.vars = groupvar, measure.vars = xvars)
  }
  temp = paste0("plyr::ddply(longdf,c(groupvar,'variable'),summarize,mean=mean(value,na.rm=TRUE))")
  df = eval(parse(text = temp))
  colnames(df)[length(df)] = "value"
  df
  groupvar
  if (is.null(groupvar)) {
    id2 = newColName(df)
    df[[id2]] = "all"
    id3 = newColName(df)
    df[[id3]] = 1:nrow(df)
    df$tooltip = paste0(df$variable, "=", round(df$value,
                                                1))
    df$tooltip2 = paste0("all")
    p <- ggplot(data = df, aes_string(x = "variable", y = "value",
                                      group = 1)) + geom_polygon_interactive(aes_string(tooltip = "tooltip2"),
                                                                             colour = colour, fill = colour, alpha = alpha) +
      geom_point_interactive(aes_string(data_id = id3,
                                        tooltip = "tooltip"), colour = colour, size = size)
  }
  else {
    if (!is.null(colorname)) {
      id2 = newColName(df)
      df[[id2]] = df[[colorname]]
    }
    id3 = newColName(df)
    df[[id3]] = 1:nrow(df)
    df$tooltip = paste0(groupvar, "=", df[[colorname]], "<br>",
                        df$variable, "=", round(df$value, 1))
    df$tooltip2 = paste0(groupvar, "=", df[[colorname]])
    p <- ggplot(data = df, aes_string(x = "variable", y = "value",
                                      colour = colorname, fill = colorname, group = colorname)) +
      geom_polygon_interactive(aes_string(tooltip = "tooltip2"),
                               alpha = alpha) + geom_point_interactive(aes_string(data_id = id3,
                                                                                  tooltip = "tooltip"), size = size)
  }
  p
  if (!is.null(facetname)) {
    formula1 = as.formula(paste0("~", facetname))
    p <- p + facet_wrap(formula1, scales = scales)
  }
  p <- p + xlab("") + ylab("") + theme(legend.position = legend.position)
  if (use.label)
    p <- p + scale_x_discrete(labels = cols)
  if (!is.null(ylim))
    p <- p + scale_y_continuous(limits = ylim, expand = c(0, 0))
  p <- p + coord_radar()
  # if (!is.null(ylim))
  #   p <- p + expand_limits(y = ylim)
  p
  if (interactive) {
    tooltip_css <- "background-color:white;font-style:italic;padding:10px;border-radius:10px 20px 10px 20px;"
    hover_css = "r:4px;cursor:pointer;stroke-width:6px;"
    selected_css = "fill:#FF3333;stroke:black;"
    p <- ggiraph(code = print(p), tooltip_extra_css = tooltip_css,
                 tooltip_opacity = 0.75, zoom_max = 10, hover_css = hover_css,
                 selected_css = selected_css)
  }
  p
  
}

roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
  if(length(x) != 1) stop("'x' must be of length 1")
  10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# MosIR .bed path
mosir_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/Documentation/pCAG-EGFP_MosIR.coordinates.bed"

######################################################## READ DATA
# read bed with coordinates
mosir_gr <- rtracklayer::import.bed(con = mosir_path)
mosir_gr <- mosir_gr[mosir_gr$name != "EGFP"]
start(mosir_gr) <- start(mosir_gr) - 1

######################################################## MAIN CODE
# ESC bams
esc_bams <- c("s_ES_Dcr_Rsc_tran_r2.SE.mis_0", "s_ES_Dcr_Rsc_utran_r2.SE.mis_0",
              "s_ES_DcrSom_Rsc_tran_r1.SE.mis_0", "s_ES_DcrSom_Rsc_utran_r1.SE.mis_0", 
              "s_ES_WT_tran_r1.SE.mis_0")

# experiment paths
experiment_paths <- 
  c("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/SOLiD/T3T_DcrTrans_2011/Data/Mapped/Shrimp_mm10", 
    "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/SOLiD/ES_DcrTrans_2012/Data/Mapped/Shrimp_mm10") %>% 
  magrittr::set_names(., str_remove(., "Data.*$") %>% basename(.))

# list of experiments
experiment_list <- names(experiment_paths)

# loop through experiments
for(experiment in experiment_list){
  
  # print experiment
  cat("Plotting experiment:", experiment, "\n")
  
  # get mapped path 
  mapped_path <- experiment_paths[experiment]
  
  # list bams 
  bam_paths <- list.files(path = mapped_path, pattern = ".*mis_0.*bam$", full.names = T)
  
  # filter bams in ESC
  # if(experiment == "ES_DcrTrans_2012"){
  #   bam_paths <- bam_paths[str_detect(bam_paths, str_c(esc_bams, collapse = "|"))]
  # }
  
  # save radar plots of 21-23 read lengths
  for(read_length in 21:23){
    
    # loop through sample in experiment 
    plot_facet_df <- purrr::map(.x = bam_paths, function(bam_path){
      
      # get bam name
      bam_name <- basename(bam_path) %>% str_remove(., ".bam")
      
      # print bam name
      cat(str_c(bam_name, read_length, sep = "."), "\n")
      
      # get coverage of both mosIR arms in one data.frame, normalize for library size
      plot_df <-
        purrr::map(.x = 1:2, function(x){
          
          # get coverage
          phasing_df <-
            phasingRegister(bam_path = bam_path, which_gr = mosir_gr[x], read_length = read_length) %>%
            dplyr::mutate(arm = mcols(mosir_gr[x])$name)
          
          return(phasing_df)
          
        }) %>%
        dplyr::bind_rows(.) %>% 
        dplyr::mutate(sample_id = bam_name) %>% 
        dplyr::select(sample_id, arm, everything())
      
    }) %>% 
      dplyr::bind_rows(.)
    
    # create and save radar plot
    ggRadar(data = plot_facet_df, aes(color = arm, facet = sample_id), 
            rescale = F, interactive = F, size = 0.1, 
            ylim = c(0, roundUpNice(max(plot_facet_df[, -c(1, 2)])))) +
      facet_wrap(~ sample_id) +
      theme_bw() +
      theme(legend.position = "bottom", 
            axis.text.x = element_text(size = 10),
            axis.text.y = element_text(size = 10), 
            panel.grid.major = element_line(size = 1, color = "gray90")) +
      ggsave(filename = file.path(outpath, str_c(experiment, "MosIR_phase", read_length, "nt.facet.png", sep = ".")), width = 20, height = 20)
    
  }
  
  cat("\n")
  
}



