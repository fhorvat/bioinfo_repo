
# --------------------------------------------------------------------------------------------------------- #
# Reads the gtf annotation
ReadGTFAnnotation = function(gtf.path, which.regions='exon', ensembl=FALSE, select=TRUE, annot=TRUE, annot.level='gene', type='file'){
    
    lib.path=file.path(Sys.getenv('MYLIB'),'RFun')
    source(file.path(lib.path, 'ScanLib.R'))
    source(file.path(lib.path, 'GeneFunctions.R'))
    require(data.table)
    require(genomation)
    require(GenomicRanges)
    require(stringr)
    
       
    if(!file.exists(gtf.path))
        stop('the gtf file does not exist')
   
    if(type=='RData'){
        Assigner(gtf.path, 'gtf')
        
    }else if(type == 'file'){
        message('Reading gtf file...')
        gtf = gffToGRanges(gtf.path, ensembl=ensembl)
        seqlevels(gtf, force=TRUE) = seqlevels(gtf)[!str_detect(seqlevels(gtf),'NT')]
        seqlevels(gtf)[seqlevels(gtf) == 'chrMT'] = 'chrM'
    }
    
    if(which.regions != 'all')
        gtf = gtf[gtf$type %in% which.regions]
    
    gtf.selected=NULL
    if(select){
        message('Selecting transcripts...')
        gtf.selected = GtfSelectTranscript(gtf[gtf$type=='exon'])
        gtf.selected = split(gtf.selected, gtf.selected$transcript_id)
    }
    
    gtf.annot = NULL
    if(annot){
        message('Constructing annotation...')
        
        gtf.exon = gtf[gtf$type == 'exon']
        gtf.annot = unique(as.data.frame(DataFrame(gtf.exon))[,c('gene_id','transcript_id','gene_name','gene_biotype')])
       
            gtf.exon.ge = split(gtf.exon, gtf.exon$gene_id)
            gtf.range.ge = unlist(range(gtf.exon.ge))
            gtf.annot$gcoord  = GCoords(gtf.range.ge[match(gtf.annot$gene_id, names(gtf.range.ge))])
            gtf.annot$gwidth = width(gtf.range.ge[match(gtf.annot$gene_id, names(gtf.range.ge))])
        
        if(annot.level == 'transcript'){
            
            gtf.exon.tr = split(gtf.exon, gtf.exon$transcript_id)
            gtf.range.tr = unlist(range(gtf.exon.tr))
            gtf.annot$tcoord  = GCoords(gtf.range.tr[match(gtf.annot$transcript_id, names(gtf.range.tr))])
            gtf.annot$twidth  = sum(width(gtf.exon.tr))[gtf.annot$transcript_id]
        }
        
    }
    
    return(list(gtf=gtf, gtf.sel = gtf.selected, annot=gtf.annot))
}

        
    



