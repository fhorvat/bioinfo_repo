#= @pepap 20200714

library(parallel)
library(DECIPHER)
library(ape)
library(phangorn)

cat("\n Created functions : \"pep.alnSeqs\", \"pep.mutRate\", \"pep.phylo\", \"pep.seq2subsRate\"\n\n",sep="")

#= (1) DNAStringSet => multi-alignment
pep.alnSeqs <- function( inpDSS,ITER=100,REFIN=100,PROC=detectCores(),VERBOSE=T ) {

#> Align sequences
 xALN <- DECIPHER::AlignSeqs( myXStringSet=inpDSS,iterations=ITER,refinements=REFIN,processors=PROC,verbose=VERBOSE )

 return( xALN )

}

#= (2) Multi-alignment => substitution matrix
pep.mutRate <- function( inpDALN,INDELS=F,NDIGITS=5 ) {

#> Create output matrix
 out.mat  <- matrix( nrow=length(inpDALN),ncol=length(inpDALN),dimnames=list(names(inpDALN),names(inpDALN)) )
#> Convert alignments to character strings
 list_CHR <- strsplit( x=as.character(inpDALN),split="",fixed=T )
#> Fill matrix by pair-sequence comparison
 for (  i in seq(from=1,to=length(inpDALN)) ) {
  for ( j in seq(from=1,to=length(inpDALN)) ) {
#> Filter out alignment-positions which equal "-" in both sequences
   not_dash  <- ( (list_CHR[[i]]!="-") | (list_CHR[[j]]!="-") )
#> Filter out matches
   N_matches <- sum(   list_CHR[[i]][not_dash]==list_CHR[[j]][not_dash] )
   if ( INDELS ) {
#> Filter out mismatches + indels
    N_substit <- sum( ( list_CHR[[i]][not_dash]!=list_CHR[[j]][not_dash] )                                                                       )
   } else        {
#> Filter out mismatches & remove indels
    N_substit <- sum( ( list_CHR[[i]][not_dash]!=list_CHR[[j]][not_dash] ) & ( list_CHR[[i]][not_dash]!="-" ) & ( list_CHR[[j]][not_dash]!="-" ) )
   }
#> Calculate substitution rate
   x_val        <- round( x=(N_substit/(N_matches+N_substit)),digits=NDIGITS )
   out.mat[i,j] <- x_val
   out.mat[j,i] <- x_val
  }
 }

 return( out.mat )

}

#= (3) Substitution matrix => shortest edge-lengths
pep.phylo <- function( inpMUTMAT,RETURN.PHYLO=F ) {

#> Transform the substitution-matrix into phylogenic-tree object
 out.upgma       <- phangorn::upgma( D=inpMUTMAT )
 if ( RETURN.PHYLO ) {
  return( out.upgma )
 }
#> Extract shortest edge-lengths
 out.dist        <- out.upgma$edge.length[                 out.upgma$edge[ , 2 ] <= nrow(inpMUTMAT)       ]
 names(out.dist) <- out.upgma$tip.label[   out.upgma$edge[ out.upgma$edge[ , 2 ] <= nrow(inpMUTMAT) , 2 ] ]

 return( out.dist )

}

#= (4) Function performing all previous transformations step-by-step
#> n.iter       ... number of iterations within the multi-alignment procedure                             ( default 100 )
#> n.refin      ... number of refinement steps within the multi-alignment procedure                       ( default 100 )
#> ncpus        ... number of processors used within the multi-alignment procedure                        ( default maximum detectable processors )
#> quiet        ... turn off verbose output from the multi-alignment procedure                            ( default FALSE : if TRUE, the alignment convergence is not reported )
#> add.indels   ... include indels into the substitution rate                                             ( default FALSE : if TRUE, try to refine your input sequences & remove all additional inserts )
#> n.digits     ... number of digits to which the value of substitution rate is rounded                   ( default 5 )
#> return.phylo ... instead of set of shortest-length values return the original phylogenetic-tree object ( default FALSE : if TRUE, you can directly plot the object )
pep.seq2subsRate <- function( inpDNAStringSet,n.iter=100,n.refin=100,ncpus=detectCores(),quiet=F,add.indels=F,n.digits=5,return.phylo=F ) {

 cat("\n ++ Input DNAStringSet => multi-alignment ++\n",sep="")
 out.aln <- pep.alnSeqs( inpDSS=inpDNAStringSet,ITER=n.iter,REFIN=n.refin,PROC=ncpus,VERBOSE=(!quiet) )
# return( out.aln )
 cat(" ++ Multi-alignment => substitution matrix ++\n",sep="")
 out.mut <- pep.mutRate( inpDALN=out.aln,INDELS=add.indels,NDIGITS=n.digits )
# return( out.mut )
 cat(" ++ Substitution matrix => shortest edge-lengths ++\n",sep="")
 out.dis <- pep.phylo( inpMUTMAT=out.mut,RETURN.PHYLO=return.phylo )
 cat(" ++ Finished ++\n\n",sep="")
 return( out.dis )

}

