

# A treemap R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","freqInDbPercent","abslog10pvalue","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0003674","molecular_function",100.000,3.1337,1.000,0.000,"molecular_function"),
c("GO:0003824","catalytic activity",27.857,11.1925,0.981,0.000,"catalytic activity"),
c("GO:0005488","binding",64.149,23.9788,0.991,0.000,"binding"),
c("GO:0008047","enzyme activator activity",2.021,4.2510,0.912,0.000,"enzyme activator activity"),
c("GO:0019869","chloride channel inhibitor activity",0.043,2.8579,0.931,0.454,"enzyme activator activity"),
c("GO:0001671","ATPase activator activity",0.091,3.0325,0.918,0.544,"enzyme activator activity"),
c("GO:0015197","peptide transporter activity",0.043,2.9928,0.953,0.000,"peptide transporter activity"),
c("GO:0008565","protein transporter activity",0.417,2.9899,0.951,0.426,"peptide transporter activity"),
c("GO:0016462","pyrophosphatase activity",3.956,9.5528,0.799,0.000,"pyrophosphatase activity"),
c("GO:0003678","DNA helicase activity",0.211,4.6021,0.842,0.660,"pyrophosphatase activity"),
c("GO:0016817","hydrolase activity, acting on acid anhydrides",3.975,9.4425,0.852,0.577,"pyrophosphatase activity"),
c("GO:0016787","hydrolase activity",12.022,5.2596,0.892,0.345,"pyrophosphatase activity"),
c("GO:0016779","nucleotidyltransferase activity",0.584,3.4908,0.916,0.220,"pyrophosphatase activity"),
c("GO:0044769","ATPase activity, coupled to transmembrane movement of ions, rotational mechanism",0.148,3.4921,0.825,0.636,"pyrophosphatase activity"),
c("GO:0016740","transferase activity",11.141,3.7423,0.893,0.428,"pyrophosphatase activity"),
c("GO:0051020","GTPase binding",1.451,8.2676,0.814,0.000,"GTPase binding"),
c("GO:0051082","unfolded protein binding",0.402,2.9853,0.885,0.118,"GTPase binding"),
c("GO:0045296","cadherin binding",1.399,5.7305,0.871,0.135,"GTPase binding"),
c("GO:0032403","protein complex binding",3.611,6.9318,0.842,0.152,"GTPase binding"),
c("GO:0031625","ubiquitin protein ligase binding",1.427,3.2557,0.815,0.525,"GTPase binding"),
c("GO:0008022","protein C-terminus binding",0.872,4.1938,0.877,0.128,"GTPase binding"),
c("GO:0019900","kinase binding",2.931,3.8928,0.802,0.573,"GTPase binding"),
c("GO:0019899","enzyme binding",8.732,19.5986,0.844,0.203,"GTPase binding"),
c("GO:0019902","phosphatase binding",0.747,3.7799,0.825,0.488,"GTPase binding"),
c("GO:0005484","SNAP receptor activity",0.177,2.9185,0.892,0.108,"GTPase binding"),
c("GO:0031072","heat shock protein binding",0.426,3.4473,0.884,0.118,"GTPase binding"),
c("GO:0044389","ubiquitin-like protein ligase binding",1.451,3.8069,0.814,0.526,"GTPase binding"),
c("GO:0050839","cell adhesion molecule binding",2.040,3.1733,0.867,0.159,"GTPase binding"),
c("GO:0031369","translation initiation factor binding",0.144,3.7100,0.893,0.106,"GTPase binding"),
c("GO:0046982","protein heterodimerization activity",2.309,3.0386,0.865,0.162,"GTPase binding"),
c("GO:0042802","identical protein binding",6.734,3.6021,0.849,0.223,"GTPase binding"),
c("GO:0031267","small GTPase binding",1.341,6.0625,0.805,0.521,"GTPase binding"),
c("GO:0008092","cytoskeletal protein binding",4.014,4.5287,0.857,0.175,"GTPase binding"),
c("GO:0008135","translation factor activity, RNA binding",0.441,7.5302,0.800,0.033,"translation factor activity, RNA binding"),
c("GO:0003676","nucleic acid binding",18.627,13.4168,0.718,0.430,"translation factor activity, RNA binding"),
c("GO:0042162","telomeric DNA binding",0.158,3.0376,0.808,0.412,"translation factor activity, RNA binding"),
c("GO:1901265","nucleoside phosphate binding",10.159,9.7799,0.742,0.461,"translation factor activity, RNA binding"),
c("GO:0070181","small ribosomal subunit rRNA binding",0.038,3.1062,0.834,0.359,"translation factor activity, RNA binding"),
c("GO:0032555","purine ribonucleotide binding",9.009,8.8297,0.602,0.683,"translation factor activity, RNA binding"),
c("GO:0003743","translation initiation factor activity",0.297,6.0768,0.806,0.424,"translation factor activity, RNA binding"),
c("GO:0019001","guanyl nucleotide binding",1.854,4.2097,0.677,0.544,"translation factor activity, RNA binding"),
c("GO:0005525","GTP binding",1.763,4.5317,0.658,0.193,"translation factor activity, RNA binding"),
c("GO:0043168","anion binding",12.803,6.8962,0.870,0.296,"translation factor activity, RNA binding"),
c("GO:0003723","RNA binding",7.582,16.9066,0.737,0.342,"translation factor activity, RNA binding"),
c("GO:0043047","single-stranded telomeric DNA binding",0.048,4.4698,0.803,0.209,"translation factor activity, RNA binding"),
c("GO:0043531","ADP binding",0.144,3.5186,0.738,0.418,"translation factor activity, RNA binding"),
c("GO:0005515","protein binding",39.060,20.9031,0.880,0.062,"protein binding"),
c("GO:0097159","organic cyclic compound binding",28.470,17.6108,0.885,0.145,"protein binding"),
c("GO:1901363","heterocyclic compound binding",28.039,18.6126,0.885,0.144,"protein binding"),
c("GO:0097367","carbohydrate derivative binding",10.676,6.8477,0.898,0.100,"protein binding"),
c("GO:0043167","ion binding",28.389,6.5952,0.885,0.145,"protein binding"),
c("GO:0036094","small molecule binding",11.491,8.7721,0.897,0.103,"protein binding"),
c("GO:0044877","macromolecular complex binding",6.643,8.4802,0.903,0.087,"macromolecular complex binding"),
c("GO:0051787","misfolded protein binding",0.057,3.1385,0.900,0.098,"misfolded protein binding"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
tmPlot(
	stuff,
	index = c("representative","description"),
	vSize = "abslog10pvalue",
	type = "categorical",
	vColor = "representative",
	title = "REVIGO Gene Ontology treemap",
	inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
	lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
	bg.labels = "#CCCCCCAA",     # define background color of group labels
												       # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
	position.legend = "none"
)

dev.off()
