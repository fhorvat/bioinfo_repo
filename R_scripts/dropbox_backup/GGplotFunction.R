ggplotScatter = function(d, x, y, diff=FALSE, plot.cor=FALSE){
	g=ggplot(d, aes(x=x,y=y)) + 
	 		 geom_point(color='black', cex=2)
	if(diff){
		if(!diff %in% colnames(d))
			stop(paste(diff, 'is not a designated colnames variable'))
		if(!length(unique(diff) == 3)
		g = geom_point(data=d[d$MII.GV>quantile(d$MII.GV,.9),], color='blue',  cex=2) +
			geom_point(data=d[d$MII.GV<quantile(d$MII.GV,.1),], color='red', cex=2) 
	}		 
	
	g = g + xlab(x) + ylab(y) 
	g = g + annotate("text", x =min(5), y = max(15), label = paste('cor:', cort), size=10)
	g = g + labs(title=name)+
	        theme(plot.title = element_text(size = rel(2)), 
		 		  axis.text  = element_text(colour = "black", size=rel(2)),
				  axis.title = element_text(size = rel(2)),
				  axis.line  = element_line(size = rel(3)),
				  panel.background = element_rect(colour = "white", fill = "white"))
	return(g)
}