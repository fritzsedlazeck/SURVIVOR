#! /usr/local/bin/Rscript

## This is the automated script for greating the plots over the summary stats using SURVIVOR. 



args<-commandArgs(TRUE)

if(length(args) != 1) {
	 cat("USAGE: produde_plots.R  summary_file_chr\n")
}else{

	if (!require("RColorBrewer")) {
		#install.packages("RColorBrewer")
		library(RColorBrewer)
	}

	file = args[[1]]
	file_chr = args[[2]]

	cols=(brewer.pal(5,"Set1"))
	#Plot2: 
	pdf(paste(file_chr,"_plot.pdf",sep="")) 
	t=read.table(file_chr,header=T)
	plot(t[,c(1,2)],ylim=c(0,max(t[,c(2:5)])),ylab="# of SVs",col=cols[1],xlab="chromosome",main=file_chr,type='points')
	points(t[,c(1,3)],col=cols[2])
	points(t[,c(1,4)],col=cols[3])
	points(t[,c(1,5)],col=cols[4])
	points(t[,c(1,6)],col=cols[5])
	legend('topright',legend=c('DEL','DUP','INV','INS','TRA'),lwd=2,col=cols)
	dev.off()
}
