#!/bin/env Rscript


## This is the automated script for greating the plots over the summary stats using SURVIVOR. 
if (!require("RColorBrewer")) {
	#install.packages("RColorBrewer")
	library(RColorBrewer)
}


args<-commandArgs(TRUE)

if(length(args) != 2) {
	 cat("USAGE: produde_plots.R summary_file summary_file_chr\n")
}else{
	file = args[[1]]
	file_chr = args[[2]]
	##paste(file,"_plot.pdf",sep="")

	cols=(brewer.pal(5,"Set1"))
	#Plot1: length of the different SVs
	pdf(paste(file,"_plot.pdf",sep=""))
	t=read.table(file,header=T)
	plot(log10(t[,1]),t[,2],xlab="log10(length(bp))",ylab="# of SVs",type='l',col=cols[1],main=file)
	lines(log10(t[,1]),t[,3],col=cols[2])
	lines(log10(t[,1]),t[,4],col=cols[3])
	lines(log10(t[,1]),t[,5],col=cols[4])
	legend('topright',legend=c('DEL','DUP','INV','INS'),lwd=2,col=cols)
	dev.off()
	
	#Plot2: 
	pdf(paste(file_chr,"_plot.pdf",sep="")) 
	t=read.table(file_chr,header=T)
	plot(c(1:length(t[,1])),t[,2],ylim=c(0,max(t[,c(2:5)])),ylab="# of SVs",col=cols[1],xlab="chromosome",main=file_chr)
	points(c(1:length(t[,1]))-0.01,t[,3],col=cols[2])
	points(c(1:length(t[,1]))+0.01,t[,4],col=cols[3])
	points(c(1:length(t[,1]))+0.05,t[,5],col=cols[4])
	points(c(1:length(t[,1]))+0.05,t[,6],col=cols[5])
	legend('topright',legend=c('DEL','DUP','INV','INS','TRA'),lwd=2,col=cols)
	dev.off()
}
