module load R/3.2.5
R

setwd("/well/mccarthy/users/agata/MAVE_HNF1A/TWIST-Jan2020/")

codon_files=list.files("codon_summary/")

twist=list()
for (f in codon_files){
	name=gsub(".summary","",f)
	twist[[name]]=read.table(paste0("codon_summary/",f)) 
}

### remove all the variants where there is no actual NT change
twist=lapply(twist, function(x) x[which(as.character(x$V4) != as.character(x$V6)),])

x_axis=1:200
col=rainbow(20)
pdf("TWIST-Jan2020.codon_coverage.pdf", width=10)
par(mar=c(5.1,5.1,5.1,1))
for (i in 1:length(twist)){
	design_cov = sapply(x_axis, function(x) nrow(twist[[i]][with(twist[[i]], V9==1 & V8>=x),]))
	nondesign_cov = sapply(x_axis, function(x) nrow(twist[[i]][with(twist[[i]], V9==0 & V8>=x),]))
	plot(x_axis, design_cov, type="l", col="darkred", xlab="Coverage depth", 
		ylab="# variants detected", cex.axis=1.5, cex.lab=1.5, main=names(twist)[i], ylim=c(0, max(c(design_cov, nondesign_cov))))
	lines(x_axis, nondesign_cov, type="l", col="coral")
	legend("topright", c("in design","not in design"), col=c("darkred","coral"), lty=1, bty="n")

}
dev.off()
