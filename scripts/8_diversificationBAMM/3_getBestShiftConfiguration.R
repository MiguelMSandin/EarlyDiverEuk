#---- 
#---- loading packages ----

library(BAMMtools)
library(coda)
library(data.table)
library(dplyr)
library(ggplot2)

#---- 
setwd("")
#----
#---- Organize the recently created files ----------------------------------------------------------

files = list(controls = grep("control_test.*ctl$", dir(recursive=TRUE), value=TRUE))

files$dirs = unique(sub("\\/bamm.*", "", files$controls))
files$clades = unique(files$controls %>% sub(".*control_test\\/", "", .) %>% sub("\\/.*", "", .))

for(d in files$dirs){
	for(g in files$clades){
		dir = paste0(d, "/bamm/diver/", g)
		if(!dir.exists(dir)){dir.create(dir)}
		system(paste0("mv ", d, "/bamm/diver/diversification_", g, "* ", dir))
	}
}; rm(d, g, dir)

#----
#---- Analyze diversification shifts from BAMM MCMC output chains ----------------------------------

files$mcmcs = grep("\\/diver\\/.*mcmc_out\\.txt\\.gz$", dir(recursive=TRUE), value=TRUE)
files$dirOutPlots = "../../plots/shifts/"
files$dirOutTables = "../../output/shifts/"
files$burn=25

for(i in 1:length(unlist(strsplit(files$dirOutPlots, "/")))){
	dir = paste(unlist(strsplit(files$dirOutPlots, "/"))[1:i], collapse="/")
	if(!dir.exists(dir)){dir.create(dir)}
};rm(i, dir)

for(i in 1:length(unlist(strsplit(files$dirOutTables, "/")))){
	dir = paste(unlist(strsplit(files$dirOutTables, "/"))[1:i], collapse="/")
	if(!dir.exists(dir)){dir.create(dir)}
};rm(i, dir)

# Start the loop for every clade and file
for(clade in files$clades){
	cat("  Working on ", clade, "\t(", grep(clade, files$clades), "/", length(files$clades), ")\n", sep="")
	pdf(paste0(files$dirOutPlots, "bestShiftConfiguration_", clade, ".pdf"), width=11.69, height=8.27, paper='special')
	outTable = data.frame()
	fileList = grep(clade, files$mcmcs, value=TRUE)
	for(file in fileList){
		cat("    Working on ", file, "\t(", grep(file, fileList), "/", length(fileList), ")\n", sep="")
		shifts = as.numeric(file %>% sub(".*shifts", "", .) %>% sub("_.*", "", .))
		div = file %>% sub(".*div", "", .) %>% sub("_.*", "", .)
		
		# Read the MCMC file
		mcmcout = read.csv(file, header=T)
		# Remove burnin
		postburn = mcmcout[((files$burn/100)*nrow(mcmcout)):nrow(mcmcout), ]
		
		# Start a plot of 4 grids
		# dev.off()
		par(mar=c(2.5,2.5,1,1))
		layout(matrix(c(1,6,2,4, 1,6,3,5), ncol=2),
			   heights=c(1, 1, 3, 3))
		plot.new()
		text(0.5,0.5, file, cex=2, font=2)
		
		# Check Effective Sampling Site
		essS = effectiveSize(postburn$N_shifts)
		essL = effectiveSize(postburn$logLik)
		
		# Check for convergence
		plot(mcmcout$logLik ~ mcmcout$generation,
			 main="Convergence", xlab="Generation", ylab="Log-Likelihood")
		
		# Plot the number of shifts
		# table(postburn$N_shifts) / nrow(postburn)
		hist(postburn$N_shifts,
			 xlim=c(0, max(postburn$N_shifts)),
			 main="Histogram of number of shifts occurrence", xlab="Number of shifts", ylab="Frequency")
		shiftsM = names(which.max(table(postburn$N_shifts)))
		
		# Compute Bayes Factor
		bayes = computeBayesFactors(mcmcout, expectedNumberOfShifts=shifts, burnin=files$burn/100)[,1]
		barplot(bayes,
				main="Bayes Factor", xlab="Number of shifts", ylab="Bayes Factor")
		bayesM = names(which.max(bayes))
		
		# Plot Prior Vs Posterior
		plotPrior(postburn, expectedNumberOfShifts=shifts)
		
		# Add subtitle with the summary of the statistics
		plot.new()
		text(0.5, 0.5, paste0("ESS shifts: ", round(essS, 1) , "; ESS LogLik:", round(essL, 1),
							  " - Most freq shift: ", shiftsM,
							  "; Shift with highest Bayes factor: ", bayesM),
			 cex=1.5, font=1)
		# Add row to table to export
		outTable = rbind(outTable, data.frame(file=file,
											  diversity=div,
											  shiftPrior=shifts,
											  mostFreqShift=shiftsM,
											  bayesShift=bayesM,
											  ESSshifts=essS,
											  ESSlogLik=essL))
	};rm(file, shifts, div, shiftsM, mcmcout, postburn, essS, essL, bayes, bayesM)
	dev.off()
	write.table(outTable, paste0(files$dirOutTables, "bestShiftConfiguration_", clade, ".tsv"), col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
}; rm(clade, fileList)

#----  
#---- Extract mean number of shifts over the 4 replicates ------------------------------------------

files$shiftFiles = grep("^bestShiftConfiguration_.*\\.tsv$", dir(files$dirOutTables, recursive=TRUE), value=TRUE)

data = data.frame()
for(file in files$shiftFiles){
	tmp = fread(paste0(files$dirOutTables, "/", file))
	tmp$clade = file %>% sub("^bestShiftConfiguration_", "", .) %>% sub("\\.tsv$", "", .)
	data = rbind(data, tmp)
}; rm(file, tmp)

summ = data %>% group_by(clade) %>% summarise(mean=mean(c(mostFreqShift, bayesShift)), sd=sd(c(mostFreqShift, bayesShift)))

write.table(summ, paste0(files$dirOutTables, "/shifts.tsv"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

datam = melt(data, id.vars=c("clade", "diversity", "shiftPrior", "file"))

(shiftsPlot = ggplot(subset(datam, variable=="mostFreqShift" | variable=="bayesShift"), aes(x=clade, y=value))+
		geom_hline(data=summ, aes(yintercept=mean), color="springgreen3", alpha=0.4)+
		# geom_point(alpha=0.6)+
		geom_jitter(width=0.2, height=0, alpha=0.6)+
		facet_wrap(~clade, scales="free")+
		theme_bw())

pdf(paste0(files$dirOutPlots, "/shiftsSummary.pdf"), width=11.69, height=8.27, paper='special'); plot(shiftsPlot); dev.off()

#----  
