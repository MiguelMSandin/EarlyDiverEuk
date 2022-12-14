
# .rs.restartR()
library(ape)
library(BAMMtools)
library(HDInterval)
library(dplyr)
library(data.table)
library(ggplot2)

# Set file names ___________________________________________________________________________________
files <- list(clade="Holozoa",
              tree="clade_Holozoa.tre",
              eventDataDir=".")

(files$eventDataFiles = grep("event_data.txt", dir(files$eventDataDir), value=TRUE))

# Open the anayzed tree ____________________________________________________________________________
tree <- read.tree(files$tree)
# plot(tree, cex=.6)
# axisPhylo()

# Start the loop ___________________________________________________________________________________
i = 0
for(eventData in files$eventDataFiles){
    
    shifts = eventData %>% sub(".*shifts", "", .) %>% sub("_.*", "", .)
    div = eventData %>% sub(".*div", "", .) %>% sub("_.*", "", .)
    
    cat("  Working on ", eventData, " (", (i = i+1), "/", length(files$eventDataFiles), ")\n", sep="")
    cat("    Diversity estimate: ", div, "%\n", sep="")
    cat("    Expected number of shifts:", shifts, "\n")
    
    outTree = paste0(files$eventDataDir, "/", sub("_event_data.txt", "_tree.pdf", eventData))
    outTable = paste0(files$eventDataDir, "/", sub("_event_data.txt", "_RTT.tsv", eventData))
    outRTT = paste0(files$eventDataDir, "/", sub("_event_data.txt", "_RTT.pdf", eventData))
    
    cat("  Reading Event Data\n") # ________________________________________________________________
    edata <- getEventData(tree, 
                          eventdata = paste0(files$eventDataDir, "/", eventData), 
                          burnin=0.25)
    
    cat("  Plotting the tree\n") # _________________________________________________________________
    pdf(outTree, width=11.69, height=8.27, paper='special')
    bamm.tree = plot.bammdata(edata, lwd=2, labels=TRUE, cex=0.2, logcolor=TRUE)
    addBAMMshifts(edata, cex=2)
    addBAMMlegend(bamm.tree)
    title(main = files$clade,
          sub = paste0("Speciation-extinction analysis with an estimated ", div, "% diversity and ", shifts, " expected shifts."))
    dev.off()
    
    cat("  Extracting RTT matrix\n") # _____________________________________________________________
    edatam <- getRateThroughTimeMatrix(edata)
    
    # Create a table for 'speciation' data
    tmp <- as.data.frame(t(edatam$lambda))
    colnames(tmp) <- paste0("replicate", seq(1:ncol(tmp)))
    edatam_spe <- cbind(time=rev(edatam$times), tmp)
    
    # Create a table for 'extinction' data
    tmp <- as.data.frame(t(edatam$mu))
    colnames(tmp) <- paste0("replicate", seq(1:ncol(tmp)))
    edatam_ext <- cbind(time=rev(edatam$times), tmp)
    
    # Create a table for 'diversification' data
    tmp <- as.data.frame(t(edatam$lambda)) - as.data.frame(t(edatam$mu))
    colnames(tmp) <- paste0("replicate", seq(1:ncol(tmp)))
    edatam_div <- cbind(time=rev(edatam$times), tmp)
    
    # first summarise 'speciation' data
    hpd = apply(edatam_spe[,-1], 1, hdi)
    edatamOut = data.frame(rate="speciation",
                           time=edatam_spe$time,
                           median=apply(edatam_spe[,-1], 1, median),
                           HPDlow=hpd[1,],
                           HPDupp=hpd[2,],
                           mean=apply(edatam_spe[,-1], 1, mean),
                           sd=apply(edatam_spe[,-1], 1, sd))
    
    # now summarise 'extinction' data
    hpd = apply(edatam_ext[,-1], 1, hdi)
    edatamOut = rbind(edatamOut, data.frame(rate="extinction",
                                            time=edatam_ext$time,
                                            median=apply(edatam_ext[,-1], 1, median),
                                            HPDlow=hpd[1,],
                                            HPDupp=hpd[2,],
                                            mean=apply(edatam_ext[,-1], 1, mean),
                                            sd=apply(edatam_ext[,-1], 1, sd)))
    
    # and finally summarise 'diversification' data
    hpd = apply(edatam_div[,-1], 1, hdi)
    edatamOut = rbind(edatamOut, data.frame(rate="diversification",
                                            time=edatam_div$time,
                                            median=apply(edatam_div[,-1], 1, median),
                                            HPDlow=hpd[1,],
                                            HPDupp=hpd[2,],
                                            mean=apply(edatam_div[,-1], 1, mean),
                                            sd=apply(edatam_div[,-1], 1, sd)))
    
    cat("  Exporting summary file\n") # ____________________________________________________________
    write.table(edatamOut, outTable, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
    
    cat("  Plotting RTT slopes\n") # _______________________________________________________________
    # Beautifying the table for plotting
    edatamOut$rate = factor(edatamOut$rate, levels=c("extinction", "speciation", "diversification"))
    
    # Plotting
    (rttPlot = ggplot(edatamOut, aes(x=-time, y=median, colour=rate))+
            geom_line()+
            geom_ribbon(aes(ymin=HPDlow, ymax=HPDupp, fill=rate, color=NULL), alpha=0.1)+
            scale_x_continuous(breaks=seq(-max(round(edatamOut$time, -2)), 0, 100),
                               minor_breaks=seq(-max(round(edatamOut$time, -2)), 0, 50)) +
            labs(title=paste0(files$clade),
                 subtitle=paste0("Speciation-extinction analysis with an estimated ", div, "% diversity and ", shifts, " expected shifts."),
                 y="Rates Through Time", x="Time before present")+
            theme_bw())
    
    # Export plot
    pdf(outRTT, width=11.69, height=8.27, paper='special')
    plot(rttPlot)
    dev.off()
    
    cat("  Done\n")
    cat("____________________\n")
}; rm(i, eventData, shifts, div, edata, outTree, outRTT, outTable, bamm.tree, edatam_spe, edatam_ext, edatam_div, hpd, edatamOut, tmp, rttPlot)


# Cat all RTT slopes into a single file

system(paste0("pdftk ", paste0(files$eventDataDir, grep("RTT\\.pdf", dir(files$eventDataDir), value=TRUE), collapse=" "),
              " cat output ",
              sub(files$clade, "", files$eventDataDir) ,"diversification_", files$clade, "_RTT.pdf", collapse=""))

# Cat all trees into a single file
system(paste0("pdftk ", paste0(files$eventDataDir, grep("tree\\.pdf", dir(files$eventDataDir), value=TRUE), collapse=" "),
              " cat output ",
              sub(files$clade, "", files$eventDataDir) ,"diversification_", files$clade, "_tree.pdf", collapse=""))

#

