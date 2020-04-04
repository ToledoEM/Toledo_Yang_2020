library(circlize)
library(dplyr)

# This file contains the code and data to replicate the plots of the figures of supplementary 2.



# Read input data and prepare data for plot
bedlxr <- read.delim("Figure_Sup2/peaks lxr control.bed", header=F,col.names = c("chr","start","end"))
bed.gw <- read.delim("Figure_Sup2/Peaks LXR Treated.bed",  header=F,col.names = c("chr","start","end"))
bed.top.gw <-read.delim("Figure_Sup2/top10peaks_withID_LXR_GW.bed", header=F,col.names = c("chr","start","end","ID")) %>% select(-ID)
bed.top.control <- read.delim("Figure_Sup2/LXR_top10_control.bed", header=F, col.names = c("chr","start","end","ID")) %>% select(-ID)


bed <- rbind(bedlxr,bed.gw) %>% mutate(value=c(rep(1,nrow(bedlxr)),rep(3,nrow(bed.gw))))


bed_list <- list(bedlxr, bed.gw)
bed_list_top <- list(bed.top.control,bed.top.gw)




### MAIN S2A

par(mar = c(1, 1, 1, 1))
circos.initializeWithIdeogram(species = "mm9") # all chr mm9
circos.genomicDensity(bedlxr,type = "h",lwd = 0.5,lty = 2,col = "gray50",border = NA,baseline = 0)
circos.genomicDensity(bed.gw,type = "h",lwd = 0.5,col = "red",lty = 2,border = NA,baseline = 0)
circos.genomicTrackPlotRegion(bed_list_top, stack = TRUE,
                              panel.fun = function(region, value, ...) {
                                i = getI(...)
                                circos.genomicLines(region, value, col = i, ...)
                              })

circos.clear()


# Chr 15
par(mar = c(1, 1, 1, 1))
circos.par("start.degree" = 150, "gap.degree" = c(240), "cell.padding" = c(0, 0, 0, 0))
circos.initializeWithIdeogram(species = "mm9",chromosome.index = "chr15") # all chr mm9
circos.genomicDensity(bedlxr,type = "h",lwd = 0.5,lty = 2,col = "gray50",border = NA,baseline = 0)
circos.genomicDensity(bed.gw,type = "h",lwd = 0.5,col = "red",lty = 2,border = NA,baseline = 0)
circos.genomicTrackPlotRegion(bed_list_top, stack = TRUE,
                            panel.fun = function(region, value, ...) {
                                i = getI(...)
                                circos.genomicLines(region, value, col = i, ...)
                              })
circos.clear()


## Chr 11
par(mar = c(1, 1, 1, 1))
circos.par("start.degree" = 150, "gap.degree" = c(240), "cell.padding" = c(0, 0, 0, 0))
circos.initializeWithIdeogram(species = "mm9",chromosome.index = "chr11") # all chr mm9
circos.genomicDensity(bedlxr,type = "h",lwd = 0.5,lty = 2,col = "gray50",border = NA,baseline = 0)
circos.genomicDensity(bed.gw,type = "h",lwd = 0.5,col = "red",lty = 2,border = NA,baseline = 0)
circos.genomicTrackPlotRegion(bed_list_top, stack = TRUE,
                              panel.fun = function(region, value, ...) {
                                i = getI(...)
                                circos.genomicLines(region, value, col = i, ...)
                              })

circos.clear()


## Chr 5
par(mar = c(1, 1, 1, 1))
circos.par("start.degree" = 150, "gap.degree" = c(240), "cell.padding" = c(0, 0, 0, 0))
circos.initializeWithIdeogram(species = "mm9",chromosome.index = "chr5") # all chr mm9
circos.genomicDensity(bedlxr,type = "h",lwd = 0.5,lty = 2,col = "gray50",border = NA,baseline = 0)
circos.genomicDensity(bed.gw,type = "h",lwd = 0.5,col = "red",lty = 2,border = NA,baseline = 0)
circos.genomicTrackPlotRegion(bed_list_top, stack = TRUE,
                              panel.fun = function(region, value, ...) {
                                i = getI(...)
                                circos.genomicLines(region, value, col = i, ...)
                              })

circos.clear()
