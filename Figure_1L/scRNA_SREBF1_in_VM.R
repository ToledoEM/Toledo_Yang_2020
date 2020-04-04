# MCMC Srebf1 expression in VM scRNAseq data from La Manno 2016 -----


# From already Calculated model ------
library(tidyverse)

MAP <- readRDS("Figure_1L/OUTPUT/MAP_Srebf1.rds")

tMAP <- MAP %>% as_tibble() %>% gather(key = "Cell_type",value = "beta")


#reorder the levels
TGTCELL <- unique(tMAP$Cell_type)

TGTCELL <- factor(TGTCELL,levels = c("Sz","mMgl","mPeric","mEndo","mEpend","mRgl3","mRgl2","mRgl1","mNProg","mNbM",
                                     "mNbDA","mDA0","mDA1","mDA2","mNbML1","mNbML2","mNbML3","mNbML4","mNbML5","mGaba1a","mGaba1b","mGaba2",
                                     "mNbL1","mNbL2","mOMTN","mRN","mSert"  ))


order <- data.frame(Cell_type=levels(TGTCELL),order=seq(1,length(TGTCELL)))
tMAP$COLOR <- ifelse(tMAP$Cell_type %in% c("mRgl1","mRgl2","mRgl3")==T,yes = 1,0)
tMAP <- tMAP %>% inner_join(order,by="Cell_type")


plot1 <- ggplot(tMAP, aes(x=reorder(Cell_type,order), y=beta, fill=factor(COLOR))) + geom_violin(trim=F,scale = "width")+theme_classic()+
  labs(x=NULL,y="Srebf1 Molecules per cell")+scale_fill_manual(values = c("grey80","red"))+geom_boxplot(fill="white",width=0.25,notch = T,outlier.size = 0,lwd=0.1)+
  theme(axis.text.x = element_text(angle = 30, hjust = 1),legend.position="none")

plot1

ggsave(plot1,filename = "Figure_1L/Srebf1_SC.pdf")

# MCMC Srebf1 calculations / need testing -----

# Srebf1 figure
library(ggplot2)
library(dplyr)
# needs
library(pipeR)
library(rstan)
library(foreach)
library(doParallel)

# You will need the fucntion from https://github.com/ToledoEM/BayesianGLM/tree/master/src 

# This file is from https://github.com/linnarsson-lab/ipynb-lamanno2016/tree/master/data
# https://github.com/linnarsson-lab/ipynb-lamanno2016/raw/master/data/Mouse_Embryo_nound.cef

# source("https://raw.githubusercontent.com/ToledoEM/BayesianGLM/master/src/BayesianGLM.R")
# source("https://raw.githubusercontent.com/ToledoEM/BayesianGLM/master/src/Ceftools.R")
# 
dir.create("Figure1L/OUTPUT")

# DATASET <- "Mouse_Embryo_nound.cef"
DATASET <- "https://raw.githubusercontent.com/linnarsson-lab/ipynb-lamanno2016/master/data/Mouse_Embryo_nound.cef"

# We need to do the MCMC for this gene.
# the details of this procedure it is in https://github.com/ToledoEM/BayesianGLM

GENE_FIG <- c("Srebf1")
BayesianGLM(DATASET, GENE_FIG, parallel=F) #one gene, so no parallel
