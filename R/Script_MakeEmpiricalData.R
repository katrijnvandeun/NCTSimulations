setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls())

#------------------------------------------------------------------------------------------------------------------------------------------------------
# Open libraries
#------------------------------------------------------------------------------------------------------------------------------------------------------
install.packages("haven")
install.packages("psych")
install.packages("lavaan")
install.packages("NetworkComparisonTest")
install.packages("qgraph")
install.packages("tidyverse")
install.packages("glasso")
install.packages("MASS")
install.packages("mgm")

library("tidyverse")
library('dplyr')
library("haven") 
library("psych")
library("lavaan")
library("NetworkComparisonTest")
library("qgraph")
library("glasso")
library("MASS")
library("mgm")

#------------------------------------------------------------------------------------------------------------------------------------------------------
# Inspect data
#------------------------------------------------------------------------------------------------------------------------------------------------------

# load dataset

PROFIEL <- read_sav("Profiel_selectie2.sav")

Data1 <- dplyr::select(PROFIEL, tumor, Age_questionnaire, gesl, aantal_comorb, tijd_dx, chemotherapie, radiotherapie, FA, CF, SL, EF, PA, 
                DY, AP, CO, NV) #note: MASS package and dplyr clash

#change names
colnames(Data1) <-c("tumor","age","sex",
                   "comorb", "dx", "chemo",
                   "radio", "FA", "CF",
                   "SL", "EF", "PA", "DY",
                   "AP", "CO", "NV")

Datatot <- dplyr::select (Data1, FA, CF, SL, EF, PA, 
                DY, AP, CO, NV)

#------------------------------------------------------------------------------------------------------------------------------------------------------
# Datasets per tumorsoort.
#------------------------------------------------------------------------------------------------------------------------------------------------------

Colorectal <- subset(Data1, tumor==1)
Colorectaldata <- dplyr::select (Colorectal, FA, CF, SL, EF, PA, DY, AP, CO, NV)

Ovarian <- subset(Data1, tumor==5)
Ovariandata <- dplyr::select (Ovarian, FA, CF, SL, EF, PA, DY, AP, CO, NV)

Thyroid <- subset(Data1, tumor==6)
Thyroiddata <- dplyr::select (Thyroid, FA, CF, SL, EF, PA, DY, AP, CO, NV) 

Endometrial <- subset(Data1, tumor==7)
Endometrialdata <- dplyr::select (Endometrial, FA, CF, SL, EF, PA, DY, AP, CO, NV)

NHL <- subset(Data1, tumor==8)
NHLdata <- dplyr::select(NHL, FA, CF, SL, EF, PA, DY, AP, CO, NV)

HL <- subset(Data1, tumor==9)
HLdata <- dplyr::select(HL, FA, CF, SL, EF, PA, DY, AP, CO, NV)

CLL <- subset(Data1, tumor==10)
CLLdata <- dplyr::select(CLL, FA, CF, SL, EF, PA, DY, AP, CO, NV)

Breast <- subset(Data1, tumor==13)
Breastdata <- dplyr::select(Breast, FA, CF, SL, EF, PA, DY, AP, CO, NV)

#------------------------------------------------------------------------------------------------------------------------------------------------------
# Correlation matrix - normal distribution
#------------------------------------------------------------------------------------------------------------------------------------------------------

corMatTot <- cor(Datatot, use = "pairwise.complete.obs")
corMatTot

corMatColorectal <- cor(Colorectaldata, use = "pairwise.complete.obs")
corMatColorectal

corMatOvarian <- cor(Ovariandata, use = "pairwise.complete.obs")
corMatOvarian

corMatThyroid <- cor(Thyroiddata, use = "pairwise.complete.obs")
corMatThyroid

corMatEndometrial <- cor(Endometrialdata, use = "pairwise.complete.obs")
corMatEndometrial

corMatNHL <- cor(NHLdata, use = "pairwise.complete.obs")
corMatNHL

corMatHL <- cor(HLdata, use = "pairwise.complete.obs")
corMatHL

corMatCLL <- cor(CLLdata, use = "pairwise.complete.obs")
corMatCLL

corMatBreast <- cor(Breastdata, use = "pairwise.complete.obs")
corMatBreast

#------------------------------------------------------------------------------------------------------------------------------------------------------
# Correlation matrix - Polychoric correlation
#------------------------------------------------------------------------------------------------------------------------------------------------------

lavcorMatTot <- lavCor(Datatot, missing = "pairwise")
lavcorMatTot

lavcorMatColorectal <- lavCor(Colorectaldata, missing = "pairwise")
lavcorMatColorectal

lavcorMatOvarian <- lavCor(Ovariandata, missing = "pairwise")
lavcorMatOvarian

lavcorMatThyroid <- lavCor(Thyroiddata, missing = "pairwise")
lavcorMatThyroid

lavcorMatEndometrial <- lavCor(Endometrialdata, missing = "pairwise")
lavcorMatEndometrial

lavcorMatNHL <- lavCor(NHLdata, missing = "pairwise")
lavcorMatNHL

lavcorMatHL <- lavCor(HLdata, missing = "pairwise")
lavcorMatHL

lavcorMatCLL <- lavCor(CLLdata, missing = "pairwise")
lavcorMatCLL

lavcorMatBreast <- lavCor(Breastdata, missing = "pairwise")
lavcorMatBreast

#------------------------------------------------------------------------------------------------------------------------------------------------------
# Network graph
#------------------------------------------------------------------------------------------------------------------------------------------------------

#total
Graph_pcorTot <- qgraph(lavcorMatTot, graph = "pcor", layout = "spring", threshold = "bonferroni", sampleSize = nrow(Colorectaldata), alpha = 0.05, edge.labels = TRUE)
Graph_lassoTot <- qgraph(corMatTot, graph = "glasso", layout = "spring", tuning = 0.25, sampleSize = nrow(Datatot), edge.labels=TRUE)

#colorectal
Graph_lassoColorectal <- qgraph(corMatColorectal, graph = "glasso", layout = "spring", tuning = 0.25, sampleSize = nrow(Colorectaldata), edge.labels=TRUE)
Graph_lassoColorectal <- qgraph(lavcorMatColorectal, graph = "glasso", layout = "spring", tuning = 0.5, sampleSize = nrow(Colorectaldata), edge.labels=TRUE)

#ovarian
Graph_lassoOvarian <- qgraph(corMatOvarian, graph = "glasso", layout = "spring", tuning = 0.25, sampleSize = nrow(Ovariandata), edge.labels=TRUE)
Graph_lassoOvarian <- qgraph(lavcorMatOvarian, graph = "glasso", layout = "spring", tuning = 0.5, sampleSize = nrow(Ovariandata), edge.labels=TRUE)

#thyroid
Graph_lassoThyroid <- qgraph(corMatThyroid, graph = "glasso", layout = "spring", tuning = 0.25, sampleSize = nrow(Thyroiddata), edge.labels=TRUE)
Graph_lassoThyroid <- qgraph(lavcorMatThyroid, graph = "glasso", layout = "spring", tuning = 0.5, sampleSize = nrow(Thyroiddata), edge.labels=TRUE)

#endometrial
Graph_lassoEndometrial <- qgraph(corMatEndometrial, graph = "glasso", layout = "spring", tuning = 0.25, sampleSize = nrow(Endometrialdata), edge.labels=TRUE)
Graph_lassoEndometrial <- qgraph(lavcorMatEndometrial, graph = "glasso", layout = "spring", tuning = 0.5, sampleSize = nrow(Endometrialdata), edge.labels=TRUE)

#NHL
Graph_lassoNHL <- qgraph(corMatNHL, graph = "glasso", layout = "spring", tuning = 0.25, sampleSize = nrow(NHLdata), edge.labels=TRUE)
Graph_lassoNHL <- qgraph(lavcorMatNHL, graph = "glasso", layout = "spring", tuning = 0.25, sampleSize = nrow(NHLdata), edge.labels=TRUE)

#HL
Graph_lassoHL <- qgraph(corMatHL, graph = "glasso", layout = "spring", tuning = 0.25, sampleSize = nrow(HLdata), edge.labels=TRUE)
Graph_lassoHL <- qgraph(lavcorMatHL, graph = "glasso", layout = "spring", tuning = 0.5, sampleSize = nrow(HLdata), edge.labels=TRUE)

#CLL
Graph_lassoCLL <- qgraph(corMatCLL, graph = "glasso", layout = "spring", tuning = 0.25, sampleSize = nrow(CLLdata), edge.labels=TRUE)
Graph_lassoCLL <- qgraph(lavcorMatCLL, graph = "glasso", layout = "spring", tuning = 0.5, sampleSize = nrow(CLLdata), edge.labels=TRUE)

#Breast
Graph_lassoBreast<- qgraph(corMatBreast, graph = "glasso", layout = "spring", tuning = 0.25, sampleSize = nrow(Breastdata), edge.labels=TRUE)
Graph_lassoBreast<- qgraph(lavcorMatBreast, graph = "glasso", layout = "spring", tuning = 0.5, sampleSize = nrow(Breastdata), edge.labels=TRUE)

mgm_Breast <- mgm(Breastdata,type=c(rep("g",9),rep("c",3)),level = c(rep(1,9),2,2,2), lambdaSel = "EBIC", k=2)

mgm(Breastdata1,type=c(rep("g",9)),level = c(rep(1,9)), lambdaSel = "EBIC",k=2)

#------------------------------------------------------------------------------------------------------------------------------------------------------
#Network comparison test 
#------------------------------------------------------------------------------------------------------------------------------------------------------

#1. invariant network structure: tests whether structure is completely identical across groups
#2. invariant edge strength: tests whether edge strength of specific edge is equal across groups
#3. invariant global strenght: tests whether the overall level of connectivity is equal across groups

Colorectaldata1 <- Colorectaldata[complete.cases(Colorectaldata), ]
Ovariandata1 <- Ovariandata[complete.cases(Ovariandata), ]
Thyroiddata1 <- Thyroiddata[complete.cases(Thyroiddata), ]
Endometrialdata1 <- Endometrialdata[complete.cases(Endometrialdata), ]
NHLdata1 <- NHLdata[complete.cases(NHLdata), ]
HLdata1 <- HLdata[complete.cases(HLdata), ]
CLLdata1 <- CLLdata[complete.cases(CLLdata), ]
Breastdata1 <- Breastdata[complete.cases(Breastdata), ]

#Put data into list to facilitate looping over the data sets (Script_MCsimulation_size_power.R)
cancerdatasets = list(Breastdata1,CLLdata1,Colorectaldata1,Endometrialdata1,HLdata1,NHLdata1,Ovariandata1,Thyroiddata1)

NCTcolorovar <- NCT(Colorectaldata1, Ovariandata1, binary.data=FALSE, paired=FALSE, test.edges=TRUE, edges="all")
NCTCLLBreast <- NCT(CLLdata1, Breastdata1, binary.data=FALSE, paired=FALSE, test.edges=TRUE, edges="all")
NCTCLLBreast