#load all packages for sorting, general exploration, analyses, and figures
library(dplyr)
library(plyr)
library(tidyr)
library(car)
library(lmerTest)
library(lme4)
library(MASS)
library(ggplot2)
library(ggExtra)
library(gridExtra)
library(scales)
library(nlme)
library(plyr)
library(rpart)
library(randomForest)
library(gbm)
library(dismo)
library(caret)
library(mlbench)
library(vip)
library(pdp)
library(randomForestExplainer)
library(vivid)
library(sna)
library(vivid)
library(intergraph)
library(cowplot)
library(VSURF) 
library(ggRandomForests)
library(randomForestSRC)

#data are available at https://doi.org/10.6073/pasta/c44c5c83d9dbd1179b4bd1f7e779a7f1


#set correct WD

getwd()
setwd('enter directory here')

#reead in CSV file
UMRS_HSI<-read.csv( "UMRS_HSI_data.csv", header = T)

#setseed
# Setting a seed
set.seed(123)

#meanHSI VSURF process
hsi.vsurf <- VSURF(UMRS_HSI[,15:29], UMRS_HSI[,14], mtry = 100)

#view selected variable under predicted category in VSURF package
plot(hsi.vsurf, step = "pred", imp.sd = FALSE, var.names = TRUE) 

#run random forest model with variables selected in VSURF
HSIo_RF <- rfsrc(MEAN_HSIo~NUM_OUTL+FDD+PCT_OUTL+CHLa+AVG_FETCH+DIS_Z_SCORE+WDP+PCT_AQVEG+ZICE+EFF_CONN, data=UMRS_HSI, importance = TRUE) 
print(HSIo_RF)


#prepaing and and producing partial dependence plots for random forest model of HSIo

#create theme for figures
my_theme = theme(axis.title.x = element_text(size = 12, face = "bold"),
                 axis.text.x = element_text(size = 10, face = "bold"),
                 axis.title.y = element_text(size = 12, face = "bold"),
                 axis.text.y = element_text(size = 10))

#prepare interface between RFSRC and GGplot
varsel_hsi<-var.select(HSIo_RF)

gg_md<-gg_minimal_vimp(varsel_hsi)

plot(gg_md)

gg_hsi<-gg_variable(HSIo_RF)

gg_v<-gg_variable(HSIo_RF)

xvar<-gg_md$names

plot(gg_v, panel = TRUE)



partial_hsio2<-plot.variable(HSIo_RF,
                             xvar = xvar,
                             partial = TRUE, sorted = TRUE,
                             show.plots = FALSE)

gg_p <- gg_partial(partial_hsio2)

#produce partial dependence plots for random forest model of HSIo

hsio_partialplot<-plot(gg_p, panel=TRUE)+my_theme+stat_smooth(se =FALSE,  size = 1, method = loess)+geom_point(size = 1.5)+ labs(x="", y = "Partial effect on HSIo")

hsio_partialplot


#Water temperature VSURF process

temp.vsurf <- VSURF(UMRS_HSI[,15:29], UMRS_HSI[,11], mtry = 100)


#view selected variable under predicted category in VSURF package

plot(temp.vsurf, step = "pred", imp.sd = FALSE, var.names = TRUE)

#run random forest model with variables selected in VSURF

water_temp_RF <- rfsrc(WATER_TEMP~FDD+PCT_OUTL+DIS_Z_SCORE, data=UMRS_HSI, importance = TRUE) 
print(water_temp_RF)

#prepare interface between RFSRC and GGplot
varsel_temp<-var.select(water_temp_RF)

gg_md_temp<-gg_minimal_vimp(varsel_temp)

plot(gg_md_temp)

gg_temp<-gg_variable(water_temp_RF)

gg_v_temp<-gg_variable(water_temp_RF)

xvar_temp<-gg_md_temp$names

plot(gg_v_temp, panel = TRUE)



partial_temp<-plot.variable(water_temp_RF,
                            xvar = xvar,
                            partial = TRUE, sorted = TRUE,
                            show.plots = FALSE)

gg_p_temp <- gg_partial(partial_temp)

#produce partial dependence plots for random forest model of water temperature

temp_partialplot<-plot(gg_p_temp, panel=TRUE)+my_theme+stat_smooth(se =FALSE,  size = 1, method = loess)+geom_point(size = 1.5)+ labs(x="", y = "Partial effect on temp. (C)")

temp_partialplot


#Dissolved oxygen VSURF process
DO.vsurf <- VSURF(UMRS_HSI[,15:29], UMRS_HSI[,12], mtry = 100)

#view selected variable under predicted category in VSURF package
plot(DO.vsurf, step = "pred", imp.sd = FALSE, var.names = TRUE)

#run random forest model with variables selected in VSURF
DO_RF <- rfsrc(DO~FDD+CHLa+ZSNOW+DIS_Z_SCORE+PCT_OUTL+ZICE+PCT_AQVEG+AVG_FETCH+WDP+SDI, data=UMRS_HSI, importance = TRUE) 
print(DO_RF)

#prepare interface between RFSRC and GGplot
varsel_DO<-var.select(DO_RF )

gg_md_DO<-gg_minimal_vimp(varsel_DO)

plot(gg_md_DO)

gg_DO<-gg_variable(DO_RF )

gg_v_DO<-gg_variable(DO_RF )

xvar_DO<-gg_md_DO$names

plot(gg_v_DO, panel = TRUE)



partial_DO<-plot.variable(DO_RF,
                          xvar = xvar_DO,
                          partial = TRUE, sorted = TRUE,
                          show.plots = FALSE)

gg_p_DO <- gg_partial(partial_DO)

#produce partial dependence plots for random forest model of dissolved oxygen
do_partialplot<-plot(gg_p_DO, panel=TRUE)+my_theme+stat_smooth(se =FALSE,  size = 1, method = loess)+geom_point(size = 1.5)+ labs(x="", y = "Partial effect on DO (mg/L)")

do_partialplot

#velocity VSURF process
vel.vsurf <- VSURF(UMRS_HSI[,15:29], UMRS_HSI[,13], mtry = 100)

#view selected variable under predicted category in VSURF package
plot(vel.vsurf, step = "pred", imp.sd = FALSE, var.names = TRUE)

#run random forest model with variables selected in VSURF
vel_RF <- rfsrc(VEL_CM_S~FDD+WDP+PCT_OUTL+DIS_Z_SCORE, data=UMRS_HSI, importance = TRUE) 
print(vel_RF)

#prepare interface between RFSRC and GGplot
varsel_VEL<-var.select(vel_RF)

gg_md_VEL<-gg_minimal_vimp(varsel_VEL)

plot(gg_md_VEL)

gg_VEL<-gg_variable(vel_RF)

gg_v_VEL<-gg_variable(vel_RF)

xvar_VEL<-gg_md_VEL$names

plot(gg_v_VEL, panel = TRUE)



partial_VEL<-plot.variable(vel_RF,
                           xvar = xvar,
                           partial = TRUE, sorted = TRUE,
                           show.plots = FALSE)

gg_p_VEL <- gg_partial(partial_VEL)

#produce partial dependence plots for random forest model of flow velocity

vel_partialplot<-plot(gg_p_VEL, panel=TRUE)+my_theme+geom_smooth(se =FALSE,  size = 1, method = loess)+geom_point(size = 1.5)+ labs(x="", y = "Partial effect on velocity (cm/s)")

vel_partialplot


