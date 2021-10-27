library(tidyverse)
library(stringi)
library(stringr)
library(pscl)
library(ROCR)
library(lmtest)
library(MASS)
library(usdm)
library(classInt)
library(dismo)
library(gbm)
library(raster)
library(rgdal)
library(gridExtra)
library(sf)
library(dggridR)
library(corrplot)
library(GGally)
library(mgcv)
library(gstat)
library(ape)
library(data.table)

options(scipen=999)

dat.use <- read.csv("CONI_CleanDataForAnalysis.csv") %>% 
  dplyr::select(-starts_with("ag")) %>% 
  unique()

#1. Set parameters----
scales <- data.frame(response=c("boom"),
                     scale=c(200))
layers <- c("conifer", "decid", "fire", "gravel", "harvest", "industry", "mixed", "moisture", "nutrient", "pine", "roads", "seismic", "water", "wells", "wetland")

#2. Load rasters----

tifs <- data.frame(file = list.files("/Volumes/SSD/GIS/Projects/Scale/3MovingWindow", pattern="*00.tif")) %>% 
  separate(file, into=c("cov", "scale", "tif"), remove=FALSE) %>% 
  mutate(var=paste0(cov, "_", scale)) %>% 
  dplyr::filter(scale %in% scales$scale) %>% 
  unique()

#100m resolution
raster.100 <- list()
for(i in 1:nrow(tifs)){
  raster.i <- raster(as.character(paste0("/Volumes/SSD/GIS/Projects/Scale/3MovingWindow/", tifs$file[i])))
  resamp <- raster(raster.i)
  res(resamp) <- c(100, 100)
  raster.low.i <- raster::resample(x=raster.i, y=resamp, method="bilinear")
  names(raster.low.i) <- tifs$var[i]
  raster.100[i] <- raster.low.i
}

raster.100.stack <- stack(raster.100)
rm(raster.100)

#30 m resolution
raster.30 <- list()
for(i in 1:nrow(tifs)){
  raster.i <- raster(as.character(paste0("/Volumes/SSD/GIS/Projects/Scale/3MovingWindow/", tifs$file[i])))
  names(raster.i) <- tifs$var[i]
  raster.30[i] <- raster.i
}

raster.30.stack <- stack(raster.30)
rm(raster.30)

#3. Set up loop----
boot <- 10
samples <- 8

set.seed(1234)

#Start dataframes
brt.overall.pred <- list()
brt.overall.eval <- data.frame()
brt.overall.eval.df <- data.frame()

for(j in 1:nrow(scales)){
  
  time <- Sys.time()
  
  scale.j <- scales$scale[j]
  response.j <- scales$response[j]
  layers.j <- c("ID", "cell", "JULIAN", "TOD", "pres.peent", "pres.boom", "pres.notboom", paste0(layers, "_", scale.j))
  dat.j <- dat.use %>% 
    dplyr::select(layers.j)
  
  for(i in 1:boot){
    
    time <- Sys.time()
    
    #4. Select data----
    
    #Randomly select one location per grid cell
    locs.i <- dat.j %>% 
      dplyr::select(ID, cell) %>% 
      unique() %>% 
      group_by(cell) %>% 
      sample_n(size = 1) %>% 
      ungroup()
    
    #Randomly select 8 recordings per selected ID
    recs.i <- dat.j %>% 
      dplyr::filter(ID %in% locs.i$ID) %>% 
      group_by(ID) %>% 
      sample_n(size = samples) %>% 
      ungroup()
    
    #Summarize data by peents & booms
    sum.i <- recs.i %>% 
      group_by(ID) %>% 
      summarize(sum.peent = ifelse(sum(pres.peent) > 0, 1, 0),
                sum.boom = ifelse(sum(pres.boom) > 0, 1, 0),
                sum.notboom = ifelse(sum(pres.notboom) > 0, 1, 0)) %>% 
      ungroup()
    
    #Join back to covariates
    dat.i <- recs.i %>% 
      dplyr::select(-pres.peent, -pres.boom, -pres.notboom, -JULIAN, -TOD) %>% 
      unique() %>% 
      right_join(sum.i) %>% 
      data.frame() %>% 
      mutate(sum.peent = as.integer(sum.peent),
             sum.boom = as.integer(sum.boom))
    
    #5. Run model----
    
    #brt
    if(response.j=="peent"){
      brt.overall.i <- dismo::gbm.step(data=dat.i, 
                                       gbm.x=3:17,
                                       gbm.y=18,
                                       family="bernoulli",
                                       tree.complexity = 3,
                                       learning.rate = 0.001,
                                       bag.fraction = 0.75,
                                       max.trees=10000,
                                       verbose=FALSE) 
    }
    else
    {
      brt.overall.i <- dismo::gbm.step(data=dat.i, 
                                       gbm.x=3:17,
                                       gbm.y=19,
                                       family="bernoulli",
                                       tree.complexity = 3,
                                       learning.rate = 0.001,
                                       bag.fraction = 0.75,
                                       max.trees=10000,
                                       verbose=FALSE) 
    }
    
    elapsed <- Sys.time() - time
    
    print(paste0("******COMPLETED BRT ", i, " OF ", boot, " BRTS FOR SCALE ", scale.j, " in ", elapsed, " seconds******"))
    
    #6. Predict----
    
    #Predict 100 m
    time <- Sys.time()
    
    brt.overall.pred.100 <- dismo::predict(raster.100.stack, brt.overall.i, n.trees=brt.overall.i$gbm.call$best.trees, type="response")
    
    writeRaster(brt.overall.pred.100, paste0("/Volumes/SSD/GIS/Projects/Scale/5Predictions/OverallBRTPredictions_100_", response.j, "_", scale.j, "_", i, ".tif"), format="GTiff", overwrite=TRUE)
    
    elapsed <- Sys.time() - time
    
    print(paste0("******COMPLETED 100 m PREDICTION ", i, " OF ", boot, " PREDICTIONS FOR SCALE ", scale.j, " in ", elapsed, " minutes******"))
    
    #Predict 30 m
    time <- Sys.time()
    
    brt.overall.pred.30 <- dismo::predict(raster.30.stack, brt.overall.i, n.trees=brt.overall.i$gbm.call$best.trees, type="response")
    
    writeRaster(brt.overall.pred.30, paste0("/Volumes/SSD/GIS/Projects/Scale/5Predictions/OverallBRTPredictions_30_", response.j, "_", scale.j, "_", i, ".tif"), format="GTiff", overwrite=TRUE)
    
    elapsed <- Sys.time() - time
    
    print(paste0("******COMPLETED 30 m PREDICTION ", i, " OF ", boot, " PREDICTIONS FOR SCALE ", scale.j, " in ", elapsed, " hours******"))
    
    time <- Sys.time()
    
    #Evaluate----
    
    #Get test data
    locs.test <- dat.use %>% 
      dplyr::select(ID, cell) %>% 
      unique() %>% 
      anti_join(locs.i) %>% 
      group_by(cell) %>% 
      sample_n(size = 1) %>% 
      ungroup()
    
    recs.test <- dat.use %>% 
      dplyr::filter(ID %in% locs.test$ID) %>% 
      group_by(ID) %>% 
      sample_n(size = samples) %>% 
      ungroup()
    
    sum.test <- recs.test %>% 
      group_by(ID) %>% 
      summarize(sum.peent = ifelse(sum(pres.peent) > 0, 1, 0),
                sum.boom = ifelse(sum(pres.boom) > 0, 1, 0)) %>% 
      ungroup()
    
    if(response.j=="peent"){
      dat.test <- recs.test %>% 
        right_join(sum.test) %>% 
        data.frame() %>% 
        mutate(sum.response = as.integer(sum.peent)) %>% 
        dplyr::select(ID, X, Y, sum.response) %>% 
        unique()
    }
    else{
      dat.test <- recs.test %>% 
        right_join(sum.test) %>% 
        data.frame() %>% 
        mutate(sum.response = as.integer(sum.boom)) %>% 
        dplyr::select(ID, X, Y, sum.response) %>% 
        unique()
    }
    
    #Evaluate
    overall.suit.100 <- data.frame(suitability = raster::extract(brt.overall.pred.100, dat.test[,c("X", "Y")]))
    overall.val.100 <- cbind(dat.test, overall.suit.100)
    overall.e.100 <- dismo::evaluate(subset(overall.val.100, sum.response==1)$suitability, subset(overall.val.100, sum.response==0)$suitability)
    overall.eval.100 <- data.frame(boot=i, mode=response.j)
    overall.eval.100$np <- overall.e.100@np
    overall.eval.100$na <- overall.e.100@na
    overall.eval.100$auc <- overall.e.100@auc
    overall.eval.100$cor <- overall.e.100@cor
    overall.eval.100$cor <- overall.e.100@pcor
    overall.eval.100$odp <- mean(overall.e.100@ODP)
    overall.eval.100$scale <- scale.j
    overall.eval.100$extent <- 100
    overall.eval.df.100 <- data.frame(t = overall.e.100@t) %>% 
      mutate(prev = overall.e.100@prevalence,
             ccr = overall.e.100@CCR,
             tpr = overall.e.100@TPR,
             tnr = overall.e.100@TNR,
             fpr = overall.e.100@FPR,
             fnr = overall.e.100@FNR,
             ppp = overall.e.100@PPP,
             npp = overall.e.100@NPP,
             mcr = overall.e.100@MCR,
             or = overall.e.100@OR,
             kappa = overall.e.100@kappa,
             boot = i,
             mode=response.j, 
             extent=100)
    
    overall.suit.30 <- data.frame(suitability = raster::extract(brt.overall.pred.30, dat.test[,c("X", "Y")]))
    overall.val.30 <- cbind(dat.test, overall.suit.30)
    overall.e.30 <- dismo::evaluate(subset(overall.val.30, sum.response==1)$suitability, subset(overall.val.30, sum.response==0)$suitability)
    overall.eval.30 <- data.frame(boot=i, mode=response.j)
    overall.eval.30$np <- overall.e.30@np
    overall.eval.30$na <- overall.e.30@na
    overall.eval.30$auc <- overall.e.30@auc
    overall.eval.30$cor <- overall.e.30@cor
    overall.eval.30$cor <- overall.e.30@pcor
    overall.eval.30$odp <- mean(overall.e.30@ODP)
    overall.eval.30$scale <- scale.j
    overall.eval.30$extent <- 30
    overall.eval.df.30 <- data.frame(t = overall.e.30@t) %>% 
      mutate(prev = overall.e.30@prevalence,
             ccr = overall.e.30@CCR,
             tpr = overall.e.30@TPR,
             tnr = overall.e.30@TNR,
             fpr = overall.e.30@FPR,
             fnr = overall.e.30@FNR,
             ppp = overall.e.30@PPP,
             npp = overall.e.30@NPP,
             mcr = overall.e.30@MCR,
             or = overall.e.30@OR,
             kappa = overall.e.30@kappa,
             boot = i,
             mode=response.j, 
             extent=30)
    
    brt.overall.eval <- rbind(brt.overall.eval, overall.eval.100, overall.eval.30)
    brt.overall.eval.df <- rbind(brt.overall.eval.df, overall.eval.df.100, overall.eval.df.30)

    print(paste0("******COMPLETED BOOTSTRAP ", i, " OF ", boot, " BOOTSTRAPS FOR SCALE ", scale.j))
  }
  
}

#7. Save results----

write.csv(brt.overall.eval, "OverallBRTEvaluation_PredictionExtentTest.csv", row.names = FALSE)
write.csv(brt.overall.eval.df, "OverallBRTEvaluationDataFrame_PredictionExtentTest.csv", row.names = FALSE)

#8. Visualize----
#AUC
brt.test.auc <- brt.overall.eval %>% 
  pivot_wider(id_cols=c(boot), names_from=extent, values_from=auc, names_prefix="auc")

ggplot(brt.test.auc) +
  geom_point(aes(x=auc100, y=auc30))

#CCR
brt.test.ccr <- brt.overall.eval.df %>% 
  group_by(boot, mode, extent) %>% 
  summarize(ccr=max(ccr),
            ccr=max(ccr)) %>% 
  ungroup() %>% 
  pivot_wider(id_cols=boot, names_from=extent, values_from=ccr, names_prefix="ccr")

ggplot(brt.test.ccr) +
  geom_point(aes(x=ccr100, y=ccr30))

#Kappa
brt.test.kappa <- brt.overall.eval.df %>% 
  group_by(boot, mode, extent) %>% 
  summarize(kappa=max(kappa),
            ccr=max(ccr)) %>% 
  ungroup() %>% 
  pivot_wider(id_cols=boot, names_from=extent, values_from=kappa, names_prefix="kappa")

ggplot(brt.test.kappa) +
  geom_point(aes(x=kappa100, y=kappa30))

#9. Correlation----
cor.test(~auc30 + auc100, data=brt.test.auc)

cor.test(~ccr30 + ccr100, data=brt.test.ccr)

cor.test(~kappa30 + kappa100, data=brt.test.kappa)
