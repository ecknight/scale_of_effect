#title: Analysis of common nighthawk scale of effect using boosted regression trees
#author: Elly C. Knight

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

my.theme <- theme_classic() +
  theme(text=element_text(size=16, family="Arial"),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(margin=margin(10,0,0,0)),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        axis.line.x=element_line(linetype=1),
        axis.line.y=element_line(linetype=1))

#load("/Users/ellyknight/Documents/UoA/Projects/Projects/Scale/Analysis/preGit/CONILAPRModel2019.RData")

#######PREP#########

#1. Load data and remove ABMI----
det <- read.csv("CONI_detection_2015.csv") %>% 
  filter(project != "ABMI") %>% 
  dplyr::select(-X)

rec <- read.csv("CONI_recording_2015.csv") %>% 
  separate(ID, into=c("project", "cluster", "site", "station"), sep="-", remove=FALSE) %>% 
  filter(project != "ABMI") %>% 
  dplyr::select(-X, -call)

loc <- read.csv("qryAllDeployments_2.csv") %>% 
  dplyr::rename(ID=StationKey,
                project=ProjectID,
                cluster=Cluster,
                site=SITE,
                station=STATION) %>% 
  dplyr::select(ID, project, cluster, site, station, Latitude, Longitude) %>% 
  unique()

off <- read.csv("Offsets.csv") %>% 
  dplyr::rename(file=file.name)

#2. Wrangle data----

#summarize by recording
det.sites <- det %>% 
  group_by(file) %>% 
  summarize(peent = n(), boom = sum(boom, na.rm=FALSE)) %>% 
  right_join(rec, by="file") %>% 
  mutate(peent = ifelse(is.na(peent), 0, peent),
         boom = ifelse(is.na(boom), 0, boom)) %>% 
  group_by(ID, year) %>% 
  summarize(peent = sum(peent, na.rm=FALSE), boom = sum(boom, na.rm=FALSE)) %>% 
  ungroup() %>% 
  mutate(site.peent = ifelse(peent>0, 1, 0),
         site.boom = ifelse(boom >0, 1, 0),
         site.notboom = ifelse(site.boom==0 & site.peent==1, 1, 0))

dat <- det %>% 
  group_by(file) %>% 
  summarize(peent = n(), boom = sum(boom)) %>% 
  right_join(rec, by="file") %>% 
  mutate(peent = ifelse(is.na(peent), 0, peent),
         boom = ifelse(is.na(boom), 0, boom)) %>% 
  dplyr::select(-cluster, -site, -station) %>% 
  left_join(loc, by=c("ID", "project")) %>% 
  left_join(off, by=c("file", "ID")) %>% 
  mutate(pres.peent = ifelse(peent>0, 1, 0),
         pres.boom = ifelse(boom >0, 1, 0)) %>% 
  left_join(det.sites) %>% 
  mutate(site.peent = ifelse(is.na(site.peent), 0, site.peent),
         site.boom = ifelse(is.na(site.boom), 0, site.boom),
         site.notboom = ifelse(is.na(site.notboom), 0, site.notboom),
         pres.notboom = ifelse(peent>0 & site.notboom==1, 1, 0))

#filter down to just p > 0.99
dat.99 <- dat %>% 
  dplyr::filter(p >0.99)

table(dat$pres.peent)
table(dat.99$pres.peent)
table(dat$pres.boom)
table(dat.99$pres.boom)
table(dat$pres.notboom)
table(dat.99$pres.notboom)

#3. Add covariates----

#get locations we need covariates for
covs <- dat.99 %>% 
  dplyr::select(ID, Latitude, Longitude) %>% 
  unique()

covs.sf <- covs %>% 
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326) %>% 
  st_transform(crs=3402)

#get covariates
setwd("/Volumes/ECK004/GIS/Projects/Scale/3MovingWindow")
tifs <- data.frame(file = list.files(pattern="*00.tif")) %>% 
  separate(file, into=c("cov", "scale", "tif"), remove=FALSE) %>% 
  mutate(name=paste0(cov, "_", scale)) %>% 
  unique()

for(i in 1:nrow(tifs)){
  raster.i <- raster(as.character(tifs$file[i]))
  covs.i <- covs.sf %>% 
    raster::extract(x=raster.i) %>% 
    data.frame()
  colnames(covs.i) <- tifs$name[i]
  covs <- covs %>% 
    cbind(covs.i)
#  layers <- stack(layers, raster.i)
}

#filter out locations that are outside our raster
covs.xy <- st_coordinates(covs.sf) %>% 
  data.frame()

covs.filter <- covs %>% 
  filter_at(vars(-ID, -Latitude, -Longitude), all_vars(is.na(.))) %>% 
  dplyr::select(ID) %>% 
  mutate(filter=1) %>% 
  full_join(covs) %>% 
  cbind(covs.xy) %>% 
  dplyr::filter(is.na(filter)) %>% 
  dplyr::select(-filter) %>% 
  filter(!is.na(ag_12800))

#4. Set up grid sampling----
grid <- dgconstruct(area=1, metric=TRUE)

covs.grid <- covs.filter %>% 
  mutate(cell = dgGEO_to_SEQNUM(grid, Longitude, Latitude)$seqnum) %>% 
  arrange(cell) %>% 
  dplyr::select(ID, cell)

table(covs.grid$cell)
length(unique(covs.grid$cell))

#5. Put it all back together----
dat.99.covs <- covs.filter %>% 
  rename(easting=Longitude, westing=Latitude) %>% 
  inner_join(dat.99) %>% 
  left_join(covs.grid)

#6. Determine # of samples----

#count # of samples per location
cov.samples <- dat.99.covs %>% 
  group_by(ID) %>% 
  summarize(n=n()) %>% 
  dplyr::filter(n > 1)

#Set minimum samples (n=8)
samples <- min(cov.samples$n)

#Filter out sites that do not have 8 samples
dat.use <- dat.99.covs %>% 
  inner_join(cov.samples)

table(dat.use$pres.peent)
table(dat.use$pres.boom)
table(dat.use$pres.notboom)

write.csv(dat.use, "CONI_CleanDataForAnalysis.csv", row.names = FALSE)

#Reload data if starting from here####
dat.use <- read.csv("CONI_CleanDataForAnalysis.csv") %>% 
  dplyr::select(-starts_with("ag")) %>% 
  unique()

samples <- 8

#7. Check covariates for covariation----
covs <- dat.use %>% 
  dplyr::select(contains("_")) %>% 
  unique()
M <- cor(covs) 

covs.100 <- dat.use %>% 
  dplyr::select(ends_with("_100")) %>% 
  filter(!is.na(nutrient_100)) %>% 
  unique()
M <- cor(covs.100)
corrplot::corrplot(M, method="circle")
vifstep(covs.100, th=10)
covs.100 %>% 
  gather(key=covariate, value=value) %>% 
  ggplot() +
  geom_histogram(aes(value)) +
  facet_wrap(~covariate, scales="free")

covs.200 <- dat.use %>% 
  dplyr::select(ends_with("_200"),
                -wetland_200) %>% 
  filter(!is.na(nutrient_200)) %>% 
  unique()
M <- cor(covs.200)
corrplot::corrplot(M, method="circle")
vifstep(covs.200, th=10)
covs.200 %>% 
  gather(key=covariate, value=value) %>% 
  ggplot() +
  geom_histogram(aes(value)) +
  facet_wrap(~covariate, scales="free")

covs.400 <- dat.use %>% 
  dplyr::select(ends_with("_400"),
                -wetland_400) %>% 
  filter(!is.na(nutrient_400)) %>% 
  unique()
M <- cor(covs.400)
corrplot::corrplot(M, method="circle")
vifstep(covs.400, th=10)
covs.400 %>% 
  gather(key=covariate, value=value) %>% 
  ggplot() +
  geom_histogram(aes(value)) +
  facet_wrap(~covariate, scales="free")

covs.800 <- dat.use %>% 
  dplyr::select(ends_with("_800"),
                -wetland_800) %>% 
  filter(!is.na(nutrient_800)) %>% 
  unique()
M <- cor(covs.800)
corrplot::corrplot(M, method="circle")
vifstep(covs.800, th=10)
covs.800 %>% 
  gather(key=covariate, value=value) %>% 
  ggplot() +
  geom_histogram(aes(value)) +
  facet_wrap(~covariate, scales="free")

covs.1600 <- dat.use %>% 
  dplyr::select(ends_with("_1600"),
                -wetland_1600) %>% 
  filter(!is.na(nutrient_1600)) %>% 
  unique()
M <- cor(covs.1600)
corrplot::corrplot(M, method="circle")
vifstep(covs.1600, th=10)
covs.1600 %>% 
  gather(key=covariate, value=value) %>% 
  ggplot() +
  geom_histogram(aes(value)) +
  facet_wrap(~covariate, scales="free")

covs.3200 <- dat.use %>% 
  dplyr::select(ends_with("_3200"),
                -wetland_3200) %>% 
  filter(!is.na(nutrient_3200),
         !is.na(pine_3200)) %>% 
  unique()
M <- cor(covs.3200)
corrplot::corrplot(M, method="circle")
vifstep(covs.3200, th=10) #Pine becomes problematic here
vif(covs.3200)
covs.3200 %>% 
  gather(key=covariate, value=value) %>% 
  ggplot() +
  geom_histogram(aes(value)) +
  facet_wrap(~covariate, scales="free")

covs.6400 <- dat.use %>% 
  dplyr::select(ends_with("_6400")) %>% 
  filter(!is.na(nutrient_6400),
         !is.na(pine_6400)) %>% 
  unique()
M <- cor(covs.6400)
corrplot::corrplot(M, method="circle")
vifstep(covs.6400, th=10)
vif(covs.6400)
covs.6400 %>% 
  gather(key=covariate, value=value) %>% 
  ggplot() +
  geom_histogram(aes(value)) +
  facet_wrap(~covariate, scales="free")

covs.12800 <- dat.use %>% 
  dplyr::select(ends_with("_12800"),
                -wetland_12800) %>% 
  filter(!is.na(nutrient_12800),
         !is.na(pine_12800)) %>% 
  unique()
M <- cor(covs.12800)
corrplot::corrplot(M, method="circle")
vifstep(covs.12800, th=10)
vif(covs.12800)
covs.12800 %>% 
  gather(key=covariate, value=value) %>% 
  ggplot() +
  geom_histogram(aes(value)) +
  facet_wrap(~covariate, scales="free")


#8. Set learning rate----
set.seed(1234)

#Randomly select one location per grid cell
locs.i <- dat.use %>% 
  dplyr::select(ID, westing, easting, cell) %>% 
  unique() %>% 
  group_by(cell) %>% 
  sample_n(size = 1) %>% 
  ungroup()

#Randomly select 8 recordings per selected ID
recs.i <- dat.use %>% 
  dplyr::filter(ID %in% locs.i$ID) %>% 
  group_by(ID) %>% 
  sample_n(size = samples) %>% 
  ungroup()

#Summarize data by peents & booms
dat.i <- recs.i %>% 
  group_by(ID) %>% 
  summarize(pres.peent = as.integer(ifelse(sum(pres.peent) > 0, 1, 0)),
            pres.boom = ifelse(sum(pres.boom) > 0, 1, 0)) %>% 
  left_join(dat.use %>% 
              dplyr::select(contains("_"), ID) %>% 
              unique()) %>% 
  ungroup() %>% 
  data.frame()

table(dat.i$pres.peent)
table(dat.i$pres.boom)

#Try looking at covariance now
covs.12800 <- dat.i %>% 
  dplyr::select(ends_with("_12800")) %>% 
  filter(!is.na(nutrient_12800),
         !is.na(pine_12800))
M <- cor(covs.12800)
corrplot::corrplot(M, method="circle")
vif(covs.12800)
vifstep(covs.12800, th=10)
covs.12800 %>% 
  gather(key=covariate, value=value) %>% 
  ggplot() +
  geom_histogram(aes(value)) +
  facet_wrap(~covariate, scales="free")

#Set parameters for optimization
Tc<-c(2,3,4,5) ##Tree Complexities

Lr<-c(0.01, 0.005, 0.001, 0.0005, 0.0001)  ## Learning Rates

pa_compare <- data.frame()

for(i in 1:length(Tc)){
  for(j in 1:length(Lr)){
    
    ##Run Models
    a <- dismo::gbm.step(data=dat.i, 
                        gbm.x=6:17,
                        gbm.y=2,
                        family="bernoulli",
                        tree.complexity = Tc[i],
                        learning.rate = Lr[j],
                        bag.fraction = 0.75,
                        max.trees=10000) 
    
    ### Extract model performance info
    b<-data.frame(model=paste("Tc=",Tc[i],"_","Lr=",Lr[j],sep=""),AUC=a$cv.statistics$discrimination.mean, Dev=a$cv.statistics$deviance.mean)
    
    ##Add to dataframe
    pa_compare<-rbind(b,pa_compare)
    
    ##Reassign model with unique name
    assign(paste0("brt.pa.Tc",Tc[i],".lr",Lr[j]),a)
    
    ##Track Progress
    print(paste("************", "Completed", i*j, "of",length(Tc)*length(Lr),"************",sep=" " ))
    print(b)
  }
}


##Reorder so that lowest deviance is first row
pa_compare<-pa_compare[order(pa_compare$Dev),]
pa_compare

#Plot optimization

opt.list.1 <- list(brt.pa.Tc2.lr0.01, brt.pa.Tc3.lr0.01, brt.pa.Tc4.lr0.01, brt.pa.Tc5.lr0.01,
                   brt.pa.Tc2.lr0.005, brt.pa.Tc3.lr0.005, brt.pa.Tc4.lr0.005, brt.pa.Tc5.lr0.005,
                   brt.pa.Tc2.lr0.001, brt.pa.Tc3.lr0.001, brt.pa.Tc4.lr0.001, brt.pa.Tc5.lr0.001,
                   brt.pa.Tc2.lr0.0005, brt.pa.Tc3.lr0.0005, brt.pa.Tc4.lr0.0005, brt.pa.Tc5.lr0.0005,
                   brt.pa.Tc2.lr0.0001, brt.pa.Tc3.lr0.0001, brt.pa.Tc4.lr0.0001, brt.pa.Tc5.lr0.0001)

par(mfrow=c(5,4), oma=c(4,4,3,3), mar=c(0.5,0.5,0.5,0.5))

for(i in 1:length(opt.list.1)){
  
  a <- opt.list.1[[i]]
  
  y.bar <- min(a$cv.values) 
  y.min <- min(a$cv.values - a$cv.loss.ses)
  y.max <- max(a$cv.values + a$cv.loss.ses)
  
  plot(a$trees.fitted, a$cv.values, type = 'l', axes=FALSE, xlim=c(0,10000), ylim=c(0.9, 1.1))
  abline(h = y.bar, col = 3)
  
  box(col="grey40")
  
  ifelse(i %in% c(17:20),
         axis(1, col="grey40", at=c(0,2000,4000,6000,8000,10000)),
         axis(1, labels=FALSE, tick=FALSE, col="grey40", at=c(0,2000,4000,6000,8000,10000)))
  
  ifelse(i %in% c(1,5,9,13,17),
         axis(2, col="grey40", at=c(0.9, 1.0, 1.1)),
         axis(2, labels=FALSE, tick=FALSE, col="grey40", at=c(0.9, 1.0, 1.1)))
  
  lines(a$trees.fitted, a$cv.values + a$cv.loss.ses, lty=2)  
  lines(a$trees.fitted, a$cv.values - a$cv.loss.ses, lty=2) 
  
  target.trees <- a$trees.fitted[match(TRUE,a$cv.values == y.bar)]
  abline(v = target.trees, col=4)
  
}

mtext("Number of trees", 1, at="center",outer=TRUE, padj=3)
mtext("Predictive deviance", 2, at="center", outer=TRUE, padj=-3)
mtext("Tree complexity", 3, at="center", outer=TRUE, padj=-2)
mtext("2", 3, at=0.125, outer=TRUE, cex=0.75)
mtext("3", 3, at=0.375, outer=TRUE, cex=0.75)
mtext("4", 3, at=0.625, outer=TRUE, cex=0.75)
mtext("5", 3, at=0.875, outer=TRUE, cex=0.75)
mtext("Learning rate", 4, at="center", outer=TRUE, padj=2)
mtext("0.01", 4, at=0.9, outer=TRUE, cex=0.75)
mtext("0.005", 4, at=0.7, outer=TRUE, cex=0.75)
mtext("0.001", 4, at=0.5, outer=TRUE, cex=0.75)
mtext("0.0005", 4, at=0.3, outer=TRUE, cex=0.75)
mtext("0.0001", 4, at=0.1, outer=TRUE, cex=0.75)

par(mfrow=c(1,1))

#Use tc = 2 and lr = 0.001
#But really, let's go with tc = 3

#9. Single scale models----
set.seed(1234)

boot <- 100

scales <- c(100, 200, 400, 800, 1600, 3200, 6400, 12800)
layers <- c("conifer", "decid", "fire", "gravel", "harvest", "industry", "mixed", "moisture", "nutrient", "pine", "roads", "seismic", "water", "wells", "wetland")

brt.perf <- data.frame()
brt.covs <- data.frame()
brt.pdp <- data.frame()
brt.int <- data.frame()

for(j in 1:length(scales)){
  
  scale.j <- scales[j]
  layers.j <- c("ID", "cell", "JULIAN", "TOD", "pres.peent", "pres.boom", "pres.notboom", paste0(layers, "_", scale.j))
  dat.j <- dat.use %>% 
    dplyr::select(layers.j)
  
  for(i in 1:boot){
    
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
    
    #peent brt
    brt.peent.i <- dismo::gbm.step(data=dat.i, 
                                   gbm.x=3:17,
                                   gbm.y=18,
                                   family="bernoulli",
                                   tree.complexity = 3,
                                   learning.rate = 0.001,
                                   bag.fraction = 0.75,
                                   max.trees=10000,
                                   verbose=FALSE) 
    
    #saveRDS(brt.peent.i, "ScaleOfEffectBRT_Peent.rds")
    
    #boom brt
    brt.boom.i <- dismo::gbm.step(data=dat.i, 
                                  gbm.x=3:17,
                                  gbm.y=19,
                                  family="bernoulli",
                                  tree.complexity = 3,
                                  learning.rate = 0.001,
                                  bag.fraction = 0.75,
                                  max.trees=10000,
                                  verbose=FALSE) 
    
    #saveRDS(brt.boom.i, "ScaleOfEffectBRT_Boom.rds")
    
    #notboom brt
    brt.notboom.i <- dismo::gbm.step(data=dat.i, 
                                  gbm.x=3:17,
                                  gbm.y=20,
                                  family="bernoulli",
                                  tree.complexity = 3,
                                  learning.rate = 0.001,
                                  bag.fraction = 0.75,
                                  max.trees=10000,
                                  verbose=FALSE) 
    
    #Save out model performance
    brt.perf.peent <- data.frame(response="peent",
                                 boot=i,
                                 scale=scale.j,
                                 total.dev = brt.peent.i$self.statistics$mean.null,
                                 train.auc = brt.peent.i$self.statistics$discrimination,
                                 train.dev = brt.peent.i$self.statistics$mean.resid,
                                 test.auc = brt.peent.i$cv.statistics$discrimination.mean,
                                 test.dev = brt.peent.i$cv.statistics$deviance.mean,
                                 trees = brt.peent.i$n.trees,
                                 n = nrow(dat.i),
                                 np = sum(dat.i$sum.peent))
    
    brt.perf.boom <- data.frame(response="boom",
                                boot=i,
                                scale=scale.j,
                                total.dev = brt.boom.i$self.statistics$mean.null,
                                train.auc = brt.boom.i$self.statistics$discrimination,
                                train.dev = brt.boom.i$self.statistics$mean.resid,
                                test.auc = brt.boom.i$cv.statistics$discrimination.mean,
                                test.dev = brt.boom.i$cv.statistics$deviance.mean,
                                trees = brt.boom.i$n.trees,
                                n = nrow(dat.i),
                                np = sum(dat.i$sum.boom))
    
    brt.perf.notboom <- data.frame(response="notboom",
                                boot=i,
                                scale=scale.j,
                                total.dev = brt.notboom.i$self.statistics$mean.null,
                                train.auc = brt.notboom.i$self.statistics$discrimination,
                                train.dev = brt.notboom.i$self.statistics$mean.resid,
                                test.auc = brt.notboom.i$cv.statistics$discrimination.mean,
                                test.dev = brt.notboom.i$cv.statistics$deviance.mean,
                                trees = brt.notboom.i$n.trees,
                                n = nrow(dat.i),
                                np = sum(dat.i$sum.notboom))
    
    brt.perf <- rbind(brt.perf, brt.perf.peent, brt.perf.boom, brt.perf.notboom)
    
    #Save out model results
    brt.covs.peent <- data.frame(summary(brt.peent.i)) %>% 
      mutate(response="peent",
             boot=i)
    row.names(brt.covs.peent) <- c()
    
    brt.covs.boom <- data.frame(summary(brt.boom.i)) %>% 
      mutate(response="boom",
             boot=i)
    row.names(brt.covs.boom) <- c()
    
    brt.covs.notboom <- data.frame(summary(brt.notboom.i)) %>% 
      mutate(response="notboom",
             boot=i)
    row.names(brt.covs.notboom) <- c()
    
    brt.covs <- rbind(brt.covs, brt.covs.peent, brt.covs.boom, brt.covs.notboom)
    
    #Save out partial dependency predictions
    brt.pdp.peent <- data.frame()
    
    for(k in 1:length(layers)){
      response.matrix.k <- gbm::plot.gbm(brt.peent.i, i.var=k, return.grid = TRUE) %>% 
        mutate(var = brt.peent.i$gbm.call$predictor.names[k],
               boot = i,
               response="peent") 
      colnames(response.matrix.k) <- c("x", "y", "var", "boot", "response")
      brt.pdp.peent<- rbind(brt.pdp.peent, response.matrix.k)
    }
    
    brt.pdp.boom <- data.frame()
    
    for(k in 1:length(layers)){
      response.matrix.k <- gbm::plot.gbm(brt.boom.i, i.var=k, return.grid = TRUE) %>% 
        mutate(var = brt.boom.i$gbm.call$predictor.names[k],
               boot = i,
               response="boom") 
      colnames(response.matrix.k) <- c("x", "y", "var", "boot", "response")
      brt.pdp.boom<- rbind(brt.pdp.boom, response.matrix.k)
    }
    
    brt.pdp.notboom <- data.frame()
    
    for(k in 1:length(layers)){
      response.matrix.k <- gbm::plot.gbm(brt.notboom.i, i.var=k, return.grid = TRUE) %>% 
        mutate(var = brt.notboom.i$gbm.call$predictor.names[k],
               boot = i,
               response="notboom") 
      colnames(response.matrix.k) <- c("x", "y", "var", "boot", "response")
      brt.pdp.notboom<- rbind(brt.pdp.notboom, response.matrix.k)
    }
    
    brt.pdp <- rbind(brt.pdp, brt.pdp.peent, brt.pdp.boom, brt.pdp.notboom)
    
    #Save out interactions
    int.peent.i <- gbm.interactions(brt.peent.i)$rank.list %>% 
      mutate(response="peent",
             boot=i)
    
    int.boom.i <- gbm.interactions(brt.boom.i)$rank.list %>% 
      mutate(response="boom",
             boot=i)
    
    int.notboom.i <- gbm.interactions(brt.notboom.i)$rank.list %>% 
      mutate(response="notboom",
             boot=i)
    
    brt.int <- rbind(brt.int, int.peent.i, int.boom.i, int.notboom.i)
    
    
    print(paste0("******COMPLETED BOOTSTRAP ", i, " OF ", boot, " BOOTSTRAPS FOR SCALE ", scale.j, "******"))
  }
  
}

#Scale and center partial dependency predictions
brt.pdp.scale <- brt.pdp %>%
  mutate(y.log = 1/(1+exp(-y)),
         y.scale = scale(y.log, center=TRUE, scale=FALSE)) %>% 
  separate(var, into=c("variable", "scale"), remove=FALSE)

write.csv(brt.covs, "BRTCovariates.csv", row.names = FALSE)
write.csv(brt.perf, "BRTPerformance.csv", row.names = FALSE)
write.csv(brt.pdp.scale, "BRTPartialPredictions.csv", row.names = FALSE)
write.csv(brt.int, "BRTInteractions.csv", row.names = FALSE)

#10. Spatial predictions for overall scale of effect----


scales <- data.frame(response=c("boom", "peent", "peent"),
                     scale=c(200, 1600, 6400))
layers <- c("conifer", "decid", "fire", "gravel", "harvest", "industry", "mixed", "moisture", "nutrient", "pine", "roads", "seismic", "water", "wells", "wetland")

tifs <- data.frame(file = list.files("/Volumes/ECK004/GIS/Projects/Scale/3MovingWindow", pattern="*00.tif")) %>% 
  separate(file, into=c("cov", "scale", "tif"), remove=FALSE) %>% 
  mutate(var=paste0(cov, "_", scale)) %>% 
  dplyr::filter(scale %in% scales$scale) %>% 
  unique()

raster <- list()
for(i in 1:nrow(tifs)){
  raster.i <- raster(as.character(paste0("/Volumes/ECK004/GIS/Projects/Scale/3MovingWindow/", tifs$file[i])))
  resamp <- raster(raster.i)
  res(resamp) <- c(100, 100)
  raster.low.i <- raster::resample(x=raster.i, y=resamp, method="bilinear")
  names(raster.low.i) <- tifs$var[i]
  raster[i] <- raster.low.i
}

raster.stack <- stack(raster)
rm(raster)

set.seed(1234)

boot <- 100

#Start dataframes
brt.overall.perf <- data.frame()
brt.overall.covs <- data.frame()
brt.overall.pdp <- data.frame()
brt.overall.int <- data.frame()
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
    
    time <- Sys.time()
    
    #Save out model performance
    brt.overall.perf.i <- data.frame(response=response.j,
                                 boot=i,
                                 scale=scale.j,
                                 total.dev = brt.overall.i$self.statistics$mean.null,
                                 train.auc = brt.overall.i$self.statistics$discrimination,
                                 train.dev = brt.overall.i$self.statistics$mean.resid,
                                 test.auc = brt.overall.i$cv.statistics$discrimination.mean,
                                 test.dev = brt.overall.i$cv.statistics$deviance.mean,
                                 trees = brt.overall.i$n.trees,
                                 n = nrow(dat.i),
                                 np = sum(dat.i$sum.peent))
    
    brt.overall.perf <- rbind(brt.overall.perf, brt.overall.perf.i)
    
    #Save out model results
    brt.overall.covs.i <- data.frame(summary(brt.overall.i)) %>% 
      mutate(response=response.j,
             boot=i)
    row.names(brt.overall.covs.i) <- c()
    
    brt.overall.covs <- rbind(brt.overall.covs, brt.overall.covs.i)
    
    #Save out partial dependency predictions
    brt.overall.pdp.i <- data.frame()
    
    for(k in 1:length(layers)){
      response.matrix.k <- gbm::plot.gbm(brt.overall.i, i.var=k, return.grid = TRUE) %>% 
        mutate(var = brt.overall.i$gbm.call$predictor.names[k],
               boot = i,
               response=response.j) 
      colnames(response.matrix.k) <- c("x", "y", "var", "boot", "response")
      brt.overall.pdp <- rbind(brt.overall.pdp, response.matrix.k)
    }
    
    brt.overall.pdp <- rbind(brt.overall.pdp, brt.overall.pdp.i)
    
    #Save out interactions
    int.overall.i <- gbm.interactions(brt.overall.i)$rank.list %>% 
      mutate(response=response.j,
             boot=i)
    
    brt.overall.int <- rbind(brt.overall.int, int.overall.i)
    
    elapsed <- Sys.time() - time
    
    print(paste0("******COMPLETED INTERACTIONS ", i, " OF ", boot, " INTERACTIONS FOR SCALE ", scale.j, " in ", elapsed, " seconds******"))
    
    time <- Sys.time()
    
    #Predict
    brt.overall.pred <- dismo::predict(raster.stack, brt.overall.i, n.trees=brt.overall.i$gbm.call$best.trees, type="response")
    writeRaster(brt.overall.pred, paste0("/Volumes/ECK004/GIS/Projects/Scale/5Predictions/OverallBRTPredictions_", response.j, "_", scale.j, "_", i, ".tif"), format="GTiff", overwrite=TRUE)
    
    elapsed <- Sys.time() - time
    #1.3 hours
    
    print(paste0("******COMPLETED PREDICTION ", i, " OF ", boot, " PREDICTIONS FOR SCALE ", scale.j, " in ", elapsed, " minutes******"))
    
    time <- Sys.time()
    
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
    overall.suit <- data.frame(suitability = raster::extract(brt.overall.pred, dat.test[,c("X", "Y")]))
    overall.val <- cbind(dat.test, overall.suit)
    overall.e <- dismo::evaluate(subset(overall.val, sum.response==1)$suitability, subset(overall.val, sum.response==0)$suitability)
    overall.eval.i <- data.frame(boot=i, mode=response.j)
    overall.eval.i$np <- overall.e@np
    overall.eval.i$na <- overall.e@na
    overall.eval.i$auc <- overall.e@auc
    overall.eval.i$cor <- overall.e@cor
    overall.eval.i$cor <- overall.e@pcor
    overall.eval.i$odp <- mean(overall.e@ODP)
    overall.eval.i$scale <- scale.j
    overall.eval.df.i <- data.frame(t = overall.e@t) %>% 
      mutate(prev = overall.e@prevalence,
             ccr = overall.e@CCR,
             tpr = overall.e@TPR,
             tnr = overall.e@TNR,
             fpr = overall.e@FPR,
             fnr = overall.e@FNR,
             ppp = overall.e@PPP,
             npp = overall.e@NPP,
             mcr = overall.e@MCR,
             or = overall.e@OR,
             kappa = overall.e@kappa,
             boot = i,
             mode=response.j)

    brt.overall.eval <- rbind(brt.overall.eval, overall.eval.i)
    brt.overall.eval.df <- rbind(brt.overall.eval.df, overall.eval.df.i)
    
    elapsed <- Sys.time() - time
    
    print(paste0("******COMPLETED BOOTSTRAP ", i, " OF ", boot, " BOOTSTRAPS FOR SCALE ", scale.j, " in ", elapsed, " minutes******"))
  }
  
}

brt.overall.pdp.scale <- brt.overall.pdp %>%
  mutate(y.log = 1/(1+exp(-y)),
         y.scale = scale(y.log, center=TRUE, scale=FALSE)) %>% 
  separate(var, into=c("variable", "scale"), remove=FALSE) %>% 
  rbind(brt.overall.pdp.scale)

write.csv(brt.overall.covs, "OverallBRTCovariates.csv", row.names = FALSE)
write.csv(brt.overall.perf, "OverallBRTPerformance.csv", row.names = FALSE)
write.csv(brt.overall.pdp.scale, "OverallBRTPartialPredictions.csv", row.names = FALSE)
write.csv(brt.overall.int, "OverallBRTInteractions.csv", row.names = FALSE)
write.csv(brt.overall.eval, "OverallBRTEvaluation.csv", row.names = FALSE)
write.csv(brt.overall.eval.df, "OverallBRTEvaluationDataFrame.csv", row.names = FALSE)

brt.overall.covs <- read.csv("OverallBRTCovariates.csv")
brt.overall.perf <- read.csv("OverallBRTPerformance.csv")
brt.overall.pdp.scale <- read.csv("OverallBRTPartialPredictions.csv")
brt.overall.int <- read.csv("OverallBRTInteractions.csv")
brt.overall.eval <- read.csv("OverallBRTEvaluation.csv")
brt.overall.eval.df <- read.csv("OverallBRTEvaluationDataFrame.csv")

#11. Choose scale of effect for each predictor----
brt.covs <- read.csv("BRTCovariates.csv")
brt.perf <- read.csv("BRTPerformance.csv")

brt.perf.scale <- brt.perf %>% 
  mutate(test.dev.exp = (total.dev - test.dev)/total.dev,
         train.dev.exp = (total.dev - train.dev)/total.dev)

brt.covs.scale <- brt.covs %>%
  separate(var, into=c("variable", "scale"), remove=FALSE) %>% 
  mutate(scale=as.numeric(scale)) %>% 
  left_join(brt.perf.scale) %>% 
  mutate(cov.test.dev = rel.inf*test.dev.exp,
         cov.train.dev = rel.inf*train.dev.exp)

brt.covs.scale.sum <- brt.covs.scale %>% 
  group_by(variable, scale, response) %>% 
  summarize(cov.test.dev.mean = mean(cov.test.dev),
            cov.test.dev.sd = sd(cov.test.dev),
            cov.train.dev.mean = mean(cov.train.dev),
            cov.train.dev.sd = sd(cov.train.dev)) %>% 
  ungroup()

brt.covs.scale.select <- brt.covs.scale.sum %>% 
  group_by(variable, response) %>% 
  summarize(cov.test.dev.mean = max(cov.test.dev.mean)) %>% 
  left_join(brt.covs.scale.sum) %>% 
  ungroup()

#12. Multiscale model----
set.seed(1234)

boot <- 100

#Choose the layers needed
layers.peent <- brt.covs.scale.select %>% 
  dplyr::filter(response=="peent") %>% 
  mutate(var = paste0(variable,"_",scale)) %>% 
  dplyr::select(var) %>% 
  arrange(var)

layers.boom <- brt.covs.scale.select %>% 
  dplyr::filter(response=="boom") %>% 
  mutate(var = paste0(variable,"_",scale)) %>% 
  dplyr::select(var) %>% 
  arrange(var)

layers.k.peent <- c("sum.peent", layers.peent$var)

layers.k.boom <- c("sum.boom", layers.boom$var)

#get covariates
setwd("/Volumes/ECK004/GIS/Projects/Scale/3MovingWindow")
tifs <- data.frame(file = list.files(pattern="*00.tif")) %>% 
  separate(file, into=c("cov", "scale", "tif"), remove=FALSE) %>% 
  mutate(var=paste0(cov, "_", scale)) %>% 
  right_join(rbind(layers.boom, layers.peent)) %>% 
  unique()

raster <- list()
for(i in 1:nrow(tifs)){
  raster.i <- raster(as.character(tifs$file[i]))
  resamp <- raster(raster.i)
  res(resamp) <- c(100, 100)
  raster.low.i <- raster::resample(x=raster.i, y=resamp, method="bilinear")
  names(raster.low.i) <- tifs$var[i]
  raster[i] <- raster.low.i
}

raster.stack <- stack(raster)
rm(raster)

#Start dataframes
brt.best.perf <- data.frame()
brt.best.covs <- data.frame()
brt.best.pdp <- data.frame()
brt.best.int <- data.frame()
brt.best.pred.peent <- list()
brt.best.pred.boom <- list()
brt.best.eval <- data.frame()
brt.best.eval.df <- data.frame()


for(i in 1:boot){
  
  time <- Sys.time()
  
  #Randomly select one location per grid cell
  locs.i <- dat.use %>% 
    dplyr::select(ID, cell) %>% 
    unique() %>% 
    group_by(cell) %>% 
    sample_n(size = 1) %>% 
    ungroup()
  
  #Randomly select 8 recordings per selected ID
  recs.i <- dat.use %>% 
    dplyr::filter(ID %in% locs.i$ID) %>% 
    group_by(ID) %>% 
    sample_n(size = samples) %>% 
    ungroup()
  
  #Summarize data by peents & booms
  sum.i <- recs.i %>% 
    group_by(ID) %>% 
    summarize(sum.peent = ifelse(sum(pres.peent) > 0, 1, 0),
              sum.boom = ifelse(sum(pres.boom) > 0, 1, 0)) %>% 
    ungroup()
  
  #Join back to covariates
  dat.peent.i <- recs.i %>% 
    right_join(sum.i) %>% 
    data.frame() %>% 
    mutate(sum.peent = as.integer(sum.peent)) %>% 
    dplyr::select(layers.k.peent) %>% 
    unique()

  dat.boom.i <- recs.i %>% 
    right_join(sum.i) %>% 
    data.frame() %>% 
    mutate(sum.boom = as.integer(sum.boom)) %>% 
    dplyr::select(layers.k.boom) %>% 
    unique()
  
  #peent brt
  brt.peent.i <- dismo::gbm.step(data=dat.peent.i, 
                                 gbm.x=2:16,
                                 gbm.y=1,
                                 family="bernoulli",
                                 tree.complexity = 3,
                                 learning.rate = 0.001,
                                 bag.fraction = 0.75,
                                 max.trees=10000,
                                 verbose=FALSE) 
  
  #boom brt
  brt.boom.i <- dismo::gbm.step(data=dat.boom.i, 
                                gbm.x=2:16,
                                gbm.y=1,
                                family="bernoulli",
                                tree.complexity = 3,
                                learning.rate = 0.001,
                                bag.fraction = 0.75,
                                max.trees=10000,
                                verbose=FALSE) 
  
  #Save out model performance
  brt.perf.peent <- data.frame(response="peent",
                               boot=i,
                               total.dev = brt.peent.i$self.statistics$mean.null,
                               train.auc = brt.peent.i$self.statistics$discrimination,
                               train.dev = brt.peent.i$self.statistics$mean.resid,
                               test.auc = brt.peent.i$cv.statistics$discrimination.mean,
                               test.dev = brt.peent.i$cv.statistics$deviance.mean,
                               trees = brt.peent.i$n.trees,
                               n = nrow(dat.boom.i),
                               np = sum(dat.boom.i$sum.boom))
  
  brt.perf.boom <- data.frame(response="boom",
                              boot=i,
                              total.dev = brt.boom.i$self.statistics$mean.null,
                              train.auc = brt.boom.i$self.statistics$discrimination,
                              train.dev = brt.boom.i$self.statistics$mean.resid,
                              test.auc = brt.boom.i$cv.statistics$discrimination.mean,
                              test.dev = brt.boom.i$cv.statistics$deviance.mean,
                              trees = brt.boom.i$n.trees,
                              n = nrow(dat.peent.i),
                              np = sum(dat.peent.i$sum.boom))
  
  brt.best.perf <- rbind(brt.best.perf, brt.perf.peent, brt.perf.boom)
  
  #Save out model results
  brt.covs.peent <- data.frame(summary(brt.peent.i)) %>% 
    mutate(response="peent",
           boot=i)
  row.names(brt.covs.peent) <- c()
  
  brt.covs.boom <- data.frame(summary(brt.boom.i)) %>% 
    mutate(response="boom",
           boot=i)
  row.names(brt.covs.boom) <- c()
  
  brt.best.covs <- rbind(brt.best.covs, brt.covs.peent, brt.covs.boom)
  
  #Save out partial dependency predictions
  brt.pdp.peent <- data.frame()
  
  for(k in 1:nrow(layers.peent)){
    response.matrix.k <- gbm::plot.gbm(brt.peent.i, i.var=k, return.grid = TRUE) %>% 
      mutate(var = brt.peent.i$gbm.call$predictor.names[k],
             boot = i,
             response="peent") 
    colnames(response.matrix.k) <- c("x", "y", "var", "boot", "response")
    brt.pdp.peent<- rbind(brt.pdp.peent, response.matrix.k)
  }
  
  brt.pdp.boom <- data.frame()
  
  for(k in 1:nrow(layers.boom)){
    response.matrix.k <- gbm::plot.gbm(brt.boom.i, i.var=k, return.grid = TRUE) %>% 
      mutate(var = brt.boom.i$gbm.call$predictor.names[k],
             boot = i,
             response="boom") 
    colnames(response.matrix.k) <- c("x", "y", "var", "boot", "response")
    brt.pdp.boom<- rbind(brt.pdp.boom, response.matrix.k)
  }
  
  brt.best.pdp <- rbind(brt.best.pdp, brt.pdp.peent, brt.pdp.boom)
  
  #Save out interactions
  int.peent.i <- gbm.interactions(brt.peent.i)$rank.list %>% 
    mutate(response="peent",
           boot=i)
  
  int.boom.i <- gbm.interactions(brt.boom.i)$rank.list %>% 
    mutate(response="boom",
           boot=i)
  
  brt.best.int <- rbind(brt.best.int, int.peent.i, int.boom.i)
  
  #Predict
  brt.best.pred.peent <- dismo::predict(raster.stack, brt.peent.i, n.trees=brt.peent.i$gbm.call$best.trees, type="response")
  writeRaster(brt.best.pred.peent, paste0("/Volumes/ECK004/GIS/Projects/Scale/5Predictions/BestBRTPredictions_Peent_", i, ".tif"), format="GTiff", overwrite=TRUE)
  
  brt.best.pred.boom <- dismo::predict(raster.stack, brt.boom.i, n.trees=brt.boom.i$gbm.call$best.trees, type="response")
  writeRaster(brt.best.pred.boom, paste0("/Volumes/ECK004/GIS/Projects/Scale/5Predictions/BestBRTPredictions_Boom_", i, ".tif"), format="GTiff", overwrite=TRUE)
  
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
  
  dat.peent.test <- recs.test %>% 
    right_join(sum.test) %>% 
    data.frame() %>% 
    mutate(sum.peent = as.integer(sum.peent)) %>% 
    dplyr::select(ID, X, Y, sum.peent) %>% 
    unique()
  
  dat.boom.test <- recs.test %>% 
    right_join(sum.test) %>% 
    data.frame() %>% 
    mutate(sum.boom = as.integer(sum.boom)) %>% 
    dplyr::select(ID, X, Y, sum.boom) %>% 
    unique()
  
  #Evaluate
  peent.suit <- data.frame(suitability = raster::extract(brt.best.pred.peent, dat.peent.test[,c("X", "Y")]))
  peent.val <- cbind(dat.peent.test, peent.suit)
  peent.e <- dismo::evaluate(subset(peent.val, sum.peent==1)$suitability, subset(peent.val, sum.peent==0)$suitability)
  peent.eval <- data.frame(boot=i, mode="peent")
  peent.eval$np <- peent.e@np
  peent.eval$na <- peent.e@na
  peent.eval$auc <- peent.e@auc
  peent.eval$cor <- peent.e@cor
  peent.eval$cor <- peent.e@pcor
  peent.eval$odp <- mean(peent.e@ODP)
  peent.eval.df <- data.frame(t = peent.e@t) %>% 
    mutate(prev = peent.e@prevalence,
           ccr = peent.e@CCR,
           tpr = peent.e@TPR,
           tnr = peent.e@TNR,
           fpr = peent.e@FPR,
           fnr = peent.e@FNR,
           ppp = peent.e@PPP,
           npp = peent.e@NPP,
           mcr = peent.e@MCR,
           or = peent.e@OR,
           kappa = peent.e@kappa,
           boot = i,
           mode="peent")
  
  boom.suit <- data.frame(suitability = raster::extract(brt.best.pred.boom, dat.boom.test[,c("X", "Y")]))
  boom.val <- cbind(dat.boom.test, boom.suit)
  boom.e <- dismo::evaluate(subset(boom.val, sum.boom==1)$suitability, subset(boom.val, sum.boom==0)$suitability)
  boom.eval <- data.frame(boot=i, mode="boom")
  boom.eval$np <- boom.e@np
  boom.eval$na <- boom.e@na
  boom.eval$auc <- boom.e@auc
  boom.eval$cor <- boom.e@cor
  boom.eval$cor <- boom.e@pcor
  boom.eval$odp <- mean(boom.e@ODP)
  boom.eval.df <- data.frame(t = boom.e@t) %>% 
    mutate(prev = boom.e@prevalence,
           ccr = boom.e@CCR,
           tpr = boom.e@TPR,
           tnr = boom.e@TNR,
           fpr = boom.e@FPR,
           fnr = boom.e@FNR,
           ppp = boom.e@PPP,
           npp = boom.e@NPP,
           mcr = boom.e@MCR,
           or = boom.e@OR,
           kappa = boom.e@kappa,
           boot = i,
           mode="boom")
  
  brt.best.eval <- rbind(brt.best.eval, peent.eval, boom.eval)
  brt.best.eval.df <- rbind(brt.best.eval.df, peent.eval.df, boom.eval.df)
  
  elapsed <- Sys.time() - time
  
  print(paste0("Completed bootstrap ", i, " in ", elapsed, " minutes"))
  
}

brt.best.pdp.scale <- brt.best.pdp %>%
  mutate(y.log = 1/(1+exp(-y)),
         y.scale = scale(y.log, center=TRUE, scale=FALSE)) %>% 
  separate(var, into=c("variable", "scale"), remove=FALSE)

write.csv(brt.best.covs, "BestBRTCovariates.csv", row.names = FALSE)
write.csv(brt.best.perf, "BestBRTPerformance.csv", row.names = FALSE)
write.csv(brt.best.pdp.scale, "BestBRTPartialPredictions.csv", row.names = FALSE)
write.csv(brt.best.int, "BestBRTInteractions.csv", row.names = FALSE)
write.csv(brt.best.eval, "BestBRTEvaluation.csv", row.names = FALSE)
write.csv(brt.best.eval.df, "BestBRTEvaluationDataFrame.csv", row.names = FALSE)

brt.best.covs <- read.csv("BestBRTCovariates.csv")
brt.best.perf <- read.csv("BestBRTPerformance.csv")
brt.best.pdp.scale <- read.csv("BestBRTPartialPredictions.csv")
brt.best.int <- read.csv("BestBRTInteractions.csv")
brt.best.eval <- read.csv("BestBRTEvaluation.csv")
brt.best.eval.df <- read.csv("BestBRTEvaluationDataFrame.csv")


#13. Merge predictions----

setwd("/Volumes/ECK004/GIS/Projects/Scale/5Predictions/")
files <- data.frame(file=list.files(pattern="*.tif", recursive=TRUE)) %>% 
  separate(file, into=c("model", "response", "i", "tif"), remove=FALSE) %>% 
  dplyr::filter(!is.na(tif))

runs <- c("Peent", "Boom")

for(i in 1:length(runs)){
  response.i <- paste0(runs[i])
  files.i <- files %>% 
    dplyr::filter(response==response.i)
  layers.i <- list()
  
  for(j in 1:nrow(files.i)){
    file <- as.character(files.i$file[j])
    layers.i[[j]] <- raster(file, band=1)
  }
  
  layers.i <- stack(layers.i)
  
#  mean.i <- stackApply(layers.i, indices =  rep(1,nlayers(layers.i)), fun = "mean", na.rm = FALSE)
  mean.i <- calc(layers.i, fun=mean)

  writeRaster(mean.i, paste0("/Volumes/ECK004/GIS/Projects/Scale/6MeanPredictions/", runs[i], "meanpredictions.tif"), format="GTiff", overwrite=TRUE)
  print(paste0("Wrote mean raster ", runs[i]))
  
  rm(layers.i)
  
}

for(i in 1:length(runs)){
  response.i <- paste0(runs[i])
  files.i <- files %>% 
    dplyr::filter(response==response.i)
  layers.i <- list()
  
  for(j in 1:nrow(files.i)){
    file <- as.character(files.i$file[j])
    layers.i[[j]] <- raster(file, band=1)
  }
  
  layers.i <- stack(layers.i)
  
  #sd.i <- stackApply(layers.i, indices =  rep(1,nlayers(layers.i)), fun = "sd", na.rm = FALSE)
  sd.i <- calc(layers.i, fun=sd)
  
  writeRaster(sd.i, paste0("/Volumes/ECK004/GIS/Projects/Scale/6MeanPredictions/", runs[i], "sdpredictions.tif"), format="GTiff", overwrite=TRUE)
  print(paste0("Wrote sd raster ", runs[i]))
  
  rm(layers.i)
  
}

#14. Covariate effects----
#Multiscale----
brt.best.pdp.scale <- read.csv("BestBRTPartialPredictions.csv")
head(brt.best.pdp.scale)

ggplot(brt.best.pdp.scale) +
  geom_smooth(aes(x=x, y=y.log)) +
  facet_wrap(response~variable, scales="free")

ggplot(brt.best.pdp.scale) +
  geom_hex(aes(x=x, y=y.log)) +
  facet_wrap(response~variable, scales="free")

vars <- brt.best.pdp.scale %>% 
  dplyr::select(variable, scale, response) %>% 
  unique()

preds.gam <- data.frame()
sum.gam <- data.frame()
for(i in 1:nrow(vars)){
  
  brt.best.pdp.scale.i <- brt.best.pdp.scale %>% 
    filter(variable==vars$variable[i],
           response==vars$response[i])
  
  gam.i <- gam(y.log ~ s(x), data=brt.best.pdp.scale.i)
  
  pred.i <- data.frame(predict(gam.i, newdata=data.frame(x=seq(min(brt.best.pdp.scale.i$x), max(brt.best.pdp.scale.i$x), by=0.001)), se.fit=TRUE)) %>% 
    cbind(data.frame(x=seq(min(brt.best.pdp.scale.i$x), max(brt.best.pdp.scale.i$x), by=0.001))) %>% 
    mutate(upr = fit + (1.96*se.fit),
           lwr = fit - (1.96*se.fit),
           variable=vars$variable[i],
           scale=vars$scale[i],
           response=vars$response[i])
  
  preds.gam <- rbind(preds.gam, pred.i)
  
  sum.gam <- rbind(sum.gam, data.frame(dev.expl = summary(gam.i)$dev.expl,
                            edf = summary(gam.i)$edf,
                            r.sq = summary(gam.i)$r.sq,
                            residual.df = summary(gam.i)$residual.df,
                            np = summary(gam.i)$np,
                            variable = vars$variable[i],
                            scale = vars$scale[i],
                            response = vars$response[i]))
  
  print(paste0("Finished number ", i, " of ", nrow(vars), " iterations"))
  
}

ggplot(preds.gam) +
  geom_line(aes(x=x, y=fit, colour=response)) +
  geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, group=response), alpha=0.3)+
  facet_wrap(variable~scale, scales="free")

write.csv(preds.gam, "BestBRTGamPredictions.csv", row.names = FALSE)
write.csv(sum.gam, "BestBRTGamSummary.csv", row.names = FALSE)

#Single scale----
brt.overall.pdp.scale <- read.csv("BRTPartialPredictions.csv")

vars <- brt.overall.pdp.scale %>% 
  dplyr::select(variable, scale, response) %>% 
  unique()

preds.gam <- data.frame()
sum.gam <- data.frame()
for(i in 1:nrow(vars)){
  
  brt.overall.pdp.scale.i <- brt.overall.pdp.scale %>% 
    filter(variable==vars$variable[i],
           response==vars$response[i],
           scale==vars$scale[i])
  
  gam.i <- gam(y.log ~ s(x), data=brt.overall.pdp.scale.i)
  
  pred.i <- data.frame(predict(gam.i, newdata=data.frame(x=seq(min(brt.overall.pdp.scale.i$x), max(brt.overall.pdp.scale.i$x), by=0.001)), se.fit=TRUE)) %>% 
    cbind(data.frame(x=seq(min(brt.overall.pdp.scale.i$x), max(brt.overall.pdp.scale.i$x), by=0.001))) %>% 
    mutate(upr = fit + (1.96*se.fit),
           lwr = fit - (1.96*se.fit),
           variable=vars$variable[i],
           scale=vars$scale[i],
           response=vars$response[i])
  
  preds.gam <- rbind(preds.gam, pred.i)
  
  sum.gam <- rbind(sum.gam, data.frame(dev.expl = summary(gam.i)$dev.expl,
                                       edf = summary(gam.i)$edf,
                                       r.sq = summary(gam.i)$r.sq,
                                       residual.df = summary(gam.i)$residual.df,
                                       np = summary(gam.i)$np,
                                       variable = vars$variable[i],
                                       scale = vars$scale[i],
                                       response = vars$response[i]))
  
  print(paste0("Finished number ", i, " of ", nrow(vars), " iterations"))
  
}

write.csv(preds.gam, "AllBRTGamPredictions.csv", row.names = FALSE)

preds.gam.peent <- preds.gam %>% 
  dplyr::filter(response=="peent")

preds.gam.boom <- preds.gam %>% 
  dplyr::filter(response=="boom")

plot.gam <- ggplot(preds.gam) +
  geom_line(aes(x=x, y=fit, colour=response)) +
  geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, group=response), alpha=0.3)+
  facet_grid(scale~variable, scales="free")

ggsave(plot.gam, file="figures/PredictionsForPDTGRevisions.jpeg", height=15, width=30)

#15. Visualize####

#15a. Single scale models----

ggplot(brt.covs.scale) +
  geom_violin(aes(x=variable, y=rel.inf, colour=response)) +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_wrap(~scale)

ggplot(brt.covs.scale) +
  geom_smooth(aes(x=log(scale), y=rel.inf, colour=response)) +
  #  geom_point(aes(x=log(scale), y=rel.inf, colour=response)) +
  facet_wrap(~variable, scales="free_x") +
  labs(x="log of scale (km)", y="relative influence")

ggplot(brt.covs.scale) +
  geom_smooth(aes(x=log(scale), y=cov.test.dev, colour=response)) +
  geom_violin(aes(y=cov.test.dev, x=log(scale), group=factor(log(scale)))) +
  geom_boxplot(aes(y=cov.test.dev, x=log(scale), group=factor(log(scale)))) +
  facet_wrap(variable~response, scales="free") +
  labs(x="log of scale (km)", y="% test deviance explained")

ggsave("CovariateTestDeviance.jpeg", device="jpeg", width=12, height=12, dpi=300, units="in")

plot.perf.1 <- ggplot(brt.perf.scale) +
  geom_boxplot(aes(y=train.dev.exp, x=factor(scale), colour=response))

plot.perf.2 <- ggplot(brt.perf.scale) +
  geom_boxplot(aes(y=test.dev.exp, x=factor(scale), colour=response))

plot.perf.3 <- ggplot(brt.perf.scale) +
  geom_boxplot(aes(y=train.auc, x=factor(scale), colour=response))

plot.perf.4 <- ggplot(brt.perf.scale) +
  geom_boxplot(aes(y=test.auc, x=factor(scale), colour=response))

grid.arrange(plot.perf.1, plot.perf.2, plot.perf.3, plot.perf.4, ncol=2, nrow=2)

ggplot(brt.pdp.scale) +
  geom_smooth(aes(x=x, y=y.scale, colour=response)) +
  facet_grid(variable~scale, scales="free")

#15b. Multiscale model----

brt.best.covs.scale <- brt.best.covs %>%
  separate(var, into=c("variable", "scale"), remove=FALSE) %>% 
  mutate(scale=as.numeric(scale)) %>% 
  left_join(brt.perf.scale) %>% 
  mutate(cov.test.dev = rel.inf*test.dev.exp,
         cov.train.dev = rel.inf*train.dev.exp)

ggplot(brt.best.covs.scale) +
  geom_boxplot(aes(y=cov.test.dev, x=variable, colour=response))

ggplot(brt.best.pdp.scale) +
  geom_smooth(aes(x=x, y=y, colour=response)) +
  facet_wrap(response~variable, scales="free")

#15c. Model performance----

#BRT performance
brt.best.perf.scale <- brt.best.perf %>% 
  mutate(test.dev.exp = (total.dev - test.dev)/total.dev,
         train.dev.exp = (total.dev - train.dev)/total.dev) %>% 
  mutate(model="multi",
         scale=NA)

brt.overall.perf.scale <- brt.overall.perf %>% 
  mutate(test.dev.exp = (total.dev - test.dev)/total.dev,
         train.dev.exp = (total.dev - train.dev)/total.dev) %>% 
  mutate(model="single")

brt.perf.scale <- rbind(brt.best.perf.scale, brt.overall.perf.scale)

#Prediction evaluation
brt.eval <- rbind(brt.best.eval %>% 
                    mutate(scale=NA,
                           model="multi") %>% 
                    rename(response=mode),
                  brt.overall.eval %>% 
                    mutate(model="single"))

brt.performance <- brt.perf.scale %>% 
  left_join(brt.eval)

ggplot(brt.perf.scale) +
  geom_boxplot(aes(y=test.dev.exp, x=response, colour=model))

#16. Spatial autocorrelation----

#Two different approaches

#16a. Try moran's I from rasters using neighbourhood = extent size for each layer----
tifs <- data.frame(file = list.files("/Volumes/ECK004/GIS/Projects/Scale/3MovingWindow", pattern="*00.tif")) %>% 
  separate(file, into=c("cov", "scale", "tif"), remove=FALSE) %>% 
  mutate(var=paste0(cov, "_", scale)) %>% 
  unique() %>% 
  dplyr::filter(cov %in% c("pine", "harvest", "fire"),
                as.numeric(scale) <= 3200)

autocorr <- data.frame()
for(i in 1:nrow(tifs)){
  
  #Resample to 100m resolution
  raster.i <- raster(as.character(tifs$file[i]))
  resamp <- raster(raster.i)
  res(resamp) <- c(100, 100)
  raster.low.i <- raster::resample(x=raster.i, y=resamp, method="bilinear")
  names(raster.low.i) <- tifs$var[i]
  
  #Create matrix to define neighbours
  neighbour <- as.numeric(tifs$scale[i])/100
  dims <- neighbour*2+1
  ones <- (dims^2-1)/2
  m <- matrix(c(rep(1,ones),0,rep(1,ones)),dims)
  
  #Calculate Moran's I
  mor <- Moran(raster.low.i, m)
  
  autocorr <- data.frame(tifs[i,]) %>% 
    mutate(moran=mor) %>% 
    rbind(autocorr)
  
  write.csv(autocorr, "MoransIFromRasters.csv", row.names = FALSE)
  
  print(paste0("Finished ", i, " of ", nrow(tifs), " layers"))
  
}

autocorr <- read.csv("MoransIFromRasters.csv")

ggplot(autocorr, aes(x=as.numeric(scale), y=moran, colour=cov)) +
  geom_point() +
  geom_line()


#16b. Calculating Moran's I from dat.use----
#https://stats.idre.ucla.edu/r/faq/how-can-i-calculate-morans-i-in-r/

dat.sites <- dat.use[,c(1,142:143,4:123)] %>% 
  unique()

dists <- as.matrix(dist(cbind(dat.sites$Longitude, dat.sites$Latitude)))
dists.inv <- 1/dists
diag(dists.inv) <- 0

mlist <- apply(dat.sites[,c(4:123)], 2, function(x) Moran.I(x, dists.inv, na.rm=TRUE))

m <- rbindlist(mlist, idcol=colnames(dat.sites[,c(4:123)])) %>% 
  rename(layer = conifer_100) %>% 
  separate(layer, into=c("variable", "distance"), remove=FALSE) %>% 
  mutate(distance = as.numeric(distance))

write.csv(m, "MoransIFromPoints.csv", row.names = FALSE)
m <- read.csv("MoransIFromPoints.csv")

ggplot(m, aes(x=distance, y=observed, colour=variable)) +
  geom_point() +
  geom_line() +
  facet_wrap(~variable)


#save.image("/Users/ellyknight/Documents/UoA/Projects/Projects/Scale/Analysis/preGit/CONILAPRModel2019.RData")
