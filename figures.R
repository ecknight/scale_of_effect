#title: Figures for analysis of common nighthawk scale of effect using boosted regression trees
#author: Elly C. Knight

library(tidyverse)
library(sf)
library(ggmap)
library(ggsn)
library(raster)
library(gridExtra)
library(rgeos)
library(maptools)
library(rgdal)
library(mapproj)

options(scipen=999)

my.theme <- theme_classic() +
  theme(text=element_text(size=16, family="Arial"),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(margin=margin(10,0,0,0)),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        axis.line.x=element_line(linetype=1),
        axis.line.y=element_line(linetype=1))

#######NEED TO SAVE OUT #PRES & ABS WITH EACH RUN########

#Figure 1. Study area----

#1a. North America----
dat.use <- read.csv("CONI_CleanDataForAnalysis.csv") %>% 
  dplyr::select(-starts_with("ag"))
arus <- dat.use %>% 
  dplyr::select(ID, westing, easting) %>% 
  unique() %>% 
  dplyr::rename(lat=westing, long=easting)

center <- arus %>% 
  summarize(long=mean(long),
            lat=mean(lat))

sa <- read_sf("/Volumes/ECK001/GIS/Projects/Scale/avie_dep_extent_sa.shp")

sa.4326 <- sa %>% 
  st_transform(crs=4326)

nam <- map_data("world", region=c("Canada", 
                                    "USA", 
                                    "Mexico")) %>% 
  dplyr::filter(!group%in%c(258:264))

plot.na <- ggplot() +
  geom_polygon(data=nam, aes(x=long, y=lat, group=group), colour = "gray85", fill = "gray75", size=0.3) +
#  geom_sf(data=sa.4326, aes(fill=Id), fill="gray75", colour="black", size=2) +
  geom_point(data=center, aes(x=long, y=lat), shape=23, colour="gray85", fill="gray35", size=4) +
  coord_sf(datum = NA) +
  coord_map(projection = "mercator") +
  xlim(c(-170, -55)) +
  ylim(c(14,72)) +
  my.theme +
  theme(plot.margin = unit(c(0,0,-1,-1), "cm"),
        panel.background=element_rect(fill = NULL, colour = "black"),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  xlab("") +
  ylab("")
#plot.na

grob.na <- ggplotGrob(plot.na)
plot(grob.na)

#1b. Study area----
register_google(key="")

map <- get_map(center, zoom=6, force=TRUE, maptype="satellite")
ggmap(map)

cities <- data.frame(city = c("Fort McMurray"),
                     lat = c(56.7267),
                     long = c(-111.379))

plot.sa <-  ggmap(map) +
  geom_sf(data=sa.4326, aes(fill=Id), fill="white", alpha=0.15, colour="black", size=2, inherit.aes = FALSE) +
  geom_point(aes(x = long, y = lat),
             shape=21,
             colour="grey85",
             fill="grey55",
             data = arus, 
             alpha = 1,
             size=2) +
  geom_point(aes(x = long, y = lat),
             shape=23,
             colour="grey85",
             fill="black",
             data = cities, 
             alpha = 1,
             size=5) +
  geom_text(data=cities, 
            aes(x = long, y = lat, label=city),
            colour="grey85",
            hjust=1.15,
            vjust=0.5) +
  my.theme +
  coord_sf(datum = NA) +
  xlab("") +
  ylab("") +
  xlim(c(st_bbox(sa.4326)$xmin-0.1, st_bbox(sa.4326)$xmax+0.1)) +
  ylim(c(st_bbox(sa.4326)$ymin, st_bbox(sa.4326)$ymax+0.1)) +
  theme(plot.margin = unit(c(0,0,0,0), "cm")) +
  north(symbol=4,
        scale=0.07,
        x.min=st_bbox(sa.4326)$xmin,
        x.max=st_bbox(sa.4326)$xmax+0.1,
        y.min=st_bbox(sa.4326)$ymin,
        y.max=st_bbox(sa.4326)$ymax+0.15) +
  blank() +
  ggsn::scalebar(dist = 50, dist_unit = "km", location="bottomleft",
           st.size=3, st.bottom=FALSE, st.color="grey85",
           transform = TRUE, model = "WGS84",
           x.min=st_bbox(sa.4326)$xmin+0.1,
           x.max=st_bbox(sa.4326)$xmax,
           y.min=st_bbox(sa.4326)$ymin+0.1,
           y.max=st_bbox(sa.4326)$ymax) +
  annotation_custom(grob=grob.na,
                    xmin=st_bbox(sa.4326)$xmax-1.2,
                    xmax=st_bbox(sa.4326)$xmax+0.3,
                    ymin=st_bbox(sa.4326)$ymin-0.55,
                    ymax=st_bbox(sa.4326)$ymin+1) 

ggsave(plot=plot.sa, filename="figures/StudyArea.jpeg", device="jpeg", width=6, height=8, units="in")

#1c. Example survey station----
dat.use <- read.csv("CONI_CleanDataForAnalysis.csv")
dat.example <- dat.use %>% 
  dplyr::select(ID, Latitude, Longitude,
                conifer_6400,
                decid_6400,
                fire_6400,
                gravel_6400,
                harvest_6400,
                industry_6400,
                mixed_6400,
                moisture_6400,
                nutrient_6400,
                pine_6400,
                roads_6400,
                seismic_6400,
                water_6400,
                wells_6400,
                wetland_6400) %>% 
  unique() %>% 
  mutate(total = wells_6400+gravel_6400) %>% 
  arrange(-total) %>% 
  top_n(1) %>% 
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326) %>% 
  st_transform(crs="+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

buf1 <- st_buffer(dat.example, 100)
buf2 <- st_buffer(dat.example, 200)
buf4 <- st_buffer(dat.example, 400)
buf8 <- st_buffer(dat.example, 800)
buf16 <- st_buffer(dat.example, 1600)
buf32 <- st_buffer(dat.example, 3200)
buf64 <- st_buffer(dat.example, 6400)
buf128 <- st_buffer(dat.example, 12800)

coords <- dat.example %>% 
  st_coordinates() %>% 
  data.frame()

pts <- matrix(c(coords$X-15000, coords$Y-15000,
              coords$X-15000, coords$Y+15000,
              coords$X+15000, coords$Y+15000,
              coords$X+15000, coords$Y-15000,
              coords$X-15000, coords$Y-15000),
              ncol =2, byrow = T)

box <- st_polygon(list(pts)) %>% 
  st_sfc() %>% 
  st_set_crs("+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs") %>% 
  st_sf()

setwd("/Volumes/ECK001/GIS/Projects/Scale/2Reclassified/")
files <- list.files(pattern="*.tif")
for(i in 2:length(files)){
  name.i <- str_sub(files[i], -100, -5)
  layer.i <- raster(files[i]) %>% 
    crop(box)
  names(layer.i) <- "p"
  layer.df <- as.data.frame(layer.i, xy=TRUE) %>% 
    dplyr::filter(!is.na(p))
  assign(name.i, layer.df)
  rm(layer.i)
}

plot.eg <- ggplot() +
  geom_raster(data = conifer, aes(x = x, y = y, fill=p), na.rm=TRUE) +
  scale_fill_gradient2(low="white", high="forestgreen") +
  geom_sf(data=buf128, fill=NA, colour="black", size=0.8) +
  geom_sf(data=buf64, fill=NA, colour="black", size=0.8) +
  geom_sf(data=buf32, fill=NA, colour="black", size=0.8) +
  geom_sf(data=buf16, fill=NA, colour="black", size=0.8) +
  geom_sf(data=buf8, fill=NA, colour="black", size=0.8) +
  geom_sf(data=buf4, fill=NA, colour="black", size=0.8) +
  geom_sf(data=buf2, fill=NA, colour="black", size=0.8) +
  geom_sf(data=buf1, fill=NA, colour="black", size=0.8) +
  geom_sf(data=box, fill=NA, colour="black", size=0.8) +
  xlab("") +
  ylab("")+
  my.theme + 
  coord_sf(datum=NA) +
  theme(legend.position = "none",
        plot.background = element_blank(),
        panel.background = element_blank()) +
  ggsn::scalebar(dist = 5, dist_unit = "km", location="bottomleft",
                 st.size=3, st.bottom=FALSE, st.color="black",
                 transform = FALSE,
                 x.min=min(conifer$x)+500,
                 x.max=max(conifer$x),
                 y.min=min(conifer$y)+500,
                 y.max=max(conifer$y))
grob.eg <- ggplotGrob(plot.eg)
plot.eg


#1d. Input layers----

plot.r1 <- ggplot() +
  geom_raster(data = decid, aes(x = x, y = y, fill=p), na.rm=TRUE) +
  scale_fill_gradient2(low="white", high="olivedrab3") +
  geom_sf(data=box, fill=NA, colour="black") +
  xlab("") +
  ylab("")+
  my.theme + 
  coord_sf(datum=NA) +
  theme(legend.position = "none",
        plot.background = element_blank(),
        panel.background = element_blank())
grob.r1 <- ggplotGrob(plot.r1)

plot.r2 <- ggplot() +
  geom_raster(data = mixed, aes(x = x, y = y, fill=p), na.rm=TRUE) +
  scale_fill_gradient2(low="white", high="darkolivegreen4") +
  geom_sf(data=box, fill=NA, colour="black") +
  xlab("") +
  ylab("")+
  my.theme + 
  coord_sf(datum=NA) +
  theme(legend.position = "none",
        plot.background = element_blank(),
        panel.background = element_blank())
grob.r2 <- ggplotGrob(plot.r2)

plot.r3 <- ggplot() +
  geom_raster(data = pine, aes(x = x, y = y, fill=p), na.rm=TRUE) +
  scale_fill_gradient2(low="white", high="limegreen") +
  geom_sf(data=box, fill=NA, colour="black") +
  xlab("") +
  ylab("")+
  my.theme + 
  coord_sf(datum=NA) +
  theme(legend.position = "none",
        plot.background = element_blank(),
        panel.background = element_blank())
grob.r3 <- ggplotGrob(plot.r3)

plot.r4 <- ggplot() +
  geom_raster(data = water, aes(x = x, y = y, fill=p), na.rm=TRUE) +
  scale_fill_gradient2(low="white", high="dodgerblue3") +
  geom_sf(data=box, fill=NA, colour="black") +
  xlab("") +
  ylab("")+
  my.theme + 
  coord_sf(datum=NA) +
  theme(legend.position = "none",
        plot.background = element_blank(),
        panel.background = element_blank())
grob.r4 <- ggplotGrob(plot.r4)

plot.r5 <- ggplot() +
  geom_raster(data = wetland, aes(x = x, y = y, fill=p), na.rm=TRUE) +
  scale_fill_gradient2(low="white", high="blue1") +
  geom_sf(data=box, fill=NA, colour="black") +
  xlab("") +
  ylab("")+
  my.theme + 
  coord_sf(datum=NA) +
  theme(legend.position = "none",
        plot.background = element_blank(),
        panel.background = element_blank())
grob.r5 <- ggplotGrob(plot.r5)

plot.r6 <- ggplot() +
  geom_raster(data = moisture, aes(x = x, y = y, fill=p), na.rm=TRUE) +
  scale_fill_gradient2(low="white", high="slateblue") +
  geom_sf(data=box, fill=NA, colour="black") +
  xlab("") +
  ylab("")+
  my.theme + 
  coord_sf(datum=NA) +
  theme(legend.position = "none",
        plot.background = element_blank(),
        panel.background = element_blank())
grob.r6 <- ggplotGrob(plot.r6)

plot.r7 <- ggplot() +
  geom_raster(data = nutrient, aes(x = x, y = y, fill=p), na.rm=TRUE) +
  scale_fill_gradient2(low="white", high="orange3") +
  geom_sf(data=box, fill=NA, colour="black") +
  xlab("") +
  ylab("")+
  my.theme + 
  coord_sf(datum=NA) +
  theme(legend.position = "none",
        plot.background = element_blank(),
        panel.background = element_blank())
grob.r7 <- ggplotGrob(plot.r7)

plot.r8 <- ggplot() +
  geom_raster(data = fire, aes(x = x, y = y, fill=p), na.rm=TRUE) +
  scale_fill_gradient2(low="white", high="firebrick4") +
  geom_sf(data=box, fill=NA, colour="black") +
  xlab("") +
  ylab("")+
  my.theme + 
  coord_sf(datum=NA) +
  theme(legend.position = "none",
        plot.background = element_blank(),
        panel.background = element_blank())
grob.r8 <- ggplotGrob(plot.r8)

plot.r9 <- ggplot() +
  geom_raster(data = harvest, aes(x = x, y = y, fill=p), na.rm=TRUE) +
  scale_fill_gradient2(low="white", high="coral4") +
  geom_sf(data=box, fill=NA, colour="black") +
  xlab("") +
  ylab("")+
  my.theme + 
  coord_sf(datum=NA) +
  theme(legend.position = "none",
        plot.background = element_blank(),
        panel.background = element_blank())
grob.r9 <- ggplotGrob(plot.r9)

plot.r10 <- ggplot() +
  geom_raster(data = wells, aes(x = x, y = y, fill=p), na.rm=TRUE) +
  scale_fill_gradient2(low="white", high="maroon4") +
  geom_sf(data=box, fill=NA, colour="black") +
  xlab("") +
  ylab("")+
  my.theme + 
  coord_sf(datum=NA) +
  theme(legend.position = "none",
        plot.background = element_blank(),
        panel.background = element_blank())
grob.r10 <- ggplotGrob(plot.r10)

plot.r11 <- ggplot() +
  geom_raster(data = industry, aes(x = x, y = y, fill=p), na.rm=TRUE) +
  scale_fill_gradient2(low="white", high="grey50") +
  geom_sf(data=box, fill=NA, colour="black") +
  xlab("") +
  ylab("")+
  my.theme + 
  coord_sf(datum=NA) +
  theme(legend.position = "none",
        plot.background = element_blank(),
        panel.background = element_blank())
grob.r11 <- ggplotGrob(plot.r11)

plot.r12 <- ggplot() +
  geom_raster(data = roads, aes(x = x, y = y, fill=p), na.rm=TRUE) +
  scale_fill_gradient2(low="white", high="grey10") +
  geom_sf(data=box, fill=NA, colour="black") +
  xlab("") +
  ylab("")+
  my.theme + 
  coord_sf(datum=NA) +
  theme(legend.position = "none",
        plot.background = element_blank(),
        panel.background = element_blank())
grob.r12 <- ggplotGrob(plot.r12)

plot.r13 <- ggplot() +
  geom_raster(data = gravel, aes(x = x, y = y, fill=p), na.rm=TRUE) +
  scale_fill_gradient2(low="white", high="grey30") +
  geom_sf(data=box, fill=NA, colour="black") +
  xlab("") +
  ylab("")+
  my.theme + 
  coord_sf(datum=NA) +
  theme(legend.position = "none",
        plot.background = element_blank(),
        panel.background = element_blank())
grob.r13 <- ggplotGrob(plot.r13)

plot.r14 <- ggplot() +
  geom_raster(data = seismic, aes(x = x, y = y, fill=p), na.rm=TRUE) +
  scale_fill_gradient2(low="white", high="chocolate4") +
  geom_sf(data=box, fill=NA, colour="black") +
  xlab("") +
  ylab("")+
  my.theme + 
  xlim(min(seismic$x)-42000, max(seismic$x)) +
  ylim(min(seismic$y)-45000, max(seismic$y)) +
  coord_sf(datum=NA) +
  theme(legend.position = "none")

plot.rs <- plot.r14 +
  annotation_custom(grob=grob.r13,
                    xmin=min(seismic$x)-10000,
                    xmax=max(seismic$x)-1000,
                    ymin=min(seismic$y)-10000,
                    ymax=max(seismic$y)-1000) +
  annotation_custom(grob=grob.r12,
                    xmin=min(seismic$x)-14000,
                    xmax=max(seismic$x)-4000,
                    ymin=min(seismic$y)-14000,
                    ymax=max(seismic$y)-4000) +
  annotation_custom(grob=grob.r11,
                    xmin=min(seismic$x)-17000,
                    xmax=max(seismic$x)-7000,
                    ymin=min(seismic$y)-17000,
                    ymax=max(seismic$y)-7000) +
  annotation_custom(grob=grob.r10,
                    xmin=min(seismic$x)-20000,
                    xmax=max(seismic$x)-10000,
                    ymin=min(seismic$y)-20000,
                    ymax=max(seismic$y)-10000) +
  annotation_custom(grob=grob.r9,
                    xmin=min(seismic$x)-23000,
                    xmax=max(seismic$x)-13000,
                    ymin=min(seismic$y)-23000,
                    ymax=max(seismic$y)-13000) +
  annotation_custom(grob=grob.r8,
                    xmin=min(seismic$x)-26000,
                    xmax=max(seismic$x)-16000,
                    ymin=min(seismic$y)-26000,
                    ymax=max(seismic$y)-16000) +
  annotation_custom(grob=grob.r7,
                    xmin=min(seismic$x)-29000,
                    xmax=max(seismic$x)-19000,
                    ymin=min(seismic$y)-29000,
                    ymax=max(seismic$y)-19000) +
  annotation_custom(grob=grob.r6,
                    xmin=min(seismic$x)-32000,
                    xmax=max(seismic$x)-22000,
                    ymin=min(seismic$y)-32000,
                    ymax=max(seismic$y)-22000) +
  annotation_custom(grob=grob.r5,
                    xmin=min(seismic$x)-35000,
                    xmax=max(seismic$x)-25000,
                    ymin=min(seismic$y)-35000,
                    ymax=max(seismic$y)-25000) +
  annotation_custom(grob=grob.r4,
                    xmin=min(seismic$x)-38000,
                    xmax=max(seismic$x)-28000,
                    ymin=min(seismic$y)-38000,
                    ymax=max(seismic$y)-28000) +
  annotation_custom(grob=grob.r3,
                    xmin=min(seismic$x)-41000,
                    xmax=max(seismic$x)-31000,
                    ymin=min(seismic$y)-41000,
                    ymax=max(seismic$y)-31000) +
  annotation_custom(grob=grob.r2,
                    xmin=min(seismic$x)-44000,
                    xmax=max(seismic$x)-34000,
                    ymin=min(seismic$y)-44000,
                    ymax=max(seismic$y)-34000) +
  annotation_custom(grob=grob.r1,
                    xmin=min(seismic$x)-47000,
                    xmax=max(seismic$x)-37000,
                    ymin=min(seismic$y)-47000,
                    ymax=max(seismic$y)-37000) +
  annotation_custom(grob=grob.eg,
                    xmin=min(seismic$x)-60000,
                    xmax=max(seismic$x)-40000,
                    ymin=min(seismic$y)-60000,
                    ymax=max(seismic$y)-40000)
plot.rs

#1e. Put together----
plot.1 <- grid.arrange(plot.sa, plot.rs,
             widths = c(3,4),
             heights = c(4),
             layout_matrix = rbind(c(1,2)))

ggsave(plot=plot.1, filename="figures/Figure1.jpeg", device="jpeg", width=12, height=6, units="in")

#Figure 2. Detection probability----



#Summary of detections results----
#Number of recordings
nrow(dat.use)

#Number of recordings with calls
table(dat.use$pres.peent)

#Number of recordings with booms
table(dat.use$pres.boom)

#Number of survey stations
det.stn <- dat.use %>% 
  group_by(ID) %>% 
  summarize(pres.peent = ifelse(sum(pres.peent) > 0, 1, 0),
            pres.boom = ifelse(sum(pres.boom) > 0, 1, 0)) %>% 
  ungroup()
nrow(det.stn)

#Number of stations with calls
table(det.stn$pres.peent)

#Number of stations with booms
table(det.stn$pres.boom)

#Mean number of stations with booms & calls after thinning
brt.perf <- read.csv("BRTPerformance.csv")
brt.perf.sum <- brt.perf %>% 
  group_by(response) %>% 
  summarize(n = mean(n),
            np = mean(np),
            perc = mean(np/(n)))
brt.perf.sum

#Figure 3. Scale of effect-----

#3a. Wrangle----
brt.perf <- read.csv("BRTPerformance.csv")

brt.perf$response <- factor(brt.perf$response,
                                  levels=c("boom", "peent"),
                                  labels=c("Territory", "Home range"))

brt.perf.scale <- brt.perf %>% 
  mutate(test.dev.exp = (total.dev - test.dev)/total.dev*100) %>% 
  dplyr::select(response, boot, scale, test.dev.exp, test.auc) %>% 
  gather(key="metric", value="value", test.dev.exp:test.auc)

brt.perf.scale$metric <- factor(brt.perf.scale$metric,
                                levels=c("test.auc", "test.dev.exp"),
                                labels=c("ROC AUC", "% deviance explained"))

brt.perf.max <- brt.perf.scale %>% 
  group_by(scale, response, metric) %>% 
  summarize(mean = mean(value)) %>% 
  group_by(response, metric) %>% 
  summarize(max = max(mean)) %>% 
  ungroup()

brt.perf.mean <- brt.perf.scale %>% 
  group_by(scale, response, metric) %>% 
  summarize(mean = mean(value),
            sd=sd(value)) %>% 
  ungroup() %>% 
  mutate(se=sd/10,
        up=mean+1.96*se,
         low=mean-1.96*se,
         scale.fact = factor(round(log(scale),2))) %>% 
  left_join(brt.perf.max) %>% 
  mutate(soe = ifelse(mean==max, "y", "n"))

brt.perf.dev <- brt.perf.mean %>% 
  dplyr::filter(metric=="% deviance explained")

brt.perf.auc <- brt.perf.mean %>% 
  dplyr::filter(metric=="ROC AUC")

brt.covs <- read.csv("BRTCovariates.csv")

brt.covs$response <- factor(brt.covs$response,
                            levels=c("boom", "peent"),
                            labels=c("Territory", "Home range"))

brt.covs.scale <- brt.covs %>%
  separate(var, into=c("variable", "scale"), remove=FALSE) %>% 
  mutate(scale=as.numeric(scale)) %>% 
  left_join(brt.perf) %>% 
  mutate(test.dev.exp = (total.dev - test.dev)/total.dev,
         cov.test.dev = rel.inf*test.dev.exp)

brt.covs.scale$variable <- factor(brt.covs.scale$variable,
                                  levels=c("pine",
                                           "fire",
                                           "harvest",
                                           "conifer", 
                                           "nutrient",
                                           "moisture",
                                           "seismic",
                                           "wetland",
                                           "wells",
                                           "decid",
                                           "mixed",
                                           "industry",
                                           "roads",
                                           "water",
                                           "gravel"),
                                  labels=c("pine",
                                           "wildfire",
                                           "harvest",
                                           "conifer",
                                           "soil nutrients",
                                           "soil moisture",
                                           "seismic lines",
                                           "wetland probability",
                                           "wellpads",
                                           "deciduous",
                                           "mixedwood",
                                           "industry",
                                           "roads",
                                           "open water",
                                           "gravel roads"))

brt.covs.max <- brt.covs.scale %>% 
  group_by(scale, response, variable) %>% 
  summarize(mean = mean(cov.test.dev)) %>% 
  group_by(response, variable) %>% 
  summarize(max = max(mean)) %>% 
  ungroup()

brt.covs.mean <- brt.covs.scale %>% 
  group_by(scale, response, variable) %>% 
  summarize(mean = mean(cov.test.dev),
            sd=sd(cov.test.dev)) %>% 
  ungroup() %>% 
  mutate(se=sd/10,
         up=mean+1.96*se,
         low=mean-1.96*se,
         scale.fact = factor(round(log(scale),2))) %>% 
  left_join(brt.covs.max) %>% 
  mutate(soe = ifelse(mean==max, "y", "n"))

brt.covs.mean2 <- brt.covs.mean %>% 
  group_by(variable) %>% 
  summarize(max = max(mean),
            mean = mean(mean)) %>% 
  ungroup() %>% 
  arrange(-mean)

#3b. AUC plot----
plot.auc <- ggplot(brt.perf.auc) +
  geom_line(aes(y=mean, x=as.numeric(scale.fact)), colour="grey30", alpha=0.5) +
  geom_errorbar(aes(x=scale.fact, ymin=low, ymax=up, alpha=soe), colour="grey30") +
  geom_point(aes(y=mean, x=scale.fact, alpha=soe), colour="grey30") +
  scale_alpha_manual(values=c(0.5, 1), name="Scale of\neffect", labels=c("No" ,"Yes")) +
  facet_wrap(~response) +
  labs(x="", y="ROC AUC") +
  ylim(c(0.7, 0.9)) +
  scale_x_discrete(labels=c("0.1", "0.2", "0.4", "0.8", "1.6", "3.2", "6.4", "12.8")) +
  my.theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
plot.auc

#3c. total % deviance plot----
plot.dev <- ggplot() + 
  geom_area(data=brt.covs.mean, aes(x=as.numeric(scale.fact), y=mean, fill=variable)) +
  scale_fill_viridis_d(name="Covariate", direction=-1) +
  geom_line(data=brt.perf.dev, aes(y=mean, x=as.numeric(scale.fact)), colour="grey30", alpha=0.5) +
  geom_errorbar(data=brt.perf.dev, aes(x=as.numeric(scale.fact), ymin=low, ymax=up, alpha=soe), colour="grey30", show.legend = FALSE) +
  geom_point(data=brt.perf.dev, aes(y=mean, x=as.numeric(scale.fact), alpha=soe), colour="grey30", show.legend = FALSE) +
  scale_alpha_manual(values=c(0.5, 1)) +
  facet_wrap(~response) +
  labs(x="Extent (km)", y="% deviance explained") +
  scale_x_continuous(breaks=c(1:8),
                     labels=c("0.1", "0.2", "0.4", "0.8", "1.6", "3.2", "6.4", "12.8")) +
  my.theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
plot.dev


#3d. Per variable scale of effect-----
clrs <- viridis::viridis(3)

plot.cov <- ggplot(brt.covs.mean) +
  geom_line(aes(y=mean, x=as.numeric(scale.fact), colour=response), alpha=0.4) +
  geom_errorbar(aes(x=scale.fact, ymin=low, ymax=up, colour=response, alpha=soe), show.legend=FALSE) +
  geom_point(aes(y=mean, x=scale.fact, colour=response, alpha=soe), show.legend=FALSE) +
  facet_wrap(~variable, scales="free_y") +
  labs(x="Extent (km)", y="% deviance explained") +
  scale_x_discrete(labels=c("0.1", "0.2", "0.4", "0.8", "1.6", "3.2", "6.4", "12.8")) +
  scale_colour_manual(values=clrs, name="") +
  scale_alpha_manual(values=c(0.4, 1)) +
  my.theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = 'bottom')
plot.cov

#3e. Put together----
plot.3 <- grid.arrange(plot.auc, plot.dev, plot.cov,
                       widths = c(4,1,4),
                       heights = c(2,2),
                       layout_matrix = rbind(c(1,NA,4),
                                             c(2,2,4)))

ggsave(plot=plot.3, filename="figures/Figure3.jpeg", device="jpeg", width=20, height=10, units="in")

#Figure 4. Spatial predictions----
map.peent <- raster("/Volumes/ECK001/GIS/Projects/Scale/6MeanPredictions/Peentmeanpredictions.tif") %>% 
  aggregate(fact=10)
names(map.peent) <- "p"
map.peent.df <- as.data.frame(map.peent, xy=TRUE) %>% 
  dplyr::filter(!is.na(p)) %>% 
  mutate(response="Home range",
         measure="Mean")

map.boom <- raster("/Volumes/ECK001/GIS/Projects/Scale/6MeanPredictions/Boommeanpredictions.tif") %>% 
  aggregate(fact=10)
names(map.boom) <- "p"
map.boom.df <- as.data.frame(map.boom, xy=TRUE) %>% 
  dplyr::filter(!is.na(p)) %>% 
  mutate(response="Territory",
         measure="Mean")

sd.peent <- raster("/Volumes/ECK001/GIS/Projects/Scale/6MeanPredictions/Peentsdpredictions.tif") %>% 
  aggregate(fact=10)
names(sd.peent) <- "p"
sd.peent.df <- as.data.frame(sd.peent, xy=TRUE) %>% 
  dplyr::filter(!is.na(p)) %>% 
  mutate(response="Home range",
         measure="Standard deviation")

sd.boom <- raster("/Volumes/ECK001/GIS/Projects/Scale/6MeanPredictions/Boomsdpredictions.tif") %>% 
  aggregate(fact=10)
names(sd.boom) <- "p"
sd.boom.df <- as.data.frame(sd.boom, xy=TRUE) %>% 
  dplyr::filter(!is.na(p)) %>% 
  mutate(response="Territory",
         measure="Standard deviation")

mean <- rbind(map.peent.df, map.boom.df)
mean$response <- factor(mean$response, levels=c("Territory", "Home range"))
sd <- rbind(sd.peent.df, sd.boom.df)
sd$response <- factor(sd$response, levels=c("Territory", "Home range"))
all <- rbind(mean, sd)

extent <- raster("/Volumes/ECK001/GIS/Projects/Scale/6MeanPredictions/Boommeanpredictions.tif") %>% 
  aggregate(fact=10) %>% 
  reclassify(cbind(-Inf, Inf, 1)) %>% 
  rasterToPolygons(dissolve=TRUE) %>% 
  st_as_sf()

plot.mean <- ggplot() +
  geom_raster(data = mean, aes(x = x, y = y, fill=p), na.rm=TRUE) +
  geom_sf(data=extent, fill=NA, colour="grey55", size=0.3) +
  scale_fill_viridis_c(name="Mean\nselection\nprobability") +
  xlab("") +
  ylab("")+
  my.theme + 
  coord_sf(datum=NA) +
  facet_wrap(~response)
#plot.mean

plot.sd <- ggplot() +
  geom_raster(data = sd, aes(x = x, y = y, fill=p), na.rm=TRUE) +
  geom_sf(data=extent, fill=NA, colour="grey55", size=0.3) +
  scale_fill_viridis_c(name="Selection\nprobability\nstandard\ndeviation") +
  xlab("") +
  ylab("")+
  my.theme + 
  coord_sf(datum=NA) +
  facet_wrap(~response) + 
  theme(strip.background = element_blank(),
    strip.text.x = element_blank())
#plot.sd

plot.4 <- grid.arrange(plot.mean, plot.sd,
                       nrow=2)

ggsave(plot=plot.4, filename="figures/Figure4.jpeg", device="jpeg", width=8, height=8, units="in")

#Summary of spatial predictions performance----
brt.best.perf <- read.csv("BestBRTPerformance.csv")

brt.best.perf.scale <- brt.best.perf %>% 
  mutate(test.dev.exp = (total.dev - test.dev)/total.dev) %>% 
  dplyr::select(response, boot, test.dev.exp, test.auc) %>% 
  gather(key="metric", value="value", test.dev.exp:test.auc) %>% 
  group_by(response, metric) %>% 
  summarize(mean = mean(value),
         sd = sd(value)) %>% 
  ungroup()
brt.best.perf.scale

brt.best.eval <- read.csv("BestBRTEvaluation.csv")
brt.best.eval.df <- read.csv("BestBRTEvaluationDataFrame.csv")

brt.best.mean <- brt.best.eval %>% 
  group_by(mode) %>% 
  summarize(auc.mean = mean(auc),
            auc.sd = sd(auc))
brt.best.mean

#Figure 5. Covariates----
preds.gam <- read.csv("BestBRTGamPredictions.csv") %>% 
  mutate(label=paste0(variable, " (", scale/1000, " km)"), 
         facet=paste0(variable, "-", response))
sum.gam <- read.csv("BestBRTGamSummary.csv")

preds.gam$variable <- factor(preds.gam$variable,
                                  levels=c("pine",
                                           "fire",
                                           "harvest",
                                           "conifer", 
                                           "nutrient",
                                           "moisture",
                                           "seismic",
                                           "wetland",
                                           "wells",
                                           "decid",
                                           "mixed",
                                           "industry",
                                           "roads",
                                           "water",
                                           "gravel"),
                                  labels=c("pine",
                                           "wildfire",
                                           "harvest",
                                           "conifer",
                                           "soil nutrients",
                                           "soil moisture",
                                           "seismic lines",
                                           "wetland probability",
                                           "wellpads",
                                           "deciduous",
                                           "mixedwood",
                                           "industry",
                                           "roads",
                                           "open water",
                                           "gravel roads"))
preds.gam$response <- factor(preds.gam$response, levels=c("boom", "peent"),
                             labels=c("Territory", "Home range"))

clrs <- viridis::viridis(3)

labels = preds.gam %>% 
  dplyr::select(variable, response, label) %>% 
  arrange(variable, response) %>% 
  unique() %>% 
  mutate(order = row_number())

labs <- labels$label
names(labs) <- labels$order

preds.gam.vars <- preds.gam %>% 
  left_join(labels)

plot.gam <- ggplot(preds.gam.vars) +
  geom_line(aes(x=x, y=fit, colour=response), size=1, alpha=0.8) +
  geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, group=response), alpha=0.3)+
  facet_wrap(~order, scales="free", labeller=labeller(order=labs), ncol=6) +
  labs(x="", y="Marginal effect on selection probability") +
  scale_colour_manual(values=clrs, name="") +
  scale_alpha_manual(values=c(0.4, 1)) +
  my.theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = 'bottom')
#plot.gam

ggsave(plot=plot.gam, filename="figures/Figure5.jpeg", device="jpeg", width=15, height=10, units="in")

#Summary of gam results----
brt.best.covs <- read.csv("BestBRTCovariates.csv")

brt.best.sum <- brt.best.covs %>% 
  separate(var, into=c("variable", "scale"), remove=FALSE) %>% 
  group_by(var, variable, scale, response) %>% 
  summarize(inf.mean = mean(rel.inf),
            inf.sd = sd(rel.inf)) %>% 
  ungroup() %>% 
  mutate(scale = as.numeric(scale)) %>% 
  left_join(sum.gam) %>% 
  arrange(response, -inf.mean) %>% 
  dplyr::select(response, variable, scale, inf.mean, inf.sd, edf, r.sq)

ggplot(brt.best.sum) +
  geom_point(aes(x=inf.mean, y=r.sq))

write.csv(brt.best.sum, "BestBRTSummary.csv", row.names = FALSE)

#Interrogate interactions----
brt.best.int <- read.csv("BestBRTInteractions.csv")
brt.int <- read.csv("BRTInteractions.csv")

brt.int.sum <- brt.int %>% 
  group_by(var1.names, var2.names, response) %>% 
  summarize(int.mean = mean(int.size),
            int.sd = sd(int.size)) %>% 
  ungroup() %>% 
  arrange(-int.mean) %>% 
  separate(var1.names, into=c("variable1", "scale1"), remove=FALSE) %>% 
  separate(var2.names, into=c("variable2", "scale2"), remove=FALSE) %>% 
  mutate(scale = scale1)

brt.int.all <- expand.grid(response=c("boom", "peent"),
                           scale=c("100", "200", "400", "800",
                                   "1600", "3200", "6400", "12800"),
                           variable1=c("pine",
                                       "fire",
                                       "harvest",
                                       "conifer", 
                                       "nutrient",
                                       "moisture",
                                       "seismic",
                                       "wetland",
                                       "wells",
                                       "decid",
                                       "mixed",
                                       "industry",
                                       "roads",
                                       "water",
                                       "gravel"),
                           variable2=c("pine",
                                       "fire",
                                       "harvest",
                                       "conifer", 
                                       "nutrient",
                                       "moisture",
                                       "seismic",
                                       "wetland",
                                       "wells",
                                       "decid",
                                       "mixed",
                                       "industry",
                                       "roads",
                                       "water",
                                       "gravel")) %>% 
  left_join(brt.int.sum) %>% 
  mutate(int.mean = ifelse(is.na(int.mean), 1, int.mean))

brt.best.int.sum <- brt.best.int %>% 
  group_by(var1.names, var2.names, response) %>% 
  summarize(int.mean = mean(int.size),
            int.sd = sd(int.size)) %>% 
  ungroup() %>% 
  arrange(-int.mean) %>% 
  separate(var1.names, into=c("variable1", "scale1"), remove=FALSE) %>% 
  separate(var2.names, into=c("variable2", "scale2"), remove=FALSE) %>% 
  mutate(scale = scale1)

plot.int <- ggplot(brt.int.sum) +
  geom_tile(aes(x=variable1, y=variable2, fill=int.mean)) +
  facet_grid(response~as.numeric(scale)) +
  my.theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_viridis_c()
plot.int

ggsave(plot=plot.int, filename="figures/Interactions.jpeg", device="jpeg", width=16, height=6, units="in")

plot.best.int <- ggplot(brt.best.int.sum) +
  geom_tile(aes(x=variable1, y=variable2, fill=int.mean)) +
  facet_wrap(~response) +
  my.theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_viridis_c()
plot.best.int

ggsave(plot=plot.best.int, filename="figures/BestBRTInteractions.jpeg", device="jpeg", width=6, height=4, units="in")

brt.int.scale <- brt.int.sum %>% 
  group_by(response, scale) %>% 
  summarize(max = max(int.mean),
            mean = mean(int.mean)) %>% 
  ungroup() %>% 
  mutate(scale = as.numeric(scale))

ggplot(brt.int.scale) +
  geom_point(aes(x=scale, y=mean, colour=response))

ggplot(brt.int.scale) +
  geom_point(aes(x=scale, y=max, colour=response))
