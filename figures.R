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

#Wrangling
library(mefa4)
int2 <- read.csv("IntervalsForFigure.csv") 
xt <- Xtab(CONI.hit ~ file.name + interval, int2)
int2.1 <- int2[rowSums(xt)>0,] %>% 
  mutate(detection=1)
int2.0 <- int2[rowSums(xt)<1,] %>% 
  mutate(detection=0)
int3 <- rbind(int2.1, int2.0) %>% 
  dplyr::select(JULIAN, TOD, p, detection) %>% 
  unique() %>% 
  mutate(Date = as.Date(JULIAN, origin="2015-01-01"))

fls2 <- read.csv("FilesForFigure.csv") %>% 
  mutate(Date = as.Date(JULIAN, origin="2015-01-01"))
pred <- read.csv("PredictionsForFigure.csv") %>% 
  mutate(Date = as.Date(JULIAN, origin="2015-01-01"))

plot.2 <- ggplot() +
  geom_raster(aes(x=Date, y=TOD*24, fill=p), data=pred, alpha=0.7) +
  scale_fill_viridis_c(name="Probability\nof detection", direction=-1) +
  geom_point(aes(x=Date, y=TOD*24, colour=factor(detection)), alpha=0.5, data=subset(int3, detection==0)) +
  geom_point(aes(x=Date, y=TOD*24, colour=factor(detection)), data=subset(int3, detection==1)) +
  scale_colour_manual(values=c("grey70", "grey30"), name="Common\nnighthawk\ndetection", labels=c("Absent", "Present")) +
  geom_contour(aes(x=Date, y=TOD*24, z=p, lty="dashed"), data=pred, breaks=c(0.99), colour="black", size=1.2) +
  scale_linetype_manual(values=c("solid"), name="Threshold\nfor habitat\nmodels", labels=c("0.99")) +  xlab("Date") +
  ylab("Hour") +
  my.theme
plot.2

ggsave(plot=plot.2, filename="figures/Figure2.jpeg", device="jpeg", width=8, height=6, units="in")

#Summary of detections results----
dat.use <- read.csv("CONI_CleanDataForAnalysis.csv") %>% 
  dplyr::select(-starts_with("ag"))

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

#Mean times for analysis
dat.use %>% 
  dplyr::select()

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
                                           "conifer", 
                                           "nutrient",
                                           "harvest",
                                           "wetland",
                                           "seismic",
                                           "moisture",
                                           "wells",
                                           "decid",
                                           "mixed",
                                           "water",
                                           "industry",
                                           "roads",
                                           "gravel"),
                                  labels=c("pine",
                                           "wildfire",
                                           "conifer",
                                           "soil nutrients",
                                           "harvest",
                                           "wetland probability",
                                           "seismic lines",
                                           "soil moisture",
                                           "wellpads",
                                           "deciduous",
                                           "mixedwood",
                                           "open water",
                                           "industry",
                                           "roads",
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
  group_by(variable, response) %>% 
  summarize(max = max(mean),
            mean = mean(mean)) %>% 
  ungroup() %>% 
  arrange(-mean)

#3b. total % deviance plot----
plot.dev <- ggplot() + 
  geom_area(data=brt.covs.mean, aes(x=as.numeric(scale.fact), y=mean, fill=variable)) +
  scale_fill_viridis_d(name="Covariate", direction=-1) +
  geom_line(data=brt.perf.dev, aes(y=mean, x=as.numeric(scale.fact)), colour="grey30", alpha=0.5) +
  geom_errorbar(data=brt.perf.dev, aes(x=as.numeric(scale.fact), ymin=low, ymax=up, alpha=soe), colour="grey30", show.legend = FALSE) +
  geom_point(data=brt.perf.dev, aes(y=mean, x=as.numeric(scale.fact), shape=soe), colour="grey30", show.legend = FALSE) +
  scale_alpha_manual(values=c(0.5, 1)) +
  scale_shape_manual(values=c(1,19)) +
  facet_wrap(~response) +
  labs(x="Extent (km)", y="% deviance explained") +
  scale_x_continuous(breaks=c(1:8),
                     labels=c("0.1", "0.2", "0.4", "0.8", "1.6", "3.2", "6.4", "12.8")) +
  my.theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
plot.dev

#3c. Save----

ggsave(plot=plot.dev, filename="figures/Figure3.jpeg", device="jpeg", width=8, height=5, units="in")


#Figure 4. Per variable scale of effect-----
#Need to run figure 3 wrangling first
num <- 20
clrs <- viridis::viridis(num)
pts <- data.frame(x=c(1:num),
                  y=c(1:num),
                  clrs=clrs)
ggplot(pts, aes(x=x, y=x, colour=clrs)) +
  geom_point(size=10) +
  scale_colour_manual(values=clrs)

plot.cov <- ggplot(brt.covs.mean) +
  geom_line(aes(y=mean, x=as.numeric(scale.fact), colour=response), alpha=0.4) +
  geom_errorbar(aes(x=scale.fact, ymin=low, ymax=up, colour=response, alpha=soe), show.legend=FALSE) +
  geom_point(aes(y=mean, x=scale.fact, colour=response, shape=soe, alpha=soe), show.legend=FALSE) +
  facet_wrap(~variable, scales="free_y", ncol=3) +
  labs(x="Extent (km)", y="% deviance explained") +
  scale_x_discrete(labels=c("0.1", "0.2", "0.4", "0.8", "1.6", "3.2", "6.4", "12.8")) +
  scale_colour_manual(values=clrs[c(2,12)], name="") +
  scale_alpha_manual(values=c(0.6, 1)) +
  scale_shape_manual(values=c(1,19)) +
  my.theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = 'bottom')
plot.cov

ggsave(plot=plot.cov, filename="figures/Figure4.jpeg", device="jpeg", width=8, height=12, units="in")

#3f. Looking at relationship between scale of effect & % deviance
brt.covs.max.terr <- brt.covs.mean %>% 
  group_by(variable, response) %>% 
  arrange(-mean) %>% 
  mutate(n=row_number()) %>% 
  ungroup()

ggplot(brt.covs.max.terr) +
  geom_point(aes(y=mean, x=n, colour=soe)) +
  facet_wrap(variable~response, scales="free")

#Figure 5. Spatial predictions----
map.peent <- raster("/Volumes/SSD/GIS/Projects/Scale/6MeanPredictions/Peentmeanpredictions.tif") %>% 
  aggregate(fact=10)
names(map.peent) <- "p"
map.peent.df <- as.data.frame(map.peent, xy=TRUE) %>% 
  dplyr::filter(!is.na(p)) %>% 
  mutate(response="Home range",
         measure="Mean")

map.boom <- raster("/Volumes/SSD/GIS/Projects/Scale/6MeanPredictions/Boommeanpredictions.tif") %>% 
  aggregate(fact=10)
names(map.boom) <- "p"
map.boom.df <- as.data.frame(map.boom, xy=TRUE) %>% 
  dplyr::filter(!is.na(p)) %>% 
  mutate(response="Territory",
         measure="Mean")

sd.peent <- raster("/Volumes/SSD/GIS/Projects/Scale/6MeanPredictions/Peentsdpredictions.tif") %>% 
  aggregate(fact=10)
names(sd.peent) <- "p"
sd.peent.df <- as.data.frame(sd.peent, xy=TRUE) %>% 
  dplyr::filter(!is.na(p)) %>% 
  mutate(response="Home range",
         measure="Standard deviation")

sd.boom <- raster("/Volumes/SSD/GIS/Projects/Scale/6MeanPredictions/Boomsdpredictions.tif") %>% 
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

extent <- raster("/Volumes/SSD/GIS/Projects/Scale/6MeanPredictions/Boommeanpredictions.tif") %>% 
  aggregate(fact=10) %>% 
  reclassify(cbind(-Inf, Inf, 1)) %>% 
  rasterToPolygons(dissolve=TRUE) %>% 
  st_as_sf()

plot.mean <- ggplot() +
  geom_raster(data = mean, aes(x = x, y = y, fill=p), na.rm=TRUE) +
#  geom_sf(data=extent, fill=NA, colour="grey55", size=0.3) +
  scale_fill_viridis_c(name="Mean\nhabitat use\nprobability", option="C") +
  xlab("") +
  ylab("")+
  my.theme + 
  coord_sf(datum=NA) +
  facet_wrap(~response) +
  theme(legend.title=element_text(size=14))
#plot.mean

plot.sd <- ggplot() +
  geom_raster(data = sd, aes(x = x, y = y, fill=p), na.rm=TRUE) +
#  geom_sf(data=extent, fill=NA, colour="grey55", size=0.3) +
  scale_fill_viridis_c(name="Habitat use\nprobability\nstandard\ndeviation", option="C") +
  xlab("") +
  ylab("")+
  my.theme + 
  coord_sf(datum=NA) +
  facet_wrap(~response) + 
  theme(strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.title = element_text(size=14))
#plot.sd

brt.perf.eval <- read.csv("Multi&singlescaleModelEvaluationSummary.csv") %>%
  dplyr::filter(metric!="test.dev.exp") %>% 
  dplyr::filter(ID!="Sing scale-6400")

brt.perf.eval$metric <- base::factor(brt.perf.eval$metric, levels=c("test.auc", "eval.auc", "ccr", "kappa"),
                          labels=c("Cross validation\nROC AUC",
                                   "Spatial prediction\nROC AUC",
                                   "Spatial prediction\ncorrect classification rate",
                                   "Spatial prediction\nCohen's kappa"))

brt.perf.eval$ID <- base::factor(brt.perf.eval$ID, levels=c("Multiscale-Multiscale", "Single scale-200", "Single scale-1600"),
                                 labels=c("Multiscale", "Single scale (0.2 km)", "Single scale (1.6 km)"))

clrs <- viridis::plasma(10, direction=-1)

plot.perf <- ggplot(brt.perf.eval) +
  geom_boxplot(aes(x=metric, y=value, colour=model), position="dodge2") +
  facet_wrap(~response) +
  my.theme +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_text(size=14),
        legend.position = "bottom") +
  ylim(c(0,1)) +
  scale_colour_manual(values=clrs[c(3,9)], name="") +
  guides(colour=guide_legend(nrow=1))
plot.perf

ggsave(plot=grid.arrange(plot.mean, plot.sd, plot.perf,
                         nrow=3,
                         heights=c(4,4,6),
                         widths=c(6.5, 1.5),
                         layout_matrix = rbind(c(1,1),
                                               c(2,2),
                                               c(3,NA))),
       filename="figures/Figure5.jpeg", device="jpeg", width=8, height=14, units="in")

#Figure 6. Covariates----
#Main effects----
preds.gam <- read.csv("AllBRTGamPredictions.csv")

preds.gam$variable <- factor(preds.gam$variable,
                             levels=c("pine",
                                      "fire",
                                      "conifer", 
                                      "nutrient",
                                      "harvest",
                                      "wetland",
                                      "seismic",
                                      "moisture",
                                      "wells",
                                      "decid",
                                      "mixed",
                                      "water",
                                      "industry",
                                      "roads",
                                      "gravel"),
                             labels=c("pine",
                                      "wildfire",
                                      "conifer",
                                      "soil nutrients",
                                      "harvest",
                                      "wetland probability",
                                      "seismic lines",
                                      "soil moisture",
                                      "wellpads",
                                      "deciduous",
                                      "mixedwood",
                                      "open water",
                                      "industry",
                                      "roads",
                                      "gravel roads"))

preds.gam$response <- factor(preds.gam$response, levels=c("boom", "peent"),
                             labels=c("Territory", "Home range"))

clrs <- viridis::viridis(20)

labels = preds.gam %>% 
  mutate(label=paste0(scale/1000, " km")) %>% 
  dplyr::select(scale, variable, label) %>% 
  arrange(variable, scale) %>% 
  unique() %>% 
  group_by(scale) %>% 
  mutate(order = row_number()) %>% 
  ungroup()

labels$label <- factor(labels$label, levels=c("0.1 km", "0.2 km", "0.4 km", "0.8 km", "1.6 km", "3.2 km", "6.4 km", "12.8 km"))

labs <- labels$label
names(labs) <- labels$order

preds.gam.peent <- preds.gam %>% 
  dplyr::filter(response=="Home range") %>% 
  left_join(labels) %>% 
  arrange(order, scale) %>% 
  dplyr::filter(order %in% c(1:6))

preds.gam.boom<- preds.gam %>% 
  dplyr::filter(response=="Territory") %>% 
  left_join(labels) %>% 
  arrange(order, scale) %>% 
  dplyr::filter(order %in% c(1:6))

vars <- preds.gam.peent %>% 
  dplyr::select(order, variable) %>% 
  unique() %>% 
  arrange(order) %>% 
  rename(var=variable)

scaleFUN <- function(x) sprintf("%.2f", x)

#Peent----
soe.peent <- list(c(0.4, 0.4, 0.4, 1, 0.4, 0.4, 0.4, 0.4),
                  c(0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 1),
                  c(0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 1, 0.4),
                  c(1, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4),
                  c(0.4, 1, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4),
                  c(1, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4))

plot.gam.peent <- list()
for(i in 1:nrow(vars)){
  
  preds.gam.i <- preds.gam.peent %>% 
    dplyr::filter(variable==vars$var[i])
  
  plot.gam.peent.i <- ggplot(preds.gam.i) +
    geom_ribbon(aes(x=x, ymin=lwr, ymax=upr), alpha=0.3)+
    geom_line(aes(x=x, y=fit, colour=response, alpha=label),  size=1) +
    facet_wrap(~label, scales="free_x", labeller=labeller(order=labs), ncol=8) +
    labs(x="", y=vars$var[i]) + 
    scale_alpha_manual(values=soe.peent[[i]]) +
    scale_colour_manual(values=clrs[12], name="") +
    scale_y_continuous(labels=scaleFUN) +
    my.theme +
    theme(plot.margin=unit(c(3, 5, -10, 10), unit="pt"),
          text=element_text(size=10, family="Arial"),
          axis.text.x=element_text(size=8, angle = 90, hjust = 1),
          axis.text.y=element_text(size=8),
          legend.position="none")
  
  
  if(i>1){
    plot.gam.peent.i <- plot.gam.peent.i +
      theme(strip.background = element_blank(),
            strip.text.x = element_blank())
  }
  
  if(i==1){
    plot.gam.peent.i <- plot.gam.peent.i +
      ggtitle("Home range") +
      theme(plot.margin=unit(c(10, 5, -10, 10), unit="pt"),
            plot.title=element_text(size=12))
  }
  
  plot.gam.peent[[i]] <- plot.gam.peent.i
  
  
}

#plot.gam.peent[[1]]
#plot.gam.peent[[2]]

plot.gam.peent.all <- grid.arrange(plot.gam.peent[[1]],
                                   plot.gam.peent[[2]],
                                   plot.gam.peent[[3]],
                                   plot.gam.peent[[4]],
                                   plot.gam.peent[[5]],
                                   plot.gam.peent[[6]],
                                   nrow=6,
                                   heights=c(1, rep(0.8,5)),
                                   left="Marginal effect on home range habitat use")

ggsave(plot=plot.gam.peent.all, filename="figures/Figure6_Peent.jpeg", device="jpeg", width=10, height=12, units="in")

#Boom----
soe.boom <- list(c(0.4, 0.4, 0.4, 1, 0.4, 0.4, 0.4, 0.4),
                 c(0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 1),
                 c(0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 1, 0.4),
                 c(1, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4),
                 c(0.4, 1, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4),
                 c(1, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4))


plot.gam.boom <- list()
for(i in 1:nrow(vars)){
  
  preds.gam.i <- preds.gam.boom %>% 
    dplyr::filter(variable==vars$var[i])
  
  plot.gam.boom.i <- ggplot(preds.gam.i) +
    geom_ribbon(aes(x=x, ymin=lwr, ymax=upr), alpha=0.3)+
    geom_line(aes(x=x, y=fit, colour=response, alpha=label),  size=1, show.legend=FALSE) +
    facet_wrap(~label, scales="free_x", labeller=labeller(order=labs), ncol=8) +
    labs(x="", y="") + 
    scale_alpha_manual(values=soe.boom[[i]]) +
    scale_colour_manual(values=clrs[2], name="") +
    scale_y_continuous(labels=scaleFUN) +
    my.theme +
    theme(plot.margin=unit(c(3, 5, -10, 10), unit="pt"),
          text=element_text(size=10, family="Arial"),
          axis.text.x=element_text(size=8, angle = 90, hjust = 1),
          axis.text.y=element_text(size=8),
          legend.position = "none")
  
  if(i>1){
    plot.gam.boom.i <- plot.gam.boom.i +
      theme(strip.background = element_blank(),
            strip.text.x = element_blank())
  }
  
  if(i==1){
    plot.gam.boom.i <- plot.gam.boom.i +
      ggtitle("Territorial") +
      theme(plot.margin=unit(c(10, 5, -10, 10), unit="pt"),
            plot.title=element_text(size=12))
  }
  
  plot.gam.boom[[i]] <- plot.gam.boom.i
  
  
}

#plot.gam.boom[[1]]
#plot.gam.boom[[6]]

plot.gam.boom.all <- grid.arrange(plot.gam.boom[[1]],
                                   plot.gam.boom[[2]],
                                   plot.gam.boom[[3]],
                                   plot.gam.boom[[4]],
                                   plot.gam.boom[[5]],
                                   plot.gam.boom[[6]],
                                   nrow=6,
                                  heights=c(1, rep(0.8,5)))

ggsave(plot=plot.gam.boom.all, filename="figures/Figure6_Boom.jpeg", device="jpeg", width=10, height=12, units="in")

ggsave(grid.arrange(plot.gam.peent.all, plot.gam.boom.all, ncol=2), filename="figures/Figure6.jpeg", width=16, height=10)

#Figure 7. Interactions----
#Wrangle
brt.best.int <- read.csv("BestBRTInteractions.csv") %>% 
  mutate(model="Multiscale")
#brt.overall.int <- read.csv("OverallBRTInteractions.csv") %>% 
#  mutate(model="Single scale")
brt.int <- read.csv("BRTInteractions.csv") %>% 
  mutate(model="Single scale")

brt.int.sum <- rbind(brt.int, brt.best.int) %>% 
  group_by(var1.names, var2.names, response, model) %>% 
  summarize(int.mean = mean(int.size),
            int.sd = sd(int.size)) %>% 
  ungroup() %>% 
  arrange(-int.mean) %>% 
  separate(var1.names, into=c("var1", "scale1"), remove=FALSE) %>% 
  separate(var2.names, into=c("var2", "scale2"), remove=FALSE) %>% 
  mutate(scale = ifelse(model=="Multiscale", "Multiscale", scale1)) %>% 
  dplyr::filter(scale=="Multiscale" |
                  scale=="200" & response=="boom" |
                  scale %in% c("1600", "6400") & response=="peent")

brt.int.all <- data.frame(expand.grid(var1=unique(brt.int.sum$var1),
                           var2=unique(brt.int.sum$var2),
                           scale=unique(brt.int.sum$scale),
                           response=unique(brt.int.sum$response))) %>% 
  left_join(brt.int.sum) %>% 
  mutate(int.mean = ifelse(is.na(int.mean), 1, int.mean)) %>% 
  dplyr::filter(scale=="Multiscale" |
                  scale=="200" & response=="boom" |
                  scale %in% c("1600", "6400") & response=="peent")

brt.int$var1 <- factor(brt.int$var1,
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
brt.int$var2 <- factor(brt.int$var2,
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
brt.int.all$response <- factor(brt.int.all$response, levels=c("boom", "peent"), labels=c("Territory", "Home range"))

#Matrix plots
ggplot(brt.int.all) +
  geom_raster(aes(var1, var2, fill=int.mean)) +
#  scale_y_discrete(limits = rev(levels(brt.int$var2))) +
  facet_grid(scale~response, scales="free") +
  scale_fill_viridis_c(name = "Interaction\nstrength") +
  my.theme +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

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

#Performance summary----
brt.best.perf <- read.csv("BestBRTPerformance.csv") %>% 
  mutate(model="Multiscale",
         scale=NA)
brt.overall.perf <- read.csv("OverallBRTPerformance.csv") %>% 
  mutate(model="Single scale")

brt.perf <- rbind(brt.best.perf, brt.overall.perf)  %>% 
  mutate(test.dev.exp = (total.dev - test.dev)/total.dev) %>% 
  dplyr::select(model, scale, response, boot, test.dev.exp, test.auc) %>% 
  gather(key="metric", value="value", test.dev.exp:test.auc) %>% 
#  group_by(model, scale, response, metric) %>% 
#  summarize(mean = mean(value),
#            sd = sd(value)) %>% 
#  ungroup() %>% 
  mutate(scale=ifelse(model=="Multiscale", "Multiscale", scale)) %>% 
  dplyr::select(response, model, scale, boot, metric, value)
brt.perf

brt.best.eval <- read.csv("BestBRTEvaluation.csv") %>% 
  mutate(model="Multiscale",
         scale="Multiscale")
brt.overall.eval <- read.csv("OverallBRTEvaluation.csv") %>% 
  mutate(model="Single scale")

brt.eval <- rbind(brt.best.eval, brt.overall.eval) %>% 
#  group_by(mode, model, scale) %>% 
#  summarize(mean = mean(auc),
#            sd = sd(auc)) %>% 
#  ungroup() %>% 
  rename(response=mode, value=auc) %>% 
  mutate(metric="eval.auc") %>% 
  dplyr::select(response, model, scale, boot, metric, value)
brt.eval

brt.best.kappa <- read.csv("BestBRTEvaluationDataFrame.csv") %>% 
  mutate(model="Multiscale",
         scale="Multiscale")
brt.overall.kappa <- read.csv("OverallBRTEvaluationDataFrame.csv") %>% 
  mutate(model="Single scale") %>% 
  mutate(scale=case_when(mode=="boom" ~ 200,
                         mode=="peent" & row_number() < 36502 ~ 1600,
                         row_number() >= 36502 ~ 6400))

brt.kappa <- rbind(brt.best.kappa, brt.overall.kappa) %>% 
  group_by(mode, model, scale, boot) %>% 
  summarize(kappa=max(kappa),
            ccr=max(ccr)) %>% 
  ungroup() %>% 
  pivot_longer(names_to="metric",
               values_to="value",
               cols=c(kappa, ccr)) %>% 
#  group_by(metric, mode, model, scale) %>% 
#  summarize(mean = mean(value),
#            sd = sd(value)) %>% 
#  ungroup() %>% 
  rename(response=mode) %>% 
  dplyr::select(response, model, scale, boot, metric, value)
brt.kappa

brt.perf.eval <- rbind(brt.eval, brt.perf, brt.kappa) %>% 
  mutate(ID = paste0(model, "-", scale))

write.csv(brt.perf.eval, "Multi&singlescaleModelEvaluationSummary.csv")

brt.perf.eval %>% 
  group_by(metric, response, model, scale) %>% 
  summarize(mean=mean(value)) 

#Interrogate interactions----
brt.best.int <- read.csv("BestBRTInteractions.csv") %>% 
  mutate(model="Multiscale")
#brt.overall.int <- read.csv("OverallBRTInteractions.csv") %>% 
#  mutate(model="Single scale")
brt.int <- read.csv("BRTInteractions.csv") %>% 
  mutate(model="Single scale")

brt.int.all <- rbind(brt.int, brt.best.int) %>% 
  separate(var1.names, into=c("var1", "scale1"), remove=FALSE) %>% 
  separate(var2.names, into=c("var2", "scale2"), remove=FALSE) %>% 
  mutate(scale = ifelse(model=="Multiscale", "Multiscale", scale1)) %>% 
  dplyr::filter(scale=="Multiscale" |
                  scale=="200" & response=="boom" |
                  scale %in% c("1600") & response=="peent")

#Mean interactions between single scale and multiscale models
brt.int.all %>% 
  group_by(response, model) %>% 
  summarize(mean=mean(int.size),
            sd=sd(int.size))

#top mean interactions in top single scale and multiscale models
brt.int.sum <- rbind(brt.int, brt.best.int) %>% 
  group_by(var1.names, var2.names, response, model) %>% 
  summarize(int.mean = mean(int.size),
            int.sd = sd(int.size)) %>% 
  ungroup() %>% 
  arrange(-int.mean) %>% 
  separate(var1.names, into=c("var1", "scale1"), remove=FALSE) %>% 
  separate(var2.names, into=c("var2", "scale2"), remove=FALSE) %>% 
  mutate(scale = ifelse(model=="Multiscale", "Multiscale", scale1)) %>% 
  dplyr::filter(scale=="Multiscale" |
                  scale=="200" & response=="boom" |
                  scale %in% c("1600") & response=="peent") %>% 
  arrange(response, scale, -int.mean) %>% 
  group_by(scale, response) %>% 
  ungroup() %>% 
  dplyr::select(response, scale, var1, var2, int.mean, int.sd) %>% 
  group_by(scale, response) %>% 
  top_n(3, int.mean) %>% 
  ungroup()
View(brt.int.sum)

#top mean interactions across single scale models
brt.int.var <- brt.int %>% 
  separate(var1.names, into=c("var1", "scale1")) %>% 
  separate(var2.names, into=c("var2", "scale2")) %>% 
  group_by(var1, var2, scale1, response) %>% 
  summarize(intmean=mean(int.size)) %>% 
  ungroup() %>% 
  group_by(var1, var2) %>% 
  mutate(varmean=mean(intmean)) %>% 
  ungroup() %>% 
  arrange(-varmean, response, as.numeric(scale1))
View(brt.int.var)

ggplot(brt.int.var) +
  geom_point(aes(x=log(as.numeric(scale1)), y=intmean, colour=response)) +
#  geom_smooth(aes(x=as.numeric(scale1), y=intmean, colour=response)) +
  facet_grid(var1~var2)

ggsave("figures/Interactions.jpeg", width=16, height=16, units="in")

#Mean strength of interactions 
brt.int.sum <- rbind(brt.int, brt.best.int) %>% 
  separate(var1.names, into=c("var1", "scale1"), remove=FALSE) %>% 
  separate(var2.names, into=c("var2", "scale2"), remove=FALSE) %>% 
  mutate(scale = ifelse(model=="Multiscale", "Multiscale", scale1)) %>% 
  group_by(scale, response, model) %>% 
  summarize(int.mean = mean(int.size),
            int.sd = sd(int.size),
            int.max = max(int.size)) %>% 
  ungroup() %>% 
  arrange(-int.mean)
View(brt.int.sum)

ggplot(brt.int.sum) +
  geom_point(aes(x=scale, y=int.mean, colour=response), size=3) +
#  geom_point(aes(x=scale, y=int.max, colour=response), size=3) +
  geom_errorbar(aes(x=scale, ymin=int.mean-int.sd, ymax=int.mean+int.sd, colour=response)) +
  facet_grid(~response, scales="free") +
  my.theme

#Interactions and performance
int.perf.auc <- brt.int.sum %>% 
  left_join(brt.perf %>% 
              mutate(scale = as.character(scale))) %>% 
  dplyr::filter(metric=="test.auc")

int.perf.dev <- brt.int.sum %>% 
  left_join(brt.perf %>% 
              mutate(scale = as.character(scale))) %>% 
  dplyr::filter(metric=="test.dev.exp")

ggplot(int.perf.auc, aes(x=int.mean, y=mean, colour=response, shape=model)) +
  geom_point(size=3)

ggplot(int.perf.dev, aes(x=int.mean, y=mean, colour=response, shape=model)) +
  geom_point(size=3)

#All response shapes----
library(mgcv)
pred.overall <- read.csv("OverallBRTPartialPredictions.csv")
#pred.best <- read.csv("BestBRTPartialPredictions.csv")

vars <- pred.overall %>% 
  dplyr::select(variable, scale, response) %>% 
  unique()

preds.gam <- data.frame()
sum.gam <- data.frame()
for(i in 1:nrow(vars)){
  
  brt.best.pdp.scale.i <- pred.overall %>% 
    filter(variable==vars$variable[i],
           response==vars$response[i],
           scale==vars$scale[i])
  
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

preds.gam$variable <- factor(preds.gam$variable,
                             levels=c("pine",
                                      "fire",
                                      "conifer", 
                                      "nutrient",
                                      "harvest",
                                      "wetland",
                                      "seismic",
                                      "moisture",
                                      "wells",
                                      "decid",
                                      "mixed",
                                      "water",
                                      "industry",
                                      "roads",
                                      "gravel"),
                             labels=c("pine",
                                      "wildfire",
                                      "conifer",
                                      "soil nutrients",
                                      "harvest",
                                      "wetland probability",
                                      "seismic lines",
                                      "soil moisture",
                                      "wellpads",
                                      "deciduous",
                                      "mixedwood",
                                      "open water",
                                      "industry",
                                      "roads",
                                      "gravel roads"))

preds.gam$response <- factor(preds.gam$response, levels=c("boom", "peent"),
                             labels=c("Territory", "Home range"))

clrs <- viridis::viridis(20)

labels = preds.gam %>% 
  mutate(label=paste0(variable, " (", scale/1000, " km)"), 
         facet=paste0(variable, "-", response)) %>% 
  dplyr::select(variable, scale, response, label) %>% 
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
  facet_wrap(~order, scales="free", labeller=labeller(order=labs), ncol=8) +
  labs(x="", y="Marginal effect on habitat use") +
  scale_colour_manual(values=clrs[c(5,20)], name="") +
  scale_alpha_manual(values=c(0.4, 1)) +
  my.theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = 'bottom')
plot.gam

ggsave(plot=plot.gam, filename="figures/PredictionsForPDTGRevisions.jpeg", device="jpeg", width=24, height=18, units="in")

#Figure R1: Spatial prediction resolution test----
brt.overall.eval <- read.csv( "OverallBRTEvaluation_PredictionExtentTest.csv")
brt.overall.eval.df <- read.csv("OverallBRTEvaluationDataFrame_PredictionExtentTest.csv")

#AUC
brt.test.auc <- brt.overall.eval %>% 
  pivot_wider(id_cols=c(boot), names_from=extent, values_from=auc, names_prefix="value") %>% 
  mutate(metric="ROC AUC")

#CCR
brt.test.ccr <- brt.overall.eval.df %>% 
  group_by(boot, mode, extent) %>% 
  summarize(ccr=max(ccr),
            ccr=max(ccr)) %>% 
  ungroup() %>% 
  pivot_wider(id_cols=boot, names_from=extent, values_from=ccr, names_prefix="value") %>% 
  mutate(metric="Correct classification rate")

#Kappa
brt.test.kappa <- brt.overall.eval.df %>% 
  group_by(boot, mode, extent) %>% 
  summarize(kappa=max(kappa),
            ccr=max(ccr)) %>% 
  ungroup() %>% 
  pivot_wider(id_cols=boot, names_from=extent, values_from=kappa, names_prefix="value") %>% 
  mutate(metric="Cohen's kappa")

brt.test <- rbind(brt.test.auc, brt.test.ccr, brt.test.kappa)

ggplot(brt.test) +
  geom_point(aes(x=value100, y=value30)) +
  facet_wrap(~metric, scales="free") +
  xlab("Predictive performance metric for 100 m resolution raster") +
  ylab("Predictive performance metric for\n30 m resolution raster") +
  my.theme

ggsave(file="figures/FigureR1.jpeg", device="jpeg", width=8.5, height=5, units="in")
