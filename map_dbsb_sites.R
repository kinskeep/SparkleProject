# NOTE: ggplot's geom_map is probably more suitable than geom_polygon
# Change to geom_map for cleaner code, and likely more options

setwd('C:/Users/Katie/OneDrive/Documents/Feder Spring 2018')
library(ggplot2)
library(maps)
library(maptools)
library(mapdata)
library(raster)
library(grid)
library(RColorBrewer)
source('MarysFiles/MarysFiles/base_map.R')
source('MarysFiles/MarysFiles/gps_coord.R')
source('MarysFiles/MarysFiles/scaleBar.R')

gps_dat <- read.csv('gpsdat.csv')
mendax <- gps_dat[6:8,]
deerb <- gps_dat[1:4,]

# make basic map with formatting. theme elements to get rid of titles, ticks, and gridlines
base <- ggplot(canada) +
  geom_polygon(aes(x = long, y = lat, group = group),
               fill = 'oldlace', color = 'black', size = 0.2) +
  plot_area +
  theme(plot.background = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank(), panel.background = element_blank())
base # display the map

# Meredith's code
base2<-
  base_map(37, 46,-91.75,-80.25,'white', 'gray30', 'gray95')
full_map<-base_map(23, 46,-95,-59,Canada = T, 'white', 'gray30', 'gray95')

# add in gps points
# uses column names and file with gps coordinates
siteMap <- full_map +
  geom_point(data=gps_dat, aes(x = lat, y = long), shape=20, size=7, color="turquoise2") + 
  geom_point(data=mendax, aes(x = lat, y = long), shape=20, size=7, color="purple") +
  geom_point(data=deerb, aes(x = lat, y = long), shape=20, size=7, color="orange") +
  scale_shape_manual(values = c(16)) +
  geom_text(data=gps_dat, aes(x = lat, y = long, label=gps_dat$site), size=3, hjust=c(-0.1), vjust=c(0.5), check_overlap = TRUE)
siteMap

#Sympatric sites only

symMap <- full_map +
  geom_point(data=deerb, aes(x = lat, y = long), shape=20, size=7, color="orange") +
  scale_shape_manual(values = c(16)) +
  geom_text(data=deerb, aes(x = lat, y = long, label=deerb$site), size=3, hjust=c(-0.1), vjust=c(0.5), check_overlap = TRUE)
symMap
# recommend exporting as scalable file, such as pdf, then converting to image file
# this makes it simpler to specify resolution (ppi) for publication
#ggsave('collection_sites.png')
