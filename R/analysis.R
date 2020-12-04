# This is an analysis script
### An R script to develop a map of
### the highest valued eBird observations in space
### based on 'completeness' of a given cell
### all data were aggregated in GEE
### prior to being read into R

# packages
library(sf)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(tidyr)
library(ggcorrplot)
library(GGally)
library(stringr)
library(mgcv)
library(mgcViz)
library(broom)
library(maptools)
library(patchwork)
library(forcats)
library(readr)
library(ggforce)
library(plotly)

# read in summarized data from GEE
summarized_dat <- st_read("Data/hexagon_grid_summary_from_GEE/ebird_hex_stats.shp") %>%
  st_set_crs(4326)

# this has all grid cells in the world
# first will want to remove any grid cell which doesn't fall on land
# read in countries dataset
data(wrld_simpl)

countries <- st_as_sf(wrld_simpl) %>%
  dplyr::select(NAME) %>%
  mutate(col.id=1:nrow(.))
# intersect those grids on land
grids_on_land <- as.data.frame(st_intersects(summarized_dat, countries))

# get list of grids on land from above intersection
list_of_grids_on_land <- summarized_dat %>%
  mutate(row.id = 1:nrow(.)) %>%
  left_join(grids_on_land, by="row.id") %>%
  dplyr::filter(complete.cases(col.id)) %>%
  dplyr::select(-row.id, -col.id) %>%
  st_set_geometry(NULL) %>%
  .$POLYFID

# filter dataset to these ones
land_dat <- summarized_dat %>%
  dplyr::filter(POLYFID %in% list_of_grids_on_land) %>%
  st_set_crs(4326)

# plot this in space
# takes a while to plot
plasma_pal <- c(viridis::plasma(n = 8))

ggplot()+
  geom_sf(data=land_dat, aes(fill=avg_compl))+
  geom_sf(data=countries, fill="transparent", color="black")+
  theme_bw()+
  scale_fill_gradientn(colours=plasma_pal)+
  theme(axis.text=element_text(color="black"))+
  theme(panel.grid=element_line(color="transparent"))+
  labs(fill="Average \ncompleteness")


# see how many nas we are dealing with for each of the columns
land_dat %>%
  st_set_geometry(NULL) %>%
  skimr::skim(.)

# so because everything is scaled to 'completeness' of grids there are none missing for that
# but, there are probably some areas in antarctica etc which don't have GEE data for the majority
# of things because for each of the variables it is about 300 grids missing data
# this isn't that many
# and so I'll remove these here
# I'll also drop the geometry for now
# as this kind of slows down things
land_dat_complete <- land_dat %>%
  st_set_geometry(NULL) %>%
  dplyr::select(-access_m) %>%
  dplyr::filter(complete.cases(.)) %>%
  mutate(avg_compl=scales::rescale(avg_compl))


# Now let's make a 'conceptual figure' which will serve as Figure 1 for the manuscript
# showing the theory behind these main points
df <- expand.grid(x=-100:100, y=-100:100) %>%
  mutate(x=scales::rescale(x)) %>%
  mutate(y=scales::rescale(y))# dataframe for all combinations

labels <- data.frame(x=c(0.20, 0.20, 0.80, 0.80), 
                     y=c(0.80, 0.15, 0.80, 0.15),
                     labels=c("Biodiversity data are urgently needed",
                              "Biodiversity data less urgently needed",
                              "Necessary biodiversity data already present",
                              "Biodiversity data established for future risks"))

pal_color <- c(viridis::viridis(n = 8))

ggplot(df, aes(x, y, fill=y-x)) +      # map fill to the sum of x & y
  geom_tile(alpha = 0.7) +      # let the grid show through a bit
  scale_fill_gradientn(colours=pal_color, name = "Sampling priority:  ", 
                       breaks=c(-1,1), labels = c("low", "high"))+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  theme(panel.grid.major=element_blank())+
  theme(panel.grid.minor=element_blank())+
  theme(axis.ticks=element_blank())+
  ylab("Risk of habitat conversion")+
  xlab("Average completeness")+
  geom_segment(aes(x=0.5, y=0, xend=0.5, yend=1), size=1, color="gray30")+
  geom_segment(aes(x=0, y=0.5, xend=1, yend=0.5), size=1, color="gray30")+
  geom_segment(aes(x=0, y=0, xend=1, yend=1), size=1.5, color="firebrick3")+
  scale_x_continuous(expand = c(0,0), limits = c(0,1), breaks=c(0,1), labels=c(0,1))+
  scale_y_continuous(expand = c(0,0), limits = c(0,1), breaks=c(0,1), labels=c(0,1))+
  geom_label(data=labels, aes(x=x, y=y, label=labels), fill="white", size=3)+
  theme(legend.position="bottom")


# plot the relationship between completeness of a grid cell
# and the risk of the grid cell
# at the 'grid cell level
ggplot(land_dat_complete, aes(x=avg_compl, y=risk_50))+
  geom_point(color="navyblue")+
  #geom_density_2d()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  theme(panel.grid.major=element_line(color="gray80"))+
  ylab("Risk of conversion")+
  xlab("Average completeness")+
  geom_segment(aes(x=0, y=0, xend=1, yend=1), size=2.5, color="gold")

# This makes sense and shows the data

# get the distance from the line for each point
# function to calculate distance
dist2d <- function(id) {
  dat <- land_dat_complete %>%
    dplyr::filter(POLYFID == id)
  a <- c(dat$avg_compl, dat$risk_50)
  z <- c(dat$avg_compl, dat$risk_80)
  b = c(0, 0)
  c = c(1, 1)
  v1 <- b - c
  v2 <- a - b
  v3 <- z - b
  m1 <- cbind(v1,v2)
  m2 <- cbind(v1,v3)
  d <- data.frame(distance_50=abs(det(m1))/sqrt(sum(v1*v1)),
                  distance_80=abs(det(m2))/sqrt(sum(v1*v1)),
                  POLYFID=id)
  return(d)
} 

distance_from_line <- bind_rows(lapply(land_dat_complete$POLYFID, dist2d))

# join these 'distances'
# and scale these to be from 0 (lowest priority) to 1 (highest priority)
land_dat_complete <- land_dat_complete %>%
  left_join(., distance_from_line) %>%
  mutate(distance_50=case_when(
    avg_compl > risk_50 ~ distance_50*-1,
    avg_compl < risk_50 ~ distance_50*1,
    avg_compl == risk_50 ~ distance_50
  )) %>%
  mutate(distance_80=case_when(
    avg_compl > risk_80 ~ distance_80*-1,
    avg_compl < risk_80 ~ distance_80*1,
    avg_compl == risk_80 ~ distance_80
  )) %>%
  mutate(risk_50_value=scales::rescale(distance_50 + abs(min(distance_50)))) %>%
  mutate(risk_80_value=scales::rescale(distance_80 + abs(min(distance_80))))

# plot again, but with the 'value' for each point
risk_50 <- ggplot(land_dat_complete, aes(x=avg_compl, y=risk_50, color=risk_50_value))+
  geom_point()+
  scale_color_gradientn(colours=pal_color, breaks=c(0,1), labels=c(0, 1))+
  #geom_density_2d()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  theme(panel.grid.major=element_line(color="gray80"))+
  ylab("Risk of conversion")+
  xlab("Average completeness")+
  geom_segment(aes(x=0, y=0, xend=1, yend=1), size=1.5, color="firebrick3")+
  ggtitle("a)")+
  labs(color="Priority: ")+
  theme(legend.position="none")

risk_50

# now make this as a map, so plotting the value of points at a hexagonal cell level
# in space
map_grids_dat <- land_dat %>%
  dplyr::select(POLYFID) %>%
  left_join(., land_dat_complete)

map_50 <- ggplot()+
  geom_sf(data=map_grids_dat, aes(fill=risk_50_value), color="transparent")+
  geom_sf(data=countries, fill="transparent", color="black")+
  theme_bw()+
  scale_fill_gradientn(colours=pal_color, breaks=c(0,1), labels=c(0, 1))+
  theme(axis.text=element_text(color="black"))+
  theme(panel.grid=element_line(color="transparent"))+
  labs(fill="Sampling priority: ")+
  theme(legend.position="bottom")+
  ggtitle("b)")

map_50

# put these two plots together to make figure 2
risk_50 + map_50 + plot_layout(ncol=1)


# The majority of grid cells with 0 completeness?
land_dat_complete %>%
  dplyr::filter(avg_compl==0) %>% 
  nrow(.)/nrow(land_dat_complete)

land_dat_complete %>%
  dplyr::filter(avg_compl==0) %>%
  ggplot(., aes(x=risk_50))+
  geom_histogram(bins=50, color="black", fill="steelblue3")+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("Mean proportion of at-risk grid cells")+
  ylab("Number of hexagonal cells")

ggsave("Final figures/Figure_S2.png", width=7, height=6, units="in")

land_dat_complete %>%
  dplyr::filter(avg_compl==0) %>%
  dplyr::filter(risk_50 %in% c(0, 1)) %>%
  .$risk_50 %>%
  table()

713/1909
309/1909
