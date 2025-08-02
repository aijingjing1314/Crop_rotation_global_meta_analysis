

library(dplyr)
library(ggplot2)
library(maps)
library(tidyverse)
library(readxl)

world.dat<-map_data("world")

ggplot() +
  geom_polygon(data=world.dat,aes(x=long,y=lat,group=group),
               fill="#dedede")+
  theme_bw()+
  scale_y_continuous(expand = expansion(mult=c(0,0)))+
  scale_x_continuous(expand = expansion(add=c(0,0)))+
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  labs(x=NULL,y=NULL)-> world.map

world.map

### Map
df<-read.csv(file.choose())

world.map + 
  geom_point(data = df, 
             aes(x = Longitude, 
                 y = Latitude, 
                 color = Type, 
                 shape = Type, 
                 size = 12), 
             stroke = 0.5,  
             fill = NA) +  
  scale_color_manual(values = c("Type2"="#000000", "Type3"="#448DCD", "Type4"="#F7AF34", "Type5"="#5B854A", "Type6"="#684D93"),
                     name = "EWEs") +
  scale_shape_manual(values = c(21, 22, 23, 24, 25),  
                     name = "EWEs") +  
  theme(legend.position = c(0.1, 0.3))
# output  8*4






