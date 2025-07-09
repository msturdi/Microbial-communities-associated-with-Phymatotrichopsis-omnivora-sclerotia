library(tidyverse)
library(readxl)
library(RColorBrewer)
library(vegan)

metadata <- read_excel(path="16S.2023.metadata.xlsx")

nmds <- read_tsv(file="final.opti_mcc.braycurtis.0.03.lt.ave.nmds.axes")

metadata_nmds <- inner_join(metadata, nmds, by=c('sample'='group'))

#NMDS projection with ellipses
ggplot(metadata_nmds, aes(x=axis1, axis2, color=trt, fill=trt)) +
  stat_ellipse(geom="polygon", level=0.75,
               alpha=0.2, show.legend=FALSE) +
  geom_point() +
  coord_fixed() +
  labs(x="NMDS Axis 1", y="NMDS Axis 2",
       title="16S rRNA Gene NMDS Projection") +
  scale_color_manual(name=NULL,
                     breaks=c("live","dead", "none"),
                     values=c('#2c7bb6', '#d7191c', 'black'),
                     labels=c("Live", "Dead", "None")) +
  scale_fill_manual(name=NULL,
                    breaks=c("live","dead", "none"),
                    values=c('dodgerblue', 'pink', 'lightgrey'),
                    labels=c("Live", "Dead", "None")) +
  theme_classic() +
  theme(legend.key.size = unit(0.4, "cm"),
        legend.position = c(1, 0.8),
        legend.background = element_rect(fill="NA", color="black"),
        legend.margin = margin(t=-2, r=3, b=3, l=3)) 

ggsave("nmds.16S.png", width=8, height=6)