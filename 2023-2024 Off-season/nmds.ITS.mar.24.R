library(tidyverse)
library(readxl)
library(RColorBrewer)
library(vegan)
library(glue)
library(ggtext)

metadata <- read_excel(path="mar.24.fungal.metadata.xlsx")%>%
  mutate(combo = glue("{location} {treatment}"))

nmds <- read_tsv(file="final.opti_mcc.braycurtis.0.03.lt.ave.mar.24.fungal.nmds.axes")

metadata_nmds <- inner_join(metadata, nmds, by=c('sample'='group'))

metadata_nmds$combo = factor(metadata_nmds$combo, levels = c("BF live", "Stiles live", "BF dead", "Stiles dead", "BF none", "Stiles none"), ordered = TRUE)
ggplot(metadata_nmds, aes(x=axis1, axis2, color=combo)) +
  stat_ellipse(aes(x=axis1, axis2, fill=combo),
               geom="polygon", level=0.75,
               alpha=0.2, show.legend=FALSE) +
  geom_point() +
  coord_fixed() +
  labs(x="NMDS Axis 1", y="NMDS Axis 2",
       title="Comparison of Fungal Communities<br>Off-season Bottom Farm vs. Stiles Farm") +
  scale_color_manual(name=NULL,
                     breaks=c("BF live", "Stiles live", "BF dead", "Stiles dead", "BF none", "Stiles none"),
                     values=c('blue', 'dodgerblue', 'red', 'pink', 'black', 'dimgrey'),
                     labels=c("Bottom Farm<br>live sclerotia", "Stiles Farm<br>live sclerotia", 
                              "Bottom Farm<br>heat-killed sclerotia", "Stiles Farm<br>heat-killed sclerotia", 
                              "Bottom Farm<br>bulk soil", "Stiles Farm<br>bulk soil")) +
  scale_fill_manual(name=NULL,
                    breaks=c("BF live", "Stiles live", "BF dead", "Stiles dead", "BF none", "Stiles none"),
                    values=c('blue', 'dodgerblue', 'red', 'pink', 'black', 'dimgrey'),
                    labels=c("Bottom Farm<br>live sclerotia", "Stiles Farm<br>live sclerotia", 
                             "Bottom Farm<br>heat-killed sclerotia", "Stiles Farm<br>heat-killed sclerotia", 
                             "Bottom Farm<br>bulk soil", "Stiles Farm<br>bulk soil")) +
  theme_classic() +
  theme(legend.key.size = unit(1, "cm"),
        legend.text = element_markdown(size=9.5),
        legend.position = "bottom",
        plot.title = element_markdown(hjust=0.5, size=12, vjust=2))

ggsave("nmds.ITS.mar.24.png", width=7, height=6)
