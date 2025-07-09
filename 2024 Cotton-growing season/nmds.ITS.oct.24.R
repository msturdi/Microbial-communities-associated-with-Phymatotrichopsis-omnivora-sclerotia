library(tidyverse)
library(readxl)
library(vegan)
library(glue)
library(ggtext)

metadata <- read_excel(path="oct.24.ITS.metadata.xlsx")%>%
  mutate(combo = glue("{location} {treatment}"))

nmds <- read_tsv(file="final.opti_mcc.braycurtis.0.03.lt.ave.oct.24.fungal.nmds.axes")

metadata_nmds <- inner_join(metadata, nmds, by=c('sample'='group'))

ggplot(metadata_nmds, aes(x=axis1, axis2, color=combo)) +
  stat_ellipse(aes(x=axis1, axis2, fill=combo),
               geom="polygon", level=0.75,
               alpha=0.2, show.legend=FALSE) +
  geom_point() +
  coord_fixed() +
  labs(x="NMDS Axis 1", y="NMDS Axis 2",
       title="Comparison of Fungal Communities<br>Cotton-growing Season Bottom Farm vs. Stiles Farm") +
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

ggsave("nmds.ITS.oct.24.png", width=7, height=6)
