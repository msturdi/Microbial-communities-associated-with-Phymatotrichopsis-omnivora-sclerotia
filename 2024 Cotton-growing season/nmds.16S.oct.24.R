library(tidyverse)
library(readxl)
library(RColorBrewer)
library(vegan)
library(glue)
library(ggtext)

metadata <- read_excel(path="oct.24.16S.metadata.xlsx")%>%
  mutate(combo = glue("{location} {treatment}"))

nmds <- read_tsv(file="final.opti_mcc.braycurtis.0.03.lt.ave.oct.24.bacterial.nmds.axes")

metadata_nmds <- inner_join(metadata, nmds, by=c('sample'='group'))

#NMDS projection with ellipses
#BF is at the top of the figure, Stiles is at the bottom
nmds <- ggplot(metadata_nmds, aes(x=axis1, axis2, color=treatment)) +
  stat_ellipse(aes(x=axis1, axis2, fill=combo),
               geom="polygon", level=0.75,
               alpha=0.2, show.legend=FALSE) +
  geom_point() +
  coord_fixed() +
  labs(x="NMDS Axis 1", y="NMDS Axis 2",
       title="Comparison of Bacterial Communities<br>Cotton-growing Season Bottom Farm vs. Stiles Farm") +
  scale_color_manual(name=NULL,
                     breaks=c("live","dead", "none"),
                     values=c('#2c7bb6', '#d7191c', 'black'),
                     labels=c("Live<br>sclerotia", "Heat-killed<br>sclerotia", "Bulk soil")) +
  scale_fill_manual(name=NULL,
                    breaks=c("BF live","BF dead", "BF none", "Stiles live", "Stiles dead", "Stiles none"),
                    values=c('dodgerblue', 'pink', 'lightgrey', 'dodgerblue', 'pink', 'lightgrey'),
                    labels=c("BF Live", "BF Dead", "BF Bulk", "Stiles Live", "Stiles Dead", "Stiles Bulk")) +
  theme_classic() +
  theme(legend.key.size = unit(1, "cm"),
        legend.text = element_markdown(size=9.5),
        legend.position = "bottom",
        plot.title = element_markdown(hjust=0.5, size=12, vjust=2)) 

nmds +
  #"BF" vertical line
  geom_line(data=tibble(x=c(0.4, 0.4), y=c(-0.4, 0.03)), 
            aes(x=x, y=y),
            linetype="solid",
            linewidth= 0.4,
            inherit.aes = FALSE) +
  #bottom "BF" horizontal line
  geom_line(data=tibble(x=c(0.35, 0.4), y=c(-0.4, -0.4)), 
            aes(x=x, y=y),
            linewidth= 0.4,
            linetype="solid",
            inherit.aes = FALSE) +
  #top "BF" horizontal line
  geom_line(data=tibble(x=c(0.35, 0.4), y=c(0.03, 0.03)), 
            aes(x=x, y=y),
            linewidth= 0.4,
            linetype="solid",
            inherit.aes = FALSE) +
  #"stiles" vertical line
  geom_line(data=tibble(x=c(0.18, 0.18), y=c(0.04, 0.33)), 
            aes(x=x, y=y),
            linetype="solid",
            linewidth= 0.4,
            inherit.aes = FALSE) +
  #bottom "stiles" horizontal line
  geom_line(data=tibble(x=c(0.13, 0.18), y=c(0.04, 0.04)), 
            aes(x=x, y=y),
            linewidth= 0.4,
            linetype="solid",
            inherit.aes = FALSE) +
  #top "stiles" vertical line
  geom_line(data=tibble(x=c(0.13, 0.18), y=c(0.33, 0.33)), 
            aes(x=x, y=y),
            linewidth= 0.4,
            linetype="solid",
            inherit.aes = FALSE) +
  #"BF" text
  geom_text(data=tibble(x=0.45, y=-0.2), 
            aes(x=x, y=y, label="Bottom\nFarm"), size=3.5, 
            inherit.aes = FALSE) +
  #"Stiles" text
  geom_text(data=tibble(x=0.23, y=0.2), 
            aes(x=x, y=y, label="Stiles\nFarm"), size=3.5, 
            inherit.aes = FALSE)

ggsave("nmds.16S.oct.24.png", width=7, height=6)
