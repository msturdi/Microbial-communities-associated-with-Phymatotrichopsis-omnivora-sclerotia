library(tidyverse)
library(readxl)
library(RColorBrewer)
library(vegan)
library(glue)
library(ggtext)

metadata <- read_excel(path="mar.24.bact.metadata.xlsx")%>%
  mutate(combo = glue("{location} {treatment}"))

nmds <- read_tsv(file="final.opti_mcc.braycurtis.0.03.lt.ave.mar.24.bact.nmds.axes")

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
       title="Comparison of Bacterial Communities<br>Off-season Bottom Farm vs. Stiles Farm") +
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
        #legend.position = c(1, 0.55),
        #legend.background = element_rect(fill="NA", color="black"),
        #legend.margin = margin(t=-2, r=3, b=1, l=2),
        plot.title = element_markdown(hjust=0.5, size=12, vjust=2)) 

nmds +
  #"BF" vertical line
  geom_line(data=tibble(x=c(0.35, 0.35), y=c(0.1, 0.3)), 
            aes(x=x, y=y),
            linetype="solid",
            linewidth= 0.4,
            inherit.aes = FALSE) +
  #bottom "BF" horizontal line
  geom_line(data=tibble(x=c(0.3, 0.35), y=c(0.1, 0.1)), 
            aes(x=x, y=y),
            linewidth= 0.4,
            linetype="solid",
            inherit.aes = FALSE) +
  #top "BF" horizontal line
  geom_line(data=tibble(x=c(0.3, 0.35), y=c(0.3, 0.3)), 
            aes(x=x, y=y),
            linewidth= 0.4,
            linetype="solid",
            inherit.aes = FALSE) +
  #"stiles" vertical line
  geom_line(data=tibble(x=c(0.4, 0.4), y=c(-0.4, -0.1)), 
            aes(x=x, y=y),
            linetype="solid",
            linewidth= 0.4,
            inherit.aes = FALSE) +
  #bottom "stiles" horizontal line
  geom_line(data=tibble(x=c(0.35, 0.4), y=c(-0.4, -0.4)), 
            aes(x=x, y=y),
            linewidth= 0.4,
            linetype="solid",
            inherit.aes = FALSE) +
  #top "stiles" vertical line
  geom_line(data=tibble(x=c(0.35, 0.4), y=c(-0.1, -0.1)), 
            aes(x=x, y=y),
            linewidth= 0.4,
            linetype="solid",
            inherit.aes = FALSE) +
  #"BF" text
  geom_text(data=tibble(x=0.40, y=0.2), 
            aes(x=x, y=y, label="Bottom\nFarm"), size=3.5, 
            inherit.aes = FALSE) +
  #"Stiles" text
  geom_text(data=tibble(x=0.44, y=-0.25), 
            aes(x=x, y=y, label="Stiles\nFarm"), size=3.5, 
            inherit.aes = FALSE)

ggsave("nmds.16S.mar.24.png", width=7, height=6)
