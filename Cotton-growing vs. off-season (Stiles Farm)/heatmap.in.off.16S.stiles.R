library(tidyverse)
library(readxl)
library(ggtext)
library(glue)
library(dplyr)
library(scales)

metadata <- read_excel("metadata.16S.in.off.stiles.xlsx") %>%
  mutate(combo = glue("{season} {treatment}"))

otu_counts <- read_tsv("final.opti_mcc.0.03.subsample.bacterial.stiles.shared") %>%
  select(Group, starts_with("Otu")) %>%
  rename(sample = Group) %>%
  pivot_longer(-sample, names_to="otu", values_to = "count")

taxonomy <- read_tsv(file="final.opti_mcc.0.03.cons.bacterial.stiles.taxonomy") %>%
  select("OTU", "Taxonomy") %>%
  rename_all(tolower) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="\\(\\d*\\)", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern=";$", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern='["]', replacement="")) %>%
  separate(taxonomy, into=c("kingdom", "phylum", "class", "order", "family", "genus"), 
           sep=";")

otu_rel_abund <- inner_join(metadata, otu_counts, by="sample") %>%
  inner_join(., taxonomy, by="otu") %>%
  group_by(sample) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  select(-count) %>%
  pivot_longer(c("kingdom", "phylum", "class", "order", "family", "genus", "otu"),
               names_to="level",
               values_to="taxon") %>%
  mutate(location = factor(season, 
                           levels=c("in",
                                    "off")))

order_rel_abund_mean <- otu_rel_abund %>%
  filter(level=="order") %>%
  group_by(combo, sample, taxon) %>%
  summarize(rel_abund = sum(rel_abund), 
            .groups="drop") %>%
  group_by(combo, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop") %>%
  mutate(taxon = str_replace(taxon,
                             "(.*)_unclassified", "Unclassified \\1"),
         taxon = str_replace_all(taxon, 
                                 "_", " "))

top_order_mean <- order_rel_abund_mean %>%
  group_by(taxon) %>%
  summarize(mean = mean(mean_rel_abund)) %>%
  arrange(mean) %>%
  top_n(11, mean) %>%
  pull(taxon)
view(top_order_mean)

order_rel_abund_mean %>%
  filter(taxon %in% top_order_mean) %>%
  mutate(taxon=factor(taxon, levels=top_order_mean)) %>%
  group_by(combo, taxon) %>%
  #print(n=100)
  ggplot(aes(x=combo, y=taxon, fill=mean_rel_abund)) +
  geom_tile(color="black") +
  scale_fill_gradient(low="white", high="red",
                      name=("Average Relative<br>Abundance (%)"),
                      limits=c(2,14), 
                      breaks = c(2, 6, 10, 14)) +
  #labels = c("0", "10", "20", "30", "40", ">50")) +
  coord_fixed() +
  scale_x_discrete(limits=c("in live", "in dead", "in none", "off live", "off dead", "off none"),
                   labels=c("In-season<br>Live Sclerotia<br>(n=12)", "In-season<br>Heat-killed<br>Sclerotia<br>(n=6)", "In-season<br>Bulk Soil<br>(n=4)", 
                            "Off-season<br>Live Sclerotia<br>(n=12)", "Off-season<br>Heat-killed<br>Sclerotia<br>(n=6)", "Off-season<br>Bulk Soil<br>(n=4)")) +
  scale_y_discrete(breaks=c("Unclassified Bacteria", "Gaiellales", "Caryophanales", "Rubrobacterales", "Hyphomicrobiales", "Solirubrobacterales", "Acidobacteria Gp6 incertae sedis",
                            "Rhodospirillales", "Unclassified Actinomycetota", "Gemmatimonadales", "Micromonosporales"),
                   labels=c("**Unclassified Bacteria (In)***", "**Gaiellales (Off)***", "**Caryophanales (In)***", "Rubrobacterales", "**Hyphomicrobiales (Off)***", 
                            "Solirubrobacterales", "**Acidobacteria Gp6 incertae sedis (Off)***", "Rhodospirillales", "**Unclassified Actinomycetota (In)***", 
                            "**Gemmatimonadales (In)***", "Micromonosporales")) +
  labs(title="Relative Abundance of the Top 11 Bacterial Orders<br>**Stiles Farm (Conducive)**",
       x=NULL,
       y=NULL) +
  theme_classic() +
  theme(axis.line = element_blank(), #gets rid of the axis lines
        axis.ticks = element_blank(), #gets rid of the axis tick marks
        axis.text.y = element_markdown(size=16),
        axis.text.x = element_markdown(size=14),
        legend.position = ("bottom"),
        legend.title = element_markdown(size=12),
        legend.text = element_text(size=9),
        legend.key.height = unit(11, "pt"),
        plot.title = element_markdown(hjust=0.5, size=18))

ggsave("heatmap.16S.order.stiles.ave.png", width=17, height=15)
