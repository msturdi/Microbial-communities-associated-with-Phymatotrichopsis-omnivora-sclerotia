library(tidyverse)
library(readxl)
library(ggtext)
library(glue)
library(dplyr)
library(scales)

metadata <- read_excel("mar.24.bact.metadata.xlsx") %>%
  mutate(combo = glue("{location} {treatment}"))

otu_counts <- read_tsv("final.opti_mcc.0.03.subsample.mar.24.bact.shared") %>%
  select(Group, starts_with("Otu")) %>%
  rename(sample = Group) %>%
  pivot_longer(-sample, names_to="otu", values_to = "count")

taxonomy <- read_tsv(file="final.opti_mcc.0.03.cons.mar.24.bact.taxonomy") %>%
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
  mutate(location = factor(location, 
                           levels=c("BF",
                                    "Stiles")))

order_rel_abund <- otu_rel_abund %>%
  filter(level=="order") %>%
  group_by(location, sample, taxon) %>%
  summarize(rel_abund = sum(rel_abund), 
            .groups="drop") %>%
  group_by(location, taxon) %>%
  mutate(taxon = str_replace(taxon,
                             "(.*)_unclassified", "Unclassified \\1"),
         taxon = str_replace_all(taxon, 
                                 "_", " "))

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

order_rel_abund %>%
  group_by(taxon) %>%
  summarize(rel_abund=sum(rel_abund)) %>%
  arrange(rel_abund) %>%
  view()

top_order <- order_rel_abund %>%
  group_by(taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  arrange(rel_abund) %>%
  top_n(10, rel_abund) %>%
  pull(taxon)
view(top_order)

top_order_mean <- order_rel_abund_mean %>%
  group_by(taxon) %>%
  summarize(mean = mean(mean_rel_abund)) %>%
  arrange(mean) %>%
  top_n(10, mean) %>%
  pull(taxon)
view(top_order_mean)

order_rel_abund %>%
  filter(taxon %in% top_order) %>%
  mutate(taxon=factor(taxon, levels=top_order)) %>%
  mutate(rel_abund=100 * rel_abund)

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
  scale_x_discrete(limits=c("BF live", "BF dead", "BF none", "Stiles live", "Stiles dead", "Stiles none"),
                   labels=c("BF<br>Live Sclerotia<br>(n=12)", "BF<br>Heat-killed<br>Sclerotia<br>(n=6)", "BF<br>Bulk Soil<br>(n=4)", "SF<br>Live Sclerotia<br>(n=12)", 
                            "SF<br>Heat-killed<br>Sclerotia<br>(n=6)", "SF<br>Bulk Soil<br>(n=4)")) +
  scale_y_discrete(breaks=c("Gaiellales", "Hyphomicrobiales", "Unclassified Bacteria", "Caryophanales", "Acidobacteria Gp6 incertae sedis", "Solirubrobacterales", "Rhodospirillales",
                            "Rubrobacterales", "Micromonosporales", "Propionibacteriales"),
                   labels=c("**Gaiellales (SF)***",  "**Hyphomicrobiales (BF)***", "**Unclassified Bacteria (SF)***", "Caryophanales", "Acidobacteria Gp6 incertae sedis", 
                            "**Solirubrobacterales (SF)***", "**Rhodospirillales (BF)***", "**Rubrobacterales (SF)***", "Micromonosporales", "**Propionibacteriales (BF)***")) +
  labs(title="Relative Abundance of the Top 10 Bacterial Orders<br>**2023-2024 Off-season**",
       x=NULL,
       y=NULL) +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_markdown(size=12),
        axis.text.x = element_markdown(size=12),
        legend.position = "bottom",
        legend.title = element_markdown(size=11),
        legend.text = element_text(size=9),
        legend.key.height = unit(11, "pt"),
        plot.title = element_markdown(hjust=0.5, size=14))

ggsave("heatmap.16S.order.ave.png", width=12, height=12)
