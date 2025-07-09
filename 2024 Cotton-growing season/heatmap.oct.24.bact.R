library(tidyverse)
library(readxl)
library(ggtext)
library(RColorBrewer)
library(MoMAColors)
library(glue)
library(dplyr)
library(scales)

display.all.moma(24, colorblind_only=F)

metadata <- read_excel("oct.24.16S.metadata.xlsx") %>%
  mutate(combo = glue("{location} {treatment}"))

shared <- read_tsv("final.opti_mcc.0.03.subsample.oct.24.bacterial.shared",
                   col_types = cols(Group = col_character(),
                                    .default = col_double()))

meta_cols <- c("label", "Group", "numOtus") #define columns to keep from shared file
otu_sums <- colSums(shared[, !(names(shared) %in% meta_cols)]) #sum OTU abundances across all samples
otus_to_keep <- names(otu_sums[otu_sums > 1]) #keep OTUs with total counts above 1

shared_filtered <- shared %>%
  select(all_of(meta_cols), all_of(otus_to_keep))

otu_counts <- shared_filtered %>%
  select(Group, starts_with("Otu")) %>%
  rename(sample = Group) %>%
  pivot_longer(-sample, names_to="otu", values_to = "count")

##
#no filtering of OTUs
##
otu_counts <- read_tsv("final.opti_mcc.0.03.subsample.oct.24.bacterial.shared") %>%
  select(Group, starts_with("Otu")) %>%
  rename(sample = Group) %>%
  pivot_longer(-sample, names_to="otu", values_to = "count")

taxonomy <- read_tsv(file="final.opti_mcc.0.03.cons.oct.24.bacterial.taxonomy") %>%
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
                                 "_", " ")) #gets ride of _ in taxon names

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
                   labels=c("BF<br>Live Sclerotia<br>(n=12)", "BF<br>Heat-killed<br>Sclerotia<br>(n=6)", "BF<br>Bulk Soil<br>(n=4)", 
                            "SF<br>Live Sclerotia<br>(n=12)", "SF<br>Heat-killed<br>Sclerotia<br>(n=6)", "SF<br>Bulk Soil<br>(n=4)")) +
  scale_y_discrete(breaks=c("Unclassified Bacteria", "Caryophanales", "Rubrobacterales", "Hyphomicrobiales", "Solirubrobacterales", "Gaiellales", "Rhodospirillales", 
                            "Acidobacteria Gp6 incertae sedis","Unclassified Actinomycetota", "Micromonosporales"),
                   labels=c("**Unclassified Bacteria (SF)***", "**Caryophanales (SF)***", "**Rubrobacterales (SF)***", "Hyphomicrobiales", "Solirubrobacterales", "**Gaiellales (SF)***", 
                            "**Rhodospirillales (BF)***", "**Acidobacteria Gp6 incertae sedis (BF)***", "**Unclassified Actinomycetota (SF)***", "Micromonosporales")) +
  labs(title="Relative Abundance of the Top 10 Bacterial Orders<br>**2024 Cotton-growing Season**",
       x=NULL,
       y=NULL) +
  theme_classic() +
  theme(axis.line = element_blank(), #gets rid of the axis lines
        axis.ticks = element_blank(), #gets rid of the axis tick marks
        axis.text.y = element_markdown(size=12),
        axis.text.x = element_markdown(size=12),
        legend.position = "bottom",
        legend.title = element_markdown(size=11),
        legend.text = element_text(size=9),
        legend.key.height = unit(11, "pt"),
        plot.title = element_markdown(hjust=0.5, size=14))

ggsave("heatmap.16S.order.ave.png", width=12, height=12)
