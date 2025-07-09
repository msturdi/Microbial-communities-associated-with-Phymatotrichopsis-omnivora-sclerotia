library(tidyverse)
library(readxl)
library(ggtext)
library(MoMAColors)
library(glue)
library(dplyr)

display.all.moma(24, colorblind_only=F)

metadata <- read_excel("oct.24.ITS.metadata.xlsx") %>%
  mutate(combo = glue("{location} {treatment}"))

otu_counts <- read_tsv("final.opti_mcc.0.03.subsample.oct.24.fungal.shared") %>%
  select(Group, starts_with("Otu")) %>%
  rename(sample = Group) %>%
  pivot_longer(-sample, names_to="otu", values_to = "count")

taxonomy <- read_tsv(file="final.opti_mcc.0.03.cons.oct.24.fungal.taxonomy") %>%
  select("OTU", "Taxonomy") %>%
  rename_all(tolower) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="\\(\\d*\\)", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern=";$", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="k__", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="p__", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="c__", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="o__", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="f__", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="g__", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="s__", replacement="")) %>%
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
         taxon = str_replace(taxon,
                             "unclassified_(.*)", "Unclassified \\1"),
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
         taxon = str_replace(taxon,
                             "unclassified_(.*)", "Unclassified \\1"),
         taxon = str_replace_all(taxon, 
                                 "_", " "))

order_pool <- order_rel_abund %>%
  group_by(taxon) %>%
  summarize(pool = max(rel_abund) < 0.029, #less than 1.9%
            mean = mean(rel_abund),
            .groups="drop")

order_pool_mean <- order_rel_abund_mean %>%
  group_by(taxon) %>%
  summarize(pool = max(mean_rel_abund) < 1, #less than 1% - this is not in decimal percentage because with the means I multiplied by 100
            mean = mean(mean_rel_abund),
            .groups="drop")

#order-level
inner_join(order_rel_abund, order_pool, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(sample, taxon) %>%
  summarize(rel_abund = 100*sum(rel_abund), 
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, mean, .desc=TRUE),
         taxon = fct_shift(taxon, n=1)) %>% #all of this mutate step is creating an "anchor" at the top and the bottom of the figure using the 2 most abundant phyla
  ggplot(aes(x=sample, y=rel_abund, fill=taxon)) +
  geom_col() +
  scale_fill_manual(name="Orders",
                    breaks=c("Hypocreales", "Pleosporales", "Sordariales", "Unclassified Ascomycota", "Ascomycota order Incertae sedis", "Capnodiales", 
                             "Sordariomycetes order Incertae sedis", "Unclassified Fungi", "Botryosphaeriales", "Unclassified Sordariomycetes", "Eurotiales", "Other"),
                    labels=c("Hypocreales", "Pleosporales", "Sordariales", "Unclassified Ascomycota", "Ascomycota order Incertae sedis", "Capnodiales", 
                             "Sordariomycetes order Incertae sedis", "Unclassified Fungi", "Botryosphaeriales", "Unclassified Sordariomycetes", "Eurotiales", "Other"),
                    values = c(moma.colors("Klein", 11, direction=1, type="continuous"), "dimgrey")) +
  scale_x_discrete(limits=c("BI1A1F", "BI1A2F", "BI1B1F", "BI1B2F", "BI2A1F", "BI2A2F", "BI2B1F", "BI2B2F", "BI3A1F", "BI3A2F", "BI3B1F", "BI3B2F",
                            "BI1C1F", "BI1C2F", "BI2C1F", "BI2C2F", "BI3C1F", "BI3C2F",
                            "BIBA1F", "BIBA2F", "BIBB1F", "BIBB2F",
                            "SI1A1F", "SI1A2F", "SI1B1F", "SI1B2F", "SI2A1F", "SI2A2F", "SI2B1F", "SI2B2F", "SI3A1F", "SI3A2F", "SI3B1F", "SI3B2F",
                            "SI1C1F", "SI1C2F", "SI2C1F", "SI2C2F", "SI3C1F", "SI3C2F",
                            "SIBA1F", "SIBA2F", "SIBB1F", "SIBB2F"),
                   labels=c("BF<br>Live<br>1", "BF<br>Live<br>2", "BF<br>Live<br>3", "BF<br>Live<br>4", "BF<br>Live<br>5","BF<br>Live<br>6", "BF<br>Live<br>7", "BF<br>Live<br>8", 
                            "BF<br>Live<br>9", "BF<br>Live<br>10", "BF<br>Live<br>11", "BF<br>Live<br>12",
                            "BF<br>Heat<br>Killed<br>1", "BF<br>Heat<br>Killed<br>2", "BF<br>Heat<br>Killed<br>3", "BF<br>Heat<br>Killed<br>4", "BF<br>Heat<br>Killed<br>5", 
                            "BF<br>Heat<br>Killed<br>6",
                            "BF<br>Bulk<br>Soil<br>1", "BF<br>Bulk<br>Soil<br>2", "BF<br>Bulk<br>Soil<br>3", "BF<br>Bulk<br>Soil<br>4",
                            "SF<br>Live<br>1", "SF<br>Live<br>2", "SF<br>Live<br>3", "SF<br>Live<br>4", "SF<br>Live<br>5","SF<br>Live<br>6", "SF<br>Live<br>7", "SF<br>Live<br>8", 
                            "SF<br>Live<br>9", "SF<br>Live<br>10", "SF<br>Live<br>11", "SF<br>Live<br>12",
                            "SF<br>Heat<br>Killed<br>1", "SF<br>Heat<br>Killed<br>2", "SF<br>Heat<br>Killed<br>3", "SF<br>Heat<br>Killed<br>4", "SF<br>Heat<br>Killed<br>5", 
                            "SF<br>Heat<br>Killed<br>6",
                            "SF<br>Bulk<br>Soil<br>1", "SF<br>Bulk<br>Soil<br>2", "SF<br>Bulk<br>Soil<br>3", "SF<br>Bulk<br>Soil<br>4")) +
  scale_y_continuous(expand=c(0,0)) +
  labs(x=NULL,
       y="Relative Abundance (%)",
       title="Order-level Fungal Alpha Diversity: **2024 Cotton-growing season**") +
  theme_classic() +
  theme(axis.text.x = element_markdown(),
        legend.text = element_markdown(),
        legend.title = element_text(size=10),
        legend.key.size = unit(10, "pt"),
        legend.position = "bottom",
        plot.title = element_markdown(hjust=0.5))

ggsave("order.stackedbar.ITS.png", width=18, height=8)
