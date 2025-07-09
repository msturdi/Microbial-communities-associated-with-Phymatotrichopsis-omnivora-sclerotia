library(tidyverse)
library(readxl)
library(ggtext)
library(RColorBrewer)
library(MoMAColors)
library(glue)
library(dplyr)

display.all.moma(24, colorblind_only=F)

metadata <- read_excel("metadata.ITS.in.off.stiles.xlsx") %>%
  mutate(combo = glue("{season} {treatment}"))

otu_counts <- read_tsv("final.opti_mcc.0.03.subsample.fungal.stiles.shared") %>%
  select(Group, starts_with("Otu")) %>%
  rename(sample = Group) %>%
  pivot_longer(-sample, names_to="otu", values_to = "count")

taxonomy <- read_tsv(file="final.opti_mcc.0.03.cons.fungal.stiles.taxonomy") %>%
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
  mutate(season = factor(season, 
                         levels=c("in",
                                  "off")))

order_rel_abund <- otu_rel_abund %>%
  filter(level=="order") %>%
  group_by(season, sample, taxon) %>%
  summarize(rel_abund = sum(rel_abund), 
            .groups="drop") %>%
  group_by(season, taxon) %>%
  mutate(taxon = str_replace(taxon,
                             "(.*)_unclassified", "Unclassified \\1"),
         taxon = str_replace(taxon, 
                             "unclassified_(.*)", "Unclassified \\1"),
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
         taxon = str_replace(taxon, 
                             "unclassified_(.*)", "Unclassified \\1"),
         taxon = str_replace_all(taxon, 
                                 "_", " "))

order_pool <- order_rel_abund %>%
  group_by(taxon) %>%
  summarize(pool = max(rel_abund) < 0.025, #less than 2.5%
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
                    breaks=c("Hypocreales", "Pleosporales", "Sordariales", "Capnodiales", "Unclassified Fungi", "Botryosphaeriales", "Xylariales", "Unclassified Ascomycota", 
                             "Unclassified Sordariomycetes", "Sordariomycetes order Incertae sedis", "Eurotiales", "Other"),
                    labels=c("Hypocreales", "Pleosporales", "Sordariales", "Capnodiales", "Unclassified Fungi", "Botryosphaeriales", "Xylariales", "Unclassified Ascomycota", 
                             "Unclassified Sordariomycetes", "Sordariomycetes order Incertae sedis", "Eurotiales", "Other"),
                    #labels=c("**Hypocreales (Off)***", "**Pleosporales (Off)***", "Sordariales", "**Capnodiales (Off)***", "**Unclassified Fungi (Off)***", "Botryosphaeriales", 
                             #"**Xylariales (Off)***", "Unclassified Ascomycota", "**Unclassified Sordariomycetes (In)***", "**Sordariomycetes order Incertae sedis (Off)***", 
                             #"**Eurotiales (In)***", "Other"),
                    values = c(moma.colors("Klein", 11, direction=1, type="continuous"), "dimgrey")) + #type could be continuous
  scale_x_discrete(limits=c("SI1A1F", "SI1A2F", "SI1B1F", "SI1B2F", "SI2A1F", "SI2A2F", "SI2B1F", "SI2B2F", "SI3A1F", "SI3A2F", "SI3B1F", "SI3B2F",
                            "SI1C1F", "SI1C2F", "SI2C1F", "SI2C2F", "SI3C1F", "SI3C2F",
                            "SIBA1F", "SIBA2F", "SIBB1F", "SIBB2F",
                            "SO1A1F", "SO1A2F", "SO1B1F", "SO1B2F", "SO2A1F", "SO2A2F", "SO2B1F", "SO2B2F", "SO3A1F", "SO3A2F", "SO3B1F", "SO3B2F",
                            "SO1C1F", "SO1C2F", "SO2C1F", "SO2C2F", "SO3C1F", "SO3C2F",
                            "SOBA1F", "SOBA2F", "SOBB1F", "SOBB2F"),
                   labels=c("In<br>Live<br>1", "In<br>Live<br>2", "In<br>Live<br>3", "In<br>Live<br>4", "In<br>Live<br>5", "In<br>Live<br>6", "In<br>Live<br>7", "In<br>Live<br>8", 
                            "In<br>Live<br>9", "In<br>Live<br>10", "In<br>Live<br>11", "In<br>Live<br>12",
                            "In<br>Heat<br>killed<br>1", "In<br>Heat<br>killed<br>2", "In<br>Heat<br>killed<br>3", "In<br>Heat<br>killed<br>4", "In<br>Heat<br>killed<br>5", 
                            "In<br>Heat<br>killed<br>6",
                            "In<br>Bulk<br>1", "In<br>Bulk<br>2", "In<br>Bulk<br>3", "In<br>Bulk<br>4",
                            "Off<br>Live<br>1", "Off<br>Live<br>2", "Off<br>Live<br>3", "Off<br>Live<br>4", "Off<br>Live<br>5", "Off<br>Live<br>6", "Off<br>Live<br>7", "Off<br>Live<br>8", 
                            "Off<br>Live<br>9", "Off<br>Live<br>10", "Off<br>Live<br>11", "Off<br>Live<br>12",
                            "Off<br>Heat<br>killed<br>1", "Off<br>Heat<br>killed<br>2", "Off<br>Heat<br>killed<br>3", "Off<br>Heat<br>killed<br>4", "Off<br>Heat<br>killed<br>5", 
                            "Off<br>Heat<br>killed<br>6",
                            "Off<br>Bulk<br>1", "Off<br>Bulk<br>2", "Off<br>Bulk<br>3", "Off<br>Bulk<br>4")) +
  scale_y_continuous(expand=c(0,0)) +
  labs(x=NULL,
       y="Relative Abundance (%)",
       title="Order-level Fungal Alpha Diversity: Off-season vs. 2024 Cotton-growing Season<br>**Stiles Farm**") +
  #labs(caption = "Other = relative abundance < 2.5%") +
  theme_classic() +
  theme(axis.text.x = element_markdown(),
        legend.text = element_markdown(),
        legend.title = element_text(size=10),
        legend.key.size = unit(10, "pt"),
        legend.position = "bottom",
        #plot.caption = element_text(hjust=0.8, size=8),
        plot.title = element_markdown(hjust=0.5))

ggsave("order.stackedbar.ITS.in.off.stiles.png", width=22, height=9)

#order-level averages
inner_join(order_rel_abund_mean, order_pool_mean, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(combo, taxon) %>%
  summarize(mean_rel_abund = sum(mean_rel_abund), 
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, mean, .desc=TRUE),
         taxon = fct_shift(taxon, n=1)) %>% #all of this mutate step is creating an "anchor" at the top and the bottom of the figure using the 2 most abundant phyla
  ggplot(aes(x=combo, y=mean_rel_abund, fill=taxon)) +
  geom_col() +
  scale_fill_manual(name="Orders",
                    breaks=c("Hypocreales", "Pleosporales", "Sordariales", "Capnodiales", "Unclassified Fungi", "Botryosphaeriales", "Xylariales", "Unclassified Ascomycota", 
                             "Unclassified Sordariomycetes", "Sordariomycetes order Incertae sedis", "Eurotiales", "Other"),
                    labels=c("**Hypocreales (Off)***", "**Pleosporales (Off)***", "Sordariales", "**Capnodiales (Off)***", "**Unclassified Fungi (Off)***", "Botryosphaeriales", 
                             "**Xylariales (Off)***", "Unclassified Ascomycota", "**Unclassified Sordariomycetes (In)***", "**Sordariomycetes order Incertae sedis (Off)***", 
                             "**Eurotiales (In)***", "Other"),
                    values = c(moma.colors("Warhol", 11, direction=1, type="continuous"), "dimgrey")) +  #type could be continuous
  scale_x_discrete(limits=c("in live", "in dead", "in none", "off live", "off dead", "off none"),
                   labels=c("In-season<br>Live (n=12)", "In-season<br>Heat-killed<br>(n=6)", "In-season<br>Bulk (n=4)",
                            "Off-season<br>Live (n=12)", "Off-season<br>Heat-killed<br>(n=6)", "Off-season<br>Bulk (n=4)")) +
  scale_y_continuous(expand=c(0,0)) + #this makes the bars touch the y axis
  labs(x=NULL,
       y="Average Relative Abundance (%)",
       title="Order-level Fungal Alpha Diversity<br>In- and Off-season Stiles Farm") +
  labs(caption = "Other = relative abundance < 1%") +
  theme_classic() +
  theme(axis.text.x = element_markdown(),
        legend.text = element_markdown(),
        legend.title = element_text(size=10),
        legend.key.size = unit(10, "pt"),
        legend.position = "bottom",
        plot.caption = element_text(hjust=0.95, size=8),
        plot.title = element_markdown(hjust=0.5))

ggsave("order.stackedbar.ITS.in.off.stiles.ave.png", width=13, height=8)
