library(tidyverse)
library(readxl)
library(ggtext)
library(RColorBrewer)
library(MoMAColors)
library(glue)
library(dplyr)

display.all.moma(24, colorblind_only=F)

metadata <- read_excel("metadata.ITS.in.off.BF.xlsx") %>%
  mutate(combo = glue("{season} {treatment}"))

otu_counts <- read_tsv("final.opti_mcc.0.03.subsample.fungal.BF.shared") %>%
  select(Group, starts_with("Otu")) %>%
  rename(sample = Group) %>%
  pivot_longer(-sample, names_to="otu", values_to = "count")

taxonomy <- read_tsv(file="final.opti_mcc.0.03.cons.fungal.BF.taxonomy") %>%
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
  summarize(pool = max(rel_abund) < 0.019, #less than 2%
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
                    breaks=c("Eurotiales", "Pleosporales", "Sordariales", "Ascomycota order Incertae sedis", "Unclassified Fungi", "Unclassified Ascomycota", "Botryosphaeriales", 
                             "Sordariomycetes order Incertae sedis","Unclassified Sordariomycetes", "Capnodiales", "Onygenales", "Coniochaetales", "Xylariales", "Unclassified Dothideomycetes", 
                             "Hypocreales", "Other"),
                    labels=c("Eurotiales", "Pleosporales", "Sordariales", "Ascomycota order Incertae sedis", "Unclassified Fungi", "Unclassified Ascomycota", "Botryosphaeriales", 
                             "Sordariomycetes order Incertae sedis","Unclassified Sordariomycetes", "Capnodiales", "Onygenales", "Coniochaetales", "Xylariales", "Unclassified Dothideomycetes", 
                             "Hypocreales", "Other"),
                    #labels=c("**Eurotiales (In)***", "**Pleosporales (Off)***", "**Sordariales (In)***", "**Ascomycota order Incertae sedis (In)***", "Unclassified Fungi", 
                             #"Unclassified Ascomycota", "Botryosphaeriales", "Sordariomycetes order Incertae sedis", "Capnodiales", "Onygenales", "**Coniochaetales (Off)***", 
                             #"**Xylariales (Off)***", "Unclassified Dothideomycetes", "**Hypocreales (Off)***", "Other"),
                    values = c(moma.colors("Klein", 15, direction=1, type="continuous"), "dimgrey")) + #type could be continuous
  scale_x_discrete(limits=c("BI1A1F", "BI1A2F", "BI1B1F", "BI1B2F", "BI2A1F", "BI2A2F", "BI2B1F", "BI2B2F", "BI3A1F", "BI3A2F", "BI3B1F", "BI3B2F",
                            "BI1C1F", "BI1C2F", "BI2C1F", "BI2C2F", "BI3C1F", "BI3C2F",
                            "BIBA1F", "BIBA2F", "BIBB1F", "BIBB2F",
                            "BFO1A1F", "BFO1A2F", "BFO1B1F", "BFO1B2F", "BFO2A1F", "BFO2A2F", "BFO2B1F", "BFO2B2F", "BFO3A1F", "BFO3A2F", "BFO3B1F", "BFO3B2F",
                            "BFO1C1F", "BFO1C2F", "BFO2C1F", "BFO2C2F", "BFO3C1F", "BFO3C2F",
                            "BFOBA1F", "BFOBA2F", "BFOBB1F", "BFOBB2F"),
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
       title="Order-level Fungal Alpha Diversity: Off-season vs. 2024 Cotton-growing Season<br>**Bottom Farm**") +
  #labs(caption = "Other = relative abundance < 2%") +
  theme_classic() +
  theme(axis.text.x = element_markdown(),
        legend.text = element_markdown(),
        legend.title = element_text(size=10),
        legend.key.size = unit(10, "pt"),
        legend.position = "bottom",
        #plot.caption = element_text(hjust=0.8, size=8),
        plot.title = element_markdown(hjust=0.5))

ggsave("order.stackedbar.ITS.in.off.BF.png", width=22, height=9)

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
                    breaks=c("Eurotiales", "Pleosporales", "Sordariales", "Ascomycota order Incertae sedis", "Unclassified Fungi", "Unclassified Ascomycota", "Botryosphaeriales", 
                             "Unclassified Sordariomycetes", "Capnodiales", "Sordariomycetes order Incertae sedis", "Unclassified Dothideomycetes", "Hypocreales", "Other"),
                    labels=c("**Eurotiales (In)***", "**Pleosporales (Off)***", "**Sordariales (In)***", "**Ascomycota order Incertae sedis (In)***", "Unclassified Fungi", 
                             "Unclassified Ascomycota", "Botryosphaeriales", "Unclassified Sordariomycetes", "Capnodiales", "Sordariomycetes order Incertae sedis", 
                             "Unclassified Dothideomycetes", "**Hypocreales (Off)***", "Other"),
                    #values = c(brewer.pal(12, "Paired"), "navyblue", "dimgrey")) +
                    values = c(moma.colors("Warhol", 12, direction=1, type="continuous"), "dimgrey")) +  #type could be continuous
  scale_x_discrete(limits=c("in live", "in dead", "in none", "off live", "off dead", "off none"),
                   labels=c("In-season<br>Live (n=12)", "In-season<br>Heat-killed<br>(n=6)", "In-season<br>Bulk (n=4)",
                            "Off-season<br>Live (n=12)", "Off-season<br>Heat-killed<br>(n=6)", "Off-season<br>Bulk (n=4)")) +
  scale_y_continuous(expand=c(0,0)) + #this makes the bars touch the y axis
  labs(x=NULL,
       y="Average Relative Abundance (%)",
       title="Order-level Fungal Alpha Diversity<br>In- and Off-season Bottom Farm") +
  labs(caption = "Other = relative abundance < 1%") +
  theme_classic() +
  theme(axis.text.x = element_markdown(),
        legend.text = element_markdown(),
        legend.title = element_text(size=10),
        legend.key.size = unit(10, "pt"),
        legend.position = "bottom",
        plot.caption = element_text(hjust=0.95, size=8),
        plot.title = element_markdown(hjust=0.5))

ggsave("order.stackedbar.ITS.in.off.BF.ave.png", width=13, height=8)
