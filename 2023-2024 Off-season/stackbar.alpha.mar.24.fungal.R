library(tidyverse)
library(readxl)
library(ggtext)
library(RColorBrewer)
library(MoMAColors)
library(glue)
library(dplyr)

display.all.moma(12, colorblind_only=F)

metadata <- read_excel("mar.24.fungal.metadata.xlsx") %>%
  mutate(combo = glue("{location} {treatment}"))

otu_counts <- read_tsv("final.opti_mcc.0.03.subsample.mar.24.fungal.shared") %>%
  select(Group, starts_with("Otu")) %>%
  rename(sample = Group) %>%
  pivot_longer(-sample, names_to="otu", values_to = "count")

taxonomy <- read_tsv(file="final.opti_mcc.0.03.cons.mar.24.fungal.taxonomy") %>%
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

#family-level analysis
family_rel_abund <- otu_rel_abund %>%
  filter(level=="family") %>%
  group_by(location, sample, taxon) %>%
  summarize(rel_abund = sum(rel_abund), 
            .groups="drop") %>%
  group_by(location, taxon) %>%
  mutate(taxon = str_replace(taxon,
                             "(.*)_unclassified", "Unclassified \\1"),
         taxon = str_replace(taxon,
                             "unclassified_(.*)", "Unclassified \\1"),
         taxon = str_replace_all(taxon, 
                                 "_", " ")) #gets ride of _ in taxon names

family_rel_abund_mean <- otu_rel_abund %>%
  filter(level=="family") %>%
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

family_pool <- family_rel_abund %>%
  group_by(taxon) %>%
  summarize(pool = max(rel_abund) < 0.04,
            mean = mean(rel_abund),
            .groups="drop")

family_pool_mean <- family_rel_abund_mean %>%
  group_by(taxon) %>%
  summarize(pool = max(mean_rel_abund) < 1, #less than 1% - this is not in decimal percentage because with the means I multiplied by 100
            mean = mean(mean_rel_abund),
            .groups="drop")

inner_join(family_rel_abund, family_pool, by="taxon") %>%
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
  scale_fill_manual(name="Families",
                    breaks=c("Trichocomaceae", "Pleosporaceae", "Pleosporales family Incertae sedis", "Sordariaceae", "Unclassified Fungi", "Unclassified Ascomycota", "Botryosphaeriaceae",
                             "Hypocreaceae", "Chaetomiaceae", "Lasiosphaeriaceae", "Ascomycota family Incertae sedis", "Unclassified Sordariomycetes", "Sporormiaceae", 
                             "Unclassified Sordariales", "Unclassified Onygenales", "Unclassified Dothideomycetes", "Nectriaceae", "Other"),
                    labels=c("Trichocomaceae", "Pleosporaceae", "Pleosporales family Incertae sedis", "Sordariaceae", "Unclassified Fungi", "Unclassified Ascomycota", "Botryosphaeriaceae",
                             "Hypocreaceae", "Chaetomiaceae", "Lasiosphaeriaceae", "Ascomycota family Incertae sedis", "Unclassified Sordariomycetes", "Sporormiaceae", 
                             "Unclassified Sordariales", "Unclassified Onygenales", "Unclassified Dothideomycetes", "Nectriaceae", "Other"),
                    #values = c(brewer.pal(12, "Paired"), "navyblue", "dimgrey")) +
                    values = c(moma.colors("Warhol", 17, direction=1, type="continuous"), "dimgrey")) + #type could be continuous
  scale_x_discrete(limits=c("BFO1A1F", "BFO1A2F", "BFO1B1F", "BFO1B2F", "BFO2A1F", "BFO2A2F", "BFO2B1F", "BFO2B2F", "BFO3A1F", "BFO3A2F", "BFO3B1F", "BFO3B2F",
                            "BFO1C1F", "BFO1C2F", "BFO2C1F", "BFO2C2F", "BFO3C1F", "BFO3C2F",
                            "BFOBA1F", "BFOBA2F", "BFOBB1F", "BFOBB2F",
                            "SO1A1F", "SO1A2F", "SO1B1F", "SO1B2F", "SO2A1F", "SO2A2F", "SO2B1F", "SO2B2F", "SO3A1F", "SO3A2F", "SO3B1F", "SO3B2F",
                            "SO1C1F", "SO1C2F", "SO2C1F", "SO2C2F", "SO3C1F", "SO3C2F",
                            "SOBA1F", "SOBA2F", "SOBB1F", "SOBB2F"),
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
       title="Family-level Fungal Alpha Diversity") +
  labs(caption = "Other = relative abundance < 4%") +
  theme_classic() +
  theme(axis.text.x = element_markdown(),
        legend.text = element_markdown(),
        legend.title = element_text(size=10),
        legend.key.size = unit(10, "pt"),
        plot.caption = element_text(hjust=1.14, size=8),
        plot.title = element_text(hjust=0.5))


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
  summarize(pool = max(rel_abund) < 0.04, #less than 4%
            mean = mean(rel_abund),
            .groups="drop")

order_pool_mean <- order_rel_abund_mean %>%
  group_by(taxon) %>%
  summarize(pool = max(mean_rel_abund) < 1, #less than 1% - this is not in decimal percentage because with the means I multiplied by 100
            mean = mean(mean_rel_abund),
            .groups="drop")

#order-level all samples
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
                    breaks=c("Eurotiales", "Pleosporales", "Sordariales", "Unclassified Fungi", "Unclassified Ascomycota", "Botryosphaeriales",
                             "Ascomycota order Incertae sedis", "Unclassified Sordariomycetes", "Xylariales", "Onygenales", "Unclassified Dothideomycetes",
                             "Hypocreales", "Other"),
                    labels=c("Eurotiales", "Pleosporales", "Sordariales", "Unclassified Fungi", "Unclassified Ascomycota", "Botryosphaeriales",
                             "Ascomycota order Incertae sedis", "Unclassified Sordariomycetes", "Xylariales", "Onygenales", "Unclassified Dothideomycetes",
                             "Hypocreales", "Other"),
                    #values = c(brewer.pal(12, "Paired"), "navyblue", "dimgrey")) +
                    values = c(moma.colors("Klein", 12, direction=1, type="continuous"), "dimgrey")) + #type could be continuous
  scale_x_discrete(limits=c("BFO1A1F", "BFO1A2F", "BFO1B1F", "BFO1B2F", "BFO2A1F", "BFO2A2F", "BFO2B1F", "BFO2B2F", "BFO3A1F", "BFO3A2F", "BFO3B1F", "BFO3B2F",
                            "BFO1C1F", "BFO1C2F", "BFO2C1F", "BFO2C2F", "BFO3C1F", "BFO3C2F",
                            "BFOBA1F", "BFOBA2F", "BFOBB1F", "BFOBB2F",
                            "SO1A1F", "SO1A2F", "SO1B1F", "SO1B2F", "SO2A1F", "SO2A2F", "SO2B1F", "SO2B2F", "SO3A1F", "SO3A2F", "SO3B1F", "SO3B2F",
                            "SO1C1F", "SO1C2F", "SO2C1F", "SO2C2F", "SO3C1F", "SO3C2F",
                            "SOBA1F", "SOBA2F", "SOBB1F", "SOBB2F"),
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
       title="Order-level Fungal Alpha Diversity: **2023-2024 Off-season**") +
  #labs(caption = "Other = relative abundance < 4%") +
  theme_classic() +
  theme(axis.text.x = element_markdown(),
        legend.text = element_markdown(),
        legend.title = element_text(size=10),
        legend.key.size = unit(10, "pt"),
        legend.position = "bottom",
        #plot.caption = element_text(hjust=0.8, size=8),
        plot.title = element_markdown(hjust=0.5))

ggsave("order.stackedbar.ITS.png", width=18, height=8)

#order-level means by combo
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
                    breaks=c("Eurotiales", "Pleosporales", "Sordariales", "Unclassified Fungi", "Unclassified Ascomycota", "Botryosphaeriales",
                             "Ascomycota order Incertae sedis", "Unclassified Sordariomycetes", "Xylariales", "Onygenales", "Unclassified Dothideomycetes",
                             "Hypocreales", "Other"),
                    labels=c("**Eurotiales (S)***", "Pleosporales", "**Sordariales (S)***", "Unclassified Fungi", "**Unclassified Ascomycota (BF)***", "Botryosphaeriales",
                             "Ascomycota order Incertae sedis", "**Unclassified Sordariomycetes (S)***", "**Xylariales (S)***", "**Onygenales (BF)***", "Unclassified Dothideomycetes",
                             "**Hypocreales (BF)***", "Other"),
                    #values = c(brewer.pal(12, "Paired"), "navyblue", "dimgrey")) +
                    values = c(moma.colors("Warhol", 12, direction=1, type="continuous"), "dimgrey")) +
  scale_x_discrete(limits=c("BF live", "BF dead", "BF none", "Stiles live", "Stiles dead", "Stiles none"),
                   labels=c("Bottom Farm<br>Live (n=12)", "Bottom Farm<br>Heat-killed<br>(n=6)", "Bottom Farm<br>Bulk (n=4)",
                            "Stiles Farm<br>Live (n=12)", "Stiles Farm<br>Heat-killed<br>(n=6)", "Stiles Farm<br>Bulk (n=4)")) +
  scale_y_continuous(expand=c(0,0)) +
  labs(x=NULL,
       y="Average Relative Abundance (%)",
       title="Order-level Fungal Alpha Diversity") +
  labs(caption = "Other = relative abundance < 1%") +
  theme_classic() +
  theme(axis.text.x = element_markdown(),
        legend.text = element_markdown(),
        legend.title = element_text(size=10),
        legend.key.size = unit(10, "pt"),
        plot.caption = element_text(hjust=1.24, size=8),
        plot.title = element_text(hjust=0.5))

ggsave("order.stackedbar.ITS.trt.ave.png", width=12, height=8)

#class-level
class_rel_abund <- otu_rel_abund %>%
  filter(level=="class") %>%
  group_by(location, sample, taxon) %>%
  summarize(rel_abund = sum(rel_abund), 
            .groups="drop") %>%
  group_by(location, taxon) %>%
  mutate(taxon = str_replace(taxon,
                             "(.*)_unclassified", "Unclassified \\1"),
         taxon = str_replace_all(taxon, 
                                 "_", " ")) #gets ride of _ in taxon names

class_rel_abund_mean <- otu_rel_abund %>%
  filter(level=="class") %>%
  group_by(combo, sample, taxon) %>%
  summarize(rel_abund = sum(rel_abund), 
            .groups="drop") %>%
  group_by(combo, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop") %>%
  mutate(taxon = str_replace(taxon,
                             "(.*)_unclassified", "Unclassified \\1"),
         taxon = str_replace_all(taxon, 
                                 "_", " "))

class_pool <- class_rel_abund %>%
  group_by(taxon) %>%
  summarize(pool = max(rel_abund) < 0.02, 
            mean = mean(rel_abund),
            .groups="drop")

class_pool_mean <- class_rel_abund_mean %>%
  group_by(taxon) %>%
  summarize(pool = max(mean_rel_abund) < 1, #less than 1% - this is not in decimal percentage because with the means I multiplied by 100
            mean = mean(mean_rel_abund),
            .groups="drop")

#showing all samples
inner_join(class_rel_abund, class_pool, by="taxon") %>%
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
  scale_fill_manual(name="Classes",
                    breaks=c("Thermoleophilia", "Actinobacteria", "Unclassified Bacteria", "Bacilli", "Acidobacteria Gp6", "Rubrobacteria",
                             "Deltaproteobacteria", "Unclassified Actinomycetota", "Betaproteobacteria", "Acidimicrobiia", "Gammaproteobacteria",
                             "Gemmatimonadia", "Acidobacteria Gp16", "Chitinophagia", "Acidobacteria Gp4", "Alphaproteobacteria", "Other"),
                    labels=c("Thermoleophilia", "Actinobacteria", "Unclassified Bacteria", "Bacilli", "Acidobacteria Gp6", "Rubrobacteria",
                             "Deltaproteobacteria", "Unclassified Actinomycetota", "Betaproteobacteria", "Acidimicrobiia", "Gammaproteobacteria",
                             "Gemmatimonadia", "Acidobacteria Gp16", "Chitinophagia", "Acidobacteria Gp4", "Alphaproteobacteria", "Other"),
                    #values = c(brewer.pal(12, "Paired"), "navyblue", "dimgrey")) +
                    values = c(moma.colors("Warhol", 16, direction=1, type="continuous"), "dimgrey")) + #type could be continuous
  scale_x_discrete(limits=c("BFO1A1B", "BFO1A2B", "BFO1B1B", "BFO1B2B", "BFO2A1B", "BFO2A2B", "BFO2B1B", "BFO2B2B", "BFO3A1B", "BFO3A2B", "BFO3B1B", "BFO3B2B",
                            "BFO1C1B", "BFO1C2B", "BFO2C1B", "BFO2C2B", "BFO3C1B", "BFO3C2B",
                            "BFOBA1B", "BFOBA2B", "BFOBB1B", "BFOBB2B",
                            "SO1A1B", "SO1A2B", "SO1B1B", "SO1B2B", "SO2A1B", "SO2A2B", "SO2B1B", "SO2B2B", "SO3A1B", "SO3A2B", "SO3B1B", "SO3B2B",
                            "SO1C1B", "SO1C2B", "SO2C1B", "SO2C2B", "SO3C1B", "SO3C2B",
                            "SOBA1B", "SOBA2B", "SOBB1B", "SOBB2B"),
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
       title="Class-level 16S rRNA Gene Alpha Diversity") +
  labs(caption = "Other = relative abundance < 2%") +
  theme_classic() +
  theme(axis.text.x = element_markdown(),
        legend.text = element_markdown(),
        legend.title = element_text(size=10),
        legend.key.size = unit(10, "pt"),
        plot.caption = element_text(hjust=1.14, size=8),
        plot.title = element_text(hjust=0.5))

ggsave("class.stackedbar.16S.png", width=18, height=8)

#showing means by combo
inner_join(class_rel_abund_mean, class_pool_mean, by="taxon") %>%
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
  scale_fill_manual(name="Classes",
                    breaks=c("Alphaproteobacteria", "Actinobacteria", "Unclassified Bacteria", "Bacilli", "Acidobacteria Gp6", "Rubrobacteria",
                             "Unclassified Actinomycetota", "Deltaproteobacteria", "Betaproteobacteria", "Acidimicrobiia", "Gemmatimonadia",
                             "Gammaproteobacteria", "Acidobacteria Gp16", "Chitinophagia", "Acidobacteria Gp4", "Thermoleophilia", "Other"),
                    labels=c("Alphaproteobacteria", "Actinobacteria", "Unclassified Bacteria", "Bacilli", "Acidobacteria Gp6", "Rubrobacteria",
                             "Unclassified Actinomycetota", "Deltaproteobacteria", "Betaproteobacteria", "Acidimicrobiia", "Gemmatimonadia",
                             "Gammaproteobacteria", "Acidobacteria Gp16", "Chitinophagia", "Acidobacteria Gp4", "Thermoleophilia", "Other"),
                    #values = c(brewer.pal(12, "Paired"), "navyblue", "dimgrey")) +
                    values = c(moma.colors("Warhol", 16, direction=1, type="continuous"), "dimgrey")) + #type could be continuous
  scale_x_discrete(limits=c("BF live", "BF dead", "BF none", "Stiles live", "Stiles dead", "Stiles none"),
                   labels=c("Bottom Farm<br>Live (n=12)", "Bottom Farm<br>Heat-killed<br>(n=6)", "Bottom Farm<br>Bulk (n=4)",
                            "Stiles Farm<br>Live (n=12)", "Stiles Farm<br>Heat-killed<br>(n=6)", "Stiles Farm<br>Bulk (n=4)")) +
  scale_y_continuous(expand=c(0,0)) +
  labs(x=NULL,
       y="Average Relative Abundance (%)",
       title="Class-level 16S rRNA Gene Alpha Diversity") +
  labs(caption = "Other = relative abundance < 1%") +
  theme_classic() +
  theme(axis.text.x = element_markdown(),
        legend.text = element_markdown(),
        legend.title = element_text(size=10),
        legend.key.size = unit(10, "pt"),
        plot.caption = element_text(hjust=1.24, size=8),
        plot.title = element_text(hjust=0.5))

ggsave("class.stackedbar.16S.trt.ave.png", width=12, height=8)
