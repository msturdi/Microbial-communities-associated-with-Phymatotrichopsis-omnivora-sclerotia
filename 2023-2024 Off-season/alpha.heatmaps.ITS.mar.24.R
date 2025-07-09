library(tidyverse)
library(readxl)
library(ggtext)
library(RColorBrewer)
library(MoMAColors)
library(glue)
library(scales)

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
  mutate(taxon=str_replace(taxon, "unclassified_(.*)", "Unclassified \\1"),
         taxon=str_replace(taxon, "(.*)_unclassified", "Unclassified \\1"),
         taxon=str_replace_all(taxon, 
                                 "_", " ")) %>%
  mutate(location = factor(location, 
                           levels=c("BF",
                                    "Stiles")))

#order-level relative abundance
order_rel_abund <- otu_rel_abund %>%
  filter(level=="order") %>%
  group_by(location, sample, taxon) %>%
  summarize(rel_abund = sum(rel_abund), 
            .groups="drop") %>%
  group_by(location, taxon)

order_rel_abund_mean <- otu_rel_abund %>%
  filter(level=="order") %>%
  group_by(combo, sample, taxon) %>%
  summarize(rel_abund = sum(rel_abund), 
            .groups="drop") %>%
  group_by(combo, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop")

order_rel_abund %>%
  group_by(taxon) %>%
  summarize(rel_abund=sum(rel_abund)) %>%
  arrange(rel_abund) %>%
  view()

top_order <- order_rel_abund %>%
  group_by(taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  arrange(rel_abund) %>%
  top_n(7, rel_abund) %>%
  pull(taxon)
view(top_order)

top_order_mean <- order_rel_abund_mean %>%
  group_by(taxon) %>%
  summarize(mean = mean(mean_rel_abund)) %>%
  arrange(mean) %>%
  top_n(7, mean) %>%
  pull(taxon)
view(top_order_mean)

order_rel_abund %>%
  filter(taxon %in% top_order) %>%
  mutate(taxon=factor(taxon, levels=top_order)) %>%
  mutate(rel_abund=100 * rel_abund)

order_rel_abund %>%
  filter(taxon %in% top_order) %>%
  mutate(taxon=factor(taxon, levels=top_order)) %>%
  mutate(rel_abund= 100 * (rel_abund + 1/20000)) %>%
  ggplot(aes(x=sample, y=taxon, fill=rel_abund)) +
  geom_tile(color="black") +
  scale_fill_gradient(low="white", high="red",
                      name=("Relative<br>Abundance (%)"),
                      limits=c(0,50), 
                      oob=squish, #part of "scales" library
                      breaks = c(0, 10, 20, 30, 40, 50), labels = c("0", "10", "20", "30", "40", ">50")) + 
  coord_fixed() +
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
  scale_y_discrete(breaks=c("Hypocreales", "Eurotiales", "Pleosporales", "Sordariales", "Unclassified Fungi", "Unclassified Ascomycota", "Botryosphaeriales"),
                   labels=c("**Hypocreales (BF)***", "**Eurotiales (SF)***", "Pleosporales", "**Sordariales (SF)***", "Unclassified Fungi", "**Unclassified Ascomycota (BF)***", 
                            "Botryosphaeriales")) +
  labs(title="Relative Abundance of the Top 7 Orders",
       x=NULL,
       y=NULL) +
  theme_classic() +
  theme(axis.line = element_blank(), #gets rid of the axis lines
        axis.ticks = element_blank(), #gets rid of the axis tick marks
        axis.text.y = element_markdown(size=8),
        axis.text.x = element_markdown(size=8),
        #legend.title.align = 0.5,
        legend.title = element_markdown(size=8),
        legend.text = element_text(size=7.5),
        legend.key.height = unit(11, "pt"),
        plot.title = element_text(hjust=0.5))

ggsave("heatmap.ITS.order.png", width=16, height=6)

order_rel_abund_mean %>%
  filter(taxon %in% top_order_mean) %>%
  mutate(taxon=factor(taxon, levels=top_order_mean)) %>%
  group_by(combo, taxon) %>%
  #print(n=100)
  ggplot(aes(x=combo, y=taxon, fill=mean_rel_abund)) +
  geom_tile(color="black") +
  scale_fill_gradient(low="white", high="red",
                      name=("Average Relative<br>Abundance (%)"),
                      limits=c(0,50), 
                      oob=squish, #part of "scales" library
                      breaks = c(0, 10, 20, 30, 40, 50), labels = c("0", "10", "20", "30", "40", ">50")) + 
  coord_fixed() +
  scale_x_discrete(limits=c("BF live", "BF dead", "BF none", "Stiles live", "Stiles dead", "Stiles none"),
                   labels=c("BF<br>Live Sclerotia<br>(n=12)", "BF<br>Heat-killed<br>Sclerotia<br>(n=6)", "BF<br>Bulk Soil<br>(n=4)", "Stiles<br>Live Sclerotia<br>(n=12)", 
                            "Stiles<br>Heat-killed<br>Sclerotia<br>(n=6)", "Stiles<br>Bulk Soil<br>(n=4)")) +
  scale_y_discrete(breaks=c("Hypocreales", "Eurotiales", "Pleosporales", "Sordariales", "Unclassified Fungi", "Unclassified Ascomycota", "Botryosphaeriales"),
                   labels=c("**Hypocreales (BF)***", "**Eurotiales (SF)***", "Pleosporales", "**Sordariales (SF)***", "Unclassified Fungi", "**Unclassified Ascomycota (BF)***", 
                            "Botryosphaeriales")) +
  labs(title="Relative Abundance of the Top 7 Fungal Orders<br>**2023-2024 Off-season**",
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

ggsave("heatmap.ITS.order.ave.png", width=12, height=10)

#genus-level relative abundance
genus_rel_abund <- otu_rel_abund %>%
  filter(level=="genus") %>%
  group_by(location, sample, taxon) %>%
  summarize(rel_abund = sum(rel_abund), 
            .groups="drop") %>%
  group_by(location, taxon) %>%
  mutate(taxon = str_replace(taxon,
                             "(.*)_unclassified", "Unclassified *\\1*"),
         taxon = str_replace(taxon,
                             "^(\\S*)$", "*\\1*"),
         taxon = str_replace_all(taxon, 
                                 "_", " ")) #gets ride of _ in taxon names

genus_rel_abund %>%
  group_by(taxon) %>%
  summarize(rel_abund=sum(rel_abund)) %>%
  arrange(rel_abund) %>%
  view()

top_genus <- genus_rel_abund %>%
  group_by(taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  arrange(rel_abund) %>%
  top_n(12, rel_abund) %>%
  pull(taxon)
view(top_genus)

genus_rel_abund %>%
  filter(taxon %in% top_genus) %>%
  mutate(taxon=factor(taxon, levels=top_genus)) %>%
  mutate(rel_abund=100 * rel_abund)

genus_rel_abund %>%
  filter(taxon %in% top_genus) %>%
  mutate(taxon=factor(taxon, levels=top_genus)) %>%
  mutate(rel_abund= 100 * (rel_abund + 1/20000)) %>%
  ggplot(aes(x=sample, y=taxon, fill=rel_abund)) +
  geom_tile(color="black") +
  scale_fill_gradient(low="white", high="red",
                      name=("Relative<br>Abundance (%)"),
                      limits=c(0,50), 
                      oob=squish, #part of "scales" library
                      breaks = c(0, 10, 20, 30, 40, 50), labels = c("0", "10", "20", "30", "40", ">50")) + 
  coord_fixed() +
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
  scale_y_discrete(breaks=c("Unclassified Nectriaceae", "*Aspergillus*", "*Alternaria*", "Unclassified Fungi", "Unclassified Ascomycota", "Unclassified Sordariaceae", "*Phoma*",
                            "*Talaromyces*", "*Penicillium*", "*Stagonosporopsis*", "*Macrophomina*", "*Trichoderma*"),
                   labels=c("**Unclassified Nectriaceae (BF)***", "*Aspergillus*", "***Alternaria* (SF)***", "Unclassified Fungi", "**Unclassified Ascomycota (BF)***", 
                            "**Unclassified Sordariaceae (S)***", "*Phoma*", "***Talaromyces* (S)***", "***Penicillium* (S)***", "***Stagonosporopsis* (BF)***", "*Macrophomina*", 
                            "***Trichoderma* (BF)***")) +
  labs(title="Relative Abundance of the Top 12 Genera",
       x=NULL,
       y=NULL) +
  theme_classic() +
  theme(axis.line = element_blank(), #gets rid of the axis lines
        axis.ticks = element_blank(), #gets rid of the axis tick marks
        axis.text.y = element_markdown(size=8),
        axis.text.x = element_markdown(size=8),
        #legend.title.align = 0.5,
        legend.title = element_markdown(size=8),
        legend.text = element_text(size=7.5),
        legend.key.height = unit(11, "pt"),
        plot.title = element_text(hjust=0.5))

ggsave("heatmap.ITS.genus.png", width=17, height=7)
