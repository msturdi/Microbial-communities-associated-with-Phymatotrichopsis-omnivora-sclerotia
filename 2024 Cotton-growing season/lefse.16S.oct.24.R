library(tidyverse)
library(readxl)
library(glue)
library(ggtext)
library(dplyr)

#creating new directory in my folder, only need to do first time with new dataset
dir.create("processed_data/")

#read in taxonomy and shared files first - metadata and shared_design change based on variables in chosen analysis
taxonomy <- read_tsv(file="final.opti_mcc.0.03.cons.oct.24.bacterial.taxonomy") %>%
  select("OTU", "Taxonomy") %>%
  rename_all(tolower) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="\\(\\d*\\)", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern=";$", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern='["]', replacement="")) %>%
  separate(taxonomy, into=c("kingdom", "phylum", "class", "order", "family", "genus"), 
           sep=";") %>%
  mutate(pretty_otu = str_replace(string=otu,
                                  pattern="tu0*",
                                  replacement = "TU "),
         genus = str_replace(string=genus,
                             pattern="(.*)",
                             replacement="*\\1*"),
         genus = str_replace(string=genus,
                             pattern="\\*(.*)_unclassified\\*",
                             replacement="Unclassified \\1"),
         taxon = glue("{genus} ({pretty_otu})"),
         taxon=str_replace_all(taxon, "_", " ")) %>%
  select(otu, taxon)

shared_file <- read_tsv("final.opti_mcc.0.03.subsample.oct.24.bacterial.shared")

#for location comparison

metadata <- read_excel("oct.24.16S.metadata.xlsx") %>%
  select(sample, location)

#shared_design <- inner_join(shared_file, metadata, by=c("Group"="sample"))
shared_design <- inner_join(shared_filtered, metadata, by=c("Group"="sample"))

run_lefse_location <- function(x, y, tag){
  x_y <- shared_design %>%
    filter(location == x | location == y)
  
  x_y %>%
    select(-location) %>%
    write_tsv(glue("processed_data/location.{tag}.shared"))
  
  x_y %>%
    select(Group, location) %>%
    write_tsv(glue("processed_data/location.{tag}.design"))
  
  command <- glue('mothur/mothur "#lefse(shared=location.{tag}.shared, design=location.{tag}.design, inputdir=processed_data)"')
  
  system(command)
  
  return(glue("processed_data/location.{tag}.0.03.lefse_summary"))
}

bf_stiles <- run_lefse_location("BF", "Stiles", "bf_stiles")

#for BF-Stiles comparison
read_tsv(bf_stiles) %>%
  drop_na(LDA) %>%
  filter(LDA > 2.5) %>%
  inner_join(., taxonomy, by=c("OTU" = "otu")) %>%
  filter(!grepl("Unclassified", taxon)) %>%
  ggplot(aes(x=LDA, y=taxon, fill=Class)) +
  geom_col() +
  labs(y=NULL, x="LDA Score(log 10)",
       title="LDA Plot (LDA > 2.5) In-season Bacterial<br>Bottom Farm vs. Stiles Farm") +
  scale_x_continuous(expand = c(0,0), limits = c(0, 4), breaks = seq(0, 4, by=2)) +
  scale_fill_manual(name=NULL, 
                    breaks = c("BF", "Stiles"),
                    labels = c("Bottom Farm", "Stiles Farm"),
                    values = c('orange2', 'midnightblue')) +
  theme_classic() +
  theme(axis.text.y = element_markdown(),
        legend.text = element_markdown(),
        plot.title = element_markdown(hjust=0.5))

ggsave("lefse.BF.Stiles.16S.png", width=12, height=22)

#for treatment comparison

metadata <- read_excel("oct.24.16S.metadata.xlsx") %>%
  mutate(combo = glue("{location} {treatment}")) %>%
  mutate(combo = str_replace_all(string=combo, pattern="BF live", replacement="blive")) %>%
  mutate(combo = str_replace_all(string=combo, pattern="BF dead", replacement="bdead")) %>%
  mutate(combo = str_replace_all(string=combo, pattern="BF none", replacement="bnone")) %>%
  mutate(combo = str_replace_all(string=combo, pattern="Stiles live", replacement="slive")) %>%
  mutate(combo = str_replace_all(string=combo, pattern="Stiles dead", replacement="sdead")) %>%
  mutate(combo = str_replace_all(string=combo, pattern="Stiles none", replacement="snone")) %>%
  select(sample, combo)

shared_design <- inner_join(shared_filtered, metadata, by=c("Group"="sample"))

run_lefse_trt <- function(x, y, tag){
  x_y <- shared_design %>%
    filter(combo == x | combo == y)
  
  x_y %>%
    select(-combo) %>%
    write_tsv(glue("processed_data/combo.{tag}.shared"))
  
  x_y %>%
    select(Group, combo) %>%
    write_tsv(glue("processed_data/combo.{tag}.design"))
  
  command <- glue('mothur/mothur "#lefse(shared=combo.{tag}.shared, design=combo.{tag}.design, inputdir=processed_data)"')
  
  system(command)
  
  return(glue("processed_data/combo.{tag}.0.03.lefse_summary"))
}

blive_bnone <- run_lefse_trt("blive", "bnone", "blive_bnone")
bdead_bnone <- run_lefse_trt("bdead", "bnone", "bdead_bnone")
slive_snone <- run_lefse_trt("slive", "snone", "slive_snone")
sdead_snone <- run_lefse_trt("sdead", "snone", "sdead_snone")
blive_slive <- run_lefse_trt("blive", "slive", "blive_slive")
bdead_sdead <- run_lefse_trt("bdead", "sdead", "bdead_sdead")

#for BF live - BF none comparison
read_tsv(blive_bnone) %>%
  drop_na(LDA) %>%
  filter(LDA > 2.5) %>%
  inner_join(., taxonomy, by=c("OTU" = "otu")) %>%
  filter(!grepl("Unclassified", taxon)) %>%
  ggplot(aes(x=LDA, y=taxon, fill=Class)) +
  geom_col() +
  labs(y=NULL, x="LDA Score(log 10)",
       title="LDA Plot (LDA > 2.5) In-season Bacterial<br>Bottom Farm") +
  scale_x_continuous(expand = c(0,0), limits = c(0, 4.2), breaks = seq(0, 4.2, by=2)) +
  scale_fill_manual(name=NULL, 
                    breaks = c("blive", "bnone"),
                    labels = c("Live Sclerotia", "Bulk Soil"),
                    values = c('#2c7bb6', "black")) +
  theme_classic() +
  theme(axis.text.y = element_markdown(),
        legend.text = element_markdown(),
        plot.title = element_markdown(hjust=0.5))

ggsave("lefse.blive.bnone.16S.png", width=10, height=8)

#for BF dead - BF none comparison
read_tsv(bdead_bnone) %>%
  drop_na(LDA) %>%
  filter(LDA > 2.5) %>%
  inner_join(., taxonomy, by=c("OTU" = "otu")) %>%
  filter(!grepl("Unclassified", taxon)) %>%
  ggplot(aes(x=LDA, y=taxon, fill=Class)) +
  geom_col() +
  labs(y=NULL, x="LDA Score(log 10)",
       title="LDA Plot (LDA > 2.5) In-season Bacterial<br>Bottom Farm") +
  scale_x_continuous(expand = c(0,0), limits = c(0, 4), breaks = seq(0, 4, by=2)) +
  scale_fill_manual(name=NULL, 
                    breaks = c("bdead", "bnone"),
                    labels = c("Heat-killed<br>Sclerotia", "Bulk Soil"),
                    values = c('#d7191c', "black")) +
  theme_classic() +
  theme(axis.text.y = element_markdown(),
        legend.text = element_markdown(),
        plot.title = element_markdown(hjust=0.5))

ggsave("lefse.bdead.bnone.16S.png", width=9, height=7)

#for Stiles live - Stiles none comparison
read_tsv(slive_snone) %>%
  drop_na(LDA) %>%
  filter(LDA > 2.5) %>%
  inner_join(., taxonomy, by=c("OTU" = "otu")) %>%
  filter(!grepl("Unclassified", taxon)) %>%
  ggplot(aes(x=LDA, y=taxon, fill=Class)) +
  geom_col() +
  labs(y=NULL, x="LDA Score(log 10)",
       title="LDA Plot (LDA > 2.5) In-season Bacterial<br>Stiles Farm") +
  scale_x_continuous(expand = c(0,0), limits = c(0, 5), breaks = seq(0, 5, by=2.5)) +
  scale_fill_manual(name=NULL, 
                    breaks = c("slive", "snone"),
                    labels = c("Live Sclerotia", "Bulk Soil"),
                    values = c('#2c7bb6', "black")) +
  theme_classic() +
  theme(axis.text.y = element_markdown(),
        legend.text = element_markdown(),
        plot.title = element_markdown(hjust=0.5))

ggsave("lefse.slive.snone.16S.png", width=9, height=6)

#for Stiles dead - Stiles none comparison
read_tsv(sdead_snone) %>%
  drop_na(LDA) %>%
  filter(LDA > 2.5) %>%
  inner_join(., taxonomy, by=c("OTU" = "otu")) %>%
  filter(!grepl("Unclassified", taxon)) %>%
  ggplot(aes(x=LDA, y=taxon, fill=Class)) +
  geom_col() +
  labs(y=NULL, x="LDA Score(log 10)",
       title="LDA Plot (LDA > 2.5) Off-season Bacterial<br>Stiles Farm") +
  scale_x_continuous(expand = c(0,0), limits = c(0, 5), breaks = seq(0, 5, by=2.5)) +
  scale_fill_manual(name=NULL, 
                    breaks = c("sdead", "snone"),
                    labels = c("Heat-killed<br>Sclerotia", "Bulk Soil"),
                    values = c('#d7191c', "black")) +
  theme_classic() +
  theme(axis.text.y = element_markdown(),
        legend.text = element_markdown(),
        plot.title = element_markdown(hjust=0.5))

ggsave("lefse.sdead.snone.16S.png", width=9, height=7)

#for BF live - Stiles live comparison
read_tsv(blive_slive) %>%
  drop_na(LDA) %>%
  filter(LDA > 2.5) %>%
  inner_join(., taxonomy, by=c("OTU" = "otu")) %>%
  filter(!grepl("Unclassified", taxon)) %>%
  ggplot(aes(x=LDA, y=taxon, fill=Class)) +
  geom_col() +
  labs(y=NULL, x="LDA Score(log 10)",
       title="LDA Plot (LDA > 2.5) In-season Bacterial<br>Bottom Farm Live vs. Stiles Farm Live") +
  scale_x_continuous(expand = c(0,0), limits = c(0, 4.2), breaks = seq(0, 4.2, by=2)) +
  scale_fill_manual(name=NULL, 
                    breaks = c("blive", "slive"),
                    labels = c("Bottom Farm<br>live", "Stiles<br>live"),
                    values = c("#E97451", "#7393B3")) +
  theme_classic() +
  theme(axis.text.y = element_markdown(),
        legend.text = element_markdown(),
        plot.title = element_markdown(hjust=0.5))

ggsave("lefse.blive.slive.16S.png", width=12, height=22)

#for BF dead - Stiles dead comparison
read_tsv(bdead_sdead) %>%
  drop_na(LDA) %>%
  filter(LDA > 2.5) %>%
  inner_join(., taxonomy, by=c("OTU" = "otu")) %>%
  filter(!grepl("Unclassified", taxon)) %>%
  ggplot(aes(x=LDA, y=taxon, fill=Class)) +
  geom_col() +
  labs(y=NULL, x="LDA Score(log 10)",
       title="LDA Plot (LDA > 2.5) In-season Bacterial<br>Bottom Farm Dead vs. Stiles Farm Dead") +
  scale_x_continuous(expand = c(0,0), limits = c(0, 4), breaks = seq(0, 4, by=2)) +
  scale_fill_manual(name=NULL, 
                    breaks = c("bdead", "sdead"),
                    labels = c("Bottom Farm<br>dead", "Stiles<br>dead"),
                    values = c("#FF5349", "#997570")) +
  theme_classic() +
  theme(axis.text.y = element_markdown(),
        legend.text = element_markdown(),
        plot.title = element_markdown(hjust=0.5))

ggsave("lefse.bdead.sdead.16S.png", width=12, height=22)
