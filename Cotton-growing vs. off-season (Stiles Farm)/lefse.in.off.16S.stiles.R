library(tidyverse)
library(readxl)
library(glue)
library(ggtext)
library(dplyr)

#creating new directory in my folder, only need to do first time with new dataset
dir.create("processed_data/")

#read in taxonomy and shared files first - metadata and shared_design change based on variables in chosen analysis
taxonomy <- read_tsv(file="final.opti_mcc.0.03.cons.bacterial.stiles.taxonomy") %>%
  select("OTU", "Taxonomy") %>%
  rename_all(tolower) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="\\(\\d*\\)", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern=";$", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern='["]', replacement="")) %>%
  separate(taxonomy, into=c("kingdom", "phylum", "class", "order", "family", "genus"), 
           sep=";") %>%
  #subset(genus !="NA") %>%
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

shared_file <- read_tsv("final.opti_mcc.0.03.subsample.bacterial.stiles.shared")

#for season comparison

metadata <- read_excel("metadata.16S.in.off.stiles.xlsx") %>%
  select(sample, season) #lefse only wants 2 columns, the sample and then the variable you want to investigate
#can switch it out for whatever column I want to look at

shared_design <- inner_join(shared_file, metadata, by=c("Group"="sample"))

run_lefse_season <- function(x, y, tag){
  x_y <- shared_design %>%
    filter(season == x | season == y)
  
  x_y %>%
    select(-season) %>%
    write_tsv(glue("processed_data/season.{tag}.shared"))
  
  x_y %>%
    select(Group, season) %>%
    write_tsv(glue("processed_data/season.{tag}.design"))
  
  command <- glue('mothur/mothur "#lefse(shared=season.{tag}.shared, design=season.{tag}.design, inputdir=processed_data)"')
  
  system(command)
  
  return(glue("processed_data/season.{tag}.0.03.lefse_summary"))
}

in_off <- run_lefse_season("in", "off", "in_off")

#for in-off comparison
read_tsv(in_off) %>%
  drop_na(LDA) %>%
  filter(LDA > 2.5) %>% #filter out the otus that have a smaller effect size
  inner_join(., taxonomy, by=c("OTU" = "otu")) %>%
  filter(!grepl("Unclassified", taxon)) %>%
  ggplot(aes(x=LDA, y=taxon, fill=Class)) +
  geom_col() +
  labs(y=NULL, x="LDA Score(log 10)",
       title="LDA Plot (LDA > 2.5) In- vs. Off-season Bacterial<br>Stiles Farm") +
  scale_x_continuous(expand = c(0,0), limits = c(0, 4), breaks = seq(0, 4, by=2)) +
  scale_fill_manual(name=NULL, 
                    breaks = c("in", "off"),
                    labels = c("In-season", "Off-season"),
                    values = c('darkgreen', 'dodgerblue3')) +
  theme_classic() +
  theme(axis.text.y = element_markdown(),
        legend.text = element_markdown(),
        plot.title = element_markdown(hjust=0.5))

ggsave("lefse.in.off.16S.season.stiles.png", width=12, height=22)


#for treatment comparison

metadata <- read_excel("metadata.16S.in.off.stiles.xlsx") %>%
  mutate(combo = glue("{season} {treatment}")) %>%
  mutate(combo = str_replace_all(string=combo, pattern="in live", replacement="ilive")) %>%
  mutate(combo = str_replace_all(string=combo, pattern="in dead", replacement="idead")) %>%
  mutate(combo = str_replace_all(string=combo, pattern="in none", replacement="inone")) %>%
  mutate(combo = str_replace_all(string=combo, pattern="off live", replacement="olive")) %>%
  mutate(combo = str_replace_all(string=combo, pattern="off dead", replacement="odead")) %>%
  mutate(combo = str_replace_all(string=combo, pattern="off none", replacement="onone")) %>%
  select(sample, combo) #lefse only wants 2 columns, the sample and then the variable you want to investigate
#can switch it out for whatever column I want to look at

shared_design <- inner_join(shared_file, metadata, by=c("Group"="sample"))

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

ilive_inone <- run_lefse_trt("ilive", "inone", "ilive_inone")
idead_inone <- run_lefse_trt("idead", "inone", "idead_inone")
olive_onone <- run_lefse_trt("olive", "onone", "olive_onone")
odead_onone <- run_lefse_trt("odead", "onone", "odead_onone")
ilive_olive <- run_lefse_trt("ilive", "olive", "ilive_olive")
idead_odead <- run_lefse_trt("idead", "odead", "idead_odead")

#for in-season live - in-season none comparison
read_tsv(ilive_inone) %>%
  drop_na(LDA) %>%
  filter(LDA > 2.5) %>% #filter out the otus that have a smaller effect size
  inner_join(., taxonomy, by=c("OTU" = "otu")) %>%
  filter(!grepl("Unclassified", taxon)) %>%
  ggplot(aes(x=LDA, y=taxon, fill=Class)) +
  geom_col() +
  labs(y=NULL, x="LDA Score(log 10)",
       title="LDA Plot (LDA > 2.5) In-season Bacterial<br>Stiles Farm") +
  scale_x_continuous(expand = c(0,0), limits = c(0, 4.2), breaks = seq(0, 4.2, by=2)) +
  scale_fill_manual(name=NULL, 
                    breaks = c("ilive", "inone"),
                    labels = c("In-season<br>Live Sclerotia", "In-season<br>Bulk Soil"),
                    values = c('#2c7bb6', "black")) +
  theme_classic() +
  theme(axis.text.y = element_markdown(),
        legend.text = element_markdown(),
        plot.title = element_markdown(hjust=0.5))

ggsave("lefse.ilive.inone.16S.stiles.png", width=10, height=8)

#for in-season dead - in-season none comparison
read_tsv(idead_inone) %>%
  drop_na(LDA) %>%
  filter(LDA > 2.5) %>% #filter out the otus that have a smaller effect size
  inner_join(., taxonomy, by=c("OTU" = "otu")) %>%
  filter(!grepl("Unclassified", taxon)) %>%
  ggplot(aes(x=LDA, y=taxon, fill=Class)) +
  geom_col() +
  labs(y=NULL, x="LDA Score(log 10)",
       title="LDA Plot (LDA > 2.5) In-season Bacterial<br>Stiles Farm") +
  scale_x_continuous(expand = c(0,0), limits = c(0, 4), breaks = seq(0, 4, by=2)) +
  scale_fill_manual(name=NULL, 
                    breaks = c("idead", "inone"),
                    labels = c("In-season<br>Heat-killed Sclerotia", "In-season<br>Bulk Soil"),
                    values = c('#d7191c', "black")) +
  theme_classic() +
  theme(axis.text.y = element_markdown(),
        legend.text = element_markdown(),
        legend.position = "bottom",
        plot.title = element_markdown(hjust=0.5))

ggsave("lefse.idead.inone.16S.stiles.png", width=10, height=8)

#for off-season live - off-season none comparison
read_tsv(olive_onone) %>%
  drop_na(LDA) %>%
  filter(LDA > 2.5) %>% #filter out the otus that have a smaller effect size
  inner_join(., taxonomy, by=c("OTU" = "otu")) %>%
  filter(!grepl("Unclassified", taxon)) %>%
  ggplot(aes(x=LDA, y=taxon, fill=Class)) +
  geom_col() +
  labs(y=NULL, x="LDA Score(log 10)",
       title="LDA Plot (LDA > 2.5) Off-season Bacterial<br>Stiles Farm") +
  scale_x_continuous(expand = c(0,0), limits = c(0, 4), breaks = seq(0, 4, by=2)) +
  scale_fill_manual(name=NULL, 
                    breaks = c("olive", "onone"),
                    labels = c("Off-season<br>Live Sclerotia", "Off-season<br>Bulk Soil"),
                    values = c('#2c7bb6', "black")) +
  theme_classic() +
  theme(axis.text.y = element_markdown(),
        legend.text = element_markdown(),
        legend.position = "bottom",
        plot.title = element_markdown(hjust=0.5))

ggsave("lefse.olive.onone.16S.stiles.png", width=12, height=15)

#for off-season dead - off-season none comparison
read_tsv(odead_onone) %>%
  drop_na(LDA) %>%
  filter(LDA > 2.5) %>% #filter out the otus that have a smaller effect size
  inner_join(., taxonomy, by=c("OTU" = "otu")) %>%
  filter(!grepl("Unclassified", taxon)) %>%
  ggplot(aes(x=LDA, y=taxon, fill=Class)) +
  geom_col() +
  labs(y=NULL, x="LDA Score(log 10)",
       title="LDA Plot (LDA > 2.5) Off-season Bacterial<br>Stiles Farm") +
  scale_x_continuous(expand = c(0,0), limits = c(0, 4), breaks = seq(0, 4, by=2)) +
  scale_fill_manual(name=NULL, 
                    breaks = c("odead", "onone"),
                    labels = c("Off-season<br>Heat-killed Sclerotia", "Off-season<br>Bulk Soil"),
                    values = c('#d7191c', "black")) +
  theme_classic() +
  theme(axis.text.y = element_markdown(),
        legend.text = element_markdown(),
        plot.title = element_markdown(hjust=0.5))

ggsave("lefse.odead.onone.16S.stiles.png", width=12, height=12)

#for in-season live - off-season live comparison
read_tsv(ilive_olive) %>%
  drop_na(LDA) %>%
  filter(LDA > 2.5) %>% #filter out the otus that have a smaller effect size
  inner_join(., taxonomy, by=c("OTU" = "otu")) %>%
  filter(!grepl("Unclassified", taxon)) %>% #dplyr library, 'grepl' matches a pattern within a string, getting rid of all instances of a word in a column, no matter what comes after it
  ggplot(aes(x=LDA, y=taxon, fill=Class)) +
  geom_col() +
  labs(y=NULL, x="LDA Score(log 10)",
       title="LDA Plot (LDA > 2.5) In- vs. Off-Season<br>Live Sclerotia Stiles Farm") +
  scale_x_continuous(expand = c(0,0), limits = c(0, 4), breaks = seq(0, 4, by=2)) +
  scale_fill_manual(name=NULL, 
                    breaks = c("ilive", "olive"),
                    labels = c("In-season<br>Live Sclerotia", "Off-season<br>Live Sclerotia"),
                    values = c("#008080", "#800020")) +
  theme_classic() +
  theme(axis.text.y = element_markdown(),
        legend.text = element_markdown(),
        plot.title = element_markdown(hjust=0.5))

ggsave("lefse.ilive.olive.16S.stiles.png", width=12, height=22)

#for in-season dead - off-season dead comparison
read_tsv(idead_odead) %>%
  drop_na(LDA) %>%
  filter(LDA > 2.5) %>% #filter out the otus that have a smaller effect size
  inner_join(., taxonomy, by=c("OTU" = "otu")) %>%
  filter(!grepl("Unclassified", taxon)) %>% #dplyr library, 'grepl' matches a pattern within a string, getting rid of all instances of a word in a column, no matter what comes after it
  ggplot(aes(x=LDA, y=taxon, fill=Class)) +
  geom_col() +
  labs(y=NULL, x="LDA Score(log 10)",
       title="LDA Plot (LDA > 2.5) In- vs. Off-Season<br>Dead Sclerotia Stiles Farm") +
  scale_x_continuous(expand = c(0,0), limits = c(0, 4), breaks = seq(0, 4, by=2)) +
  scale_fill_manual(name=NULL, 
                    breaks = c("idead", "odead"),
                    labels = c("In-season<br>Dead Sclerotia", "Off-season<br>Dead Sclerotia"),
                    values = c("#759116", "#8f00ff")) +
  theme_classic() +
  theme(axis.text.y = element_markdown(),
        legend.text = element_markdown(),
        plot.title = element_markdown(hjust=0.5))

ggsave("lefse.idead.odead.16S.stiles.png", width=12, height=22)
