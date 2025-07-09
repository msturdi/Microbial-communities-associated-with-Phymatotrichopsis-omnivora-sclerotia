library(tidyverse)
library(readxl)
library(glue)
library(ggtext)

taxonomy <- read_tsv(file="final.opti_mcc.0.03.cons.taxonomy") %>%
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
                             replacement="Unclassified<br>*\\1*"),
         taxon = glue("{genus}<br>({pretty_otu})")) %>%
  select(otu, taxon)

metadata <- read_excel("16S.2023.metadata.xlsx")

shared_file <- read_tsv("final.opti_mcc.0.03.subsample.shared")

shared_design <- inner_join(shared_file, metadata, by=c("Group"="sample"))

#creating new directory in my 16S folder
dir.create("processed_data/")

run_lefse <- function(x, y, tag){
  x_y <- shared_design %>%
    filter(trt == x | trt == y)
  
  x_y %>%
    select(-trt) %>%
    write_tsv(glue("processed_data/sclerotia.{tag}.shared"))
  
  x_y %>%
    select(Group, trt) %>%
    write_tsv(glue("processed_data/sclerotia.{tag}.design"))
  
  command <- glue('mothur/mothur "#lefse(shared=sclerotia.{tag}.shared, design=sclerotia.{tag}.design, inputdir=processed_data)"')
  
  system(command)
  
  return(glue("processed_data/sclerotia.{tag}.0.03.lefse_summary"))
}

live_dead <- run_lefse("live", "dead", "live_dead")
live_none <- run_lefse("live", "none", "live_none")
dead_none <- run_lefse("dead", "none", "dead_none")

#for live-dead comparison
read_tsv("processed_data/sclerotia.live_dead.0.03.lefse_summary")

read_tsv(live_dead) %>%
  drop_na(LDA) %>%
  filter(LDA > 2.5) %>% #filter out the otus that have a smaller effect size
  inner_join(., taxonomy, by=c("OTU" = "otu")) %>%
  mutate(LDA = if_else(Class == "live", -1 * LDA, LDA),
         taxon = fct_reorder(taxon, LDA)) %>%
  ggplot(aes(x=LDA, y=taxon, fill=Class)) +
  geom_col() +
  labs(y=NULL, x="LDA Score(log 10)",
       title="16S rRNA Gene LDA Plot (LDA > 2.5) Live vs. Dead") +
  scale_x_continuous(limits = c(-5, 5), breaks = seq(-5, 5, by=2.5)) +
  scale_fill_manual(name=NULL, 
                    breaks = c("dead", "live"),
                    labels = c("Heat-killed<br>Sclerotia", "Live Sclerotia"),
                    values = c('#d7191c', '#2c7bb6')) +
  theme_classic() +
  theme(axis.text.y = element_markdown(),
        legend.text = element_markdown())

ggsave("lefse.live.dead.16S.png", width=7, height=7)

#for live-none comparison
read_tsv("processed_data/sclerotia.live_none.0.03.lefse_summary")

read_tsv(live_none) %>%
  drop_na(LDA) %>%
  filter(LDA > 2.5) %>% #filter out the otus that have a smaller effect size
  #subset(OTU !="Otu0049") %>% #this is the otu that is NA in the figure
  inner_join(., taxonomy, by=c("OTU" = "otu")) %>%
  mutate(LDA = if_else(Class == "live", -1 * LDA, LDA),
         taxon = fct_reorder(taxon, LDA)) %>%
  ggplot(aes(x=LDA, y=taxon, fill=Class)) +
  geom_col() +
  labs(y=NULL, x="LDA Score(log 10)",
       title="16S rRNA Gene LDA Plot (LDA > 2.5)") +
  scale_x_continuous(limits = c(-5, 5), breaks = seq(-5, 5, by=2.5)) +
  scale_fill_manual(name=NULL, 
                    breaks = c("live", "none"),
                    labels = c("Live Sclerotia", "Bulk Soil"),
                    values = c('#2c7bb6', "black")) +
  theme_classic() +
  theme(axis.text.y = element_markdown(size = 7),
        legend.text = element_markdown())

ggsave("lefse.live.none.16S.png", width=7, height=12)

#for dead-none comparison
read_tsv("processed_data/sclerotia.dead_none.0.03.lefse_summary")

read_tsv(dead_none) %>%
  drop_na(LDA) %>%
  filter(LDA > 2.5) %>% #filter out the otus that have a smaller effect size
  inner_join(., taxonomy, by=c("OTU" = "otu")) %>%
  mutate(LDA = if_else(Class == "dead", -1 * LDA, LDA),
         taxon = fct_reorder(taxon, LDA)) %>%
  ggplot(aes(x=LDA, y=taxon, fill=Class)) +
  geom_col() +
  labs(y=NULL, x="LDA Score(log 10)",
       title="16S rRNA Gene LDA Plot (LDA > 2.5)") +
  scale_x_continuous(limits = c(-5, 5), breaks = seq(-5, 5, by=2.5)) +
  scale_fill_manual(name=NULL, 
                    breaks = c("dead", "none"),
                    labels = c("Heat-killed<br>Sclerotia", "Bulk Soil"),
                    values = c('#d7191c', "black")) +
  theme_classic() +
  theme(axis.text.y = element_markdown(size = 6),
        legend.text = element_markdown())

ggsave("lefse.dead.none.16S.png", width=7, height=13)
