library(tidyverse)
library(broom)
library(ggtext)
library(readxl)
library(glue)

set.seed(19980609)

shared <- read_tsv("final.opti_mcc.0.03.subsample.mar.24.bact.shared",
                   col_types = cols(Group = col_character(),
                                    .default = col_double())) %>%
  rename_all(tolower) %>%
  select(group, starts_with("otu")) %>%
  pivot_longer(-group, names_to="otu", values_to="count") %>%
  mutate(otu = str_replace(string=otu,
                           pattern="otu",
                           replacement = "Otu"))

##
#genus-level
##
taxonomy <- read_tsv(file="final.opti_mcc.0.03.cons.mar.24.bact.taxonomy") %>%
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
         taxon = glue("{genus}<br>({pretty_otu})"))%>%
  select(otu, genus)

##
#order-level
##
taxonomy <- read_tsv(file="final.opti_mcc.0.03.cons.mar.24.bact.taxonomy") %>%
  select("OTU", "Taxonomy") %>%
  rename_all(tolower) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="\\(\\d*\\)", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern=";$", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern='["]', replacement="")) %>%
  separate(taxonomy, into=c("kingdom", "phylum", "class", "order", "family", "genus"), 
           sep=";") %>%
  #mutation code for family or genus level below
  mutate(pretty_otu = str_replace(string=otu,
                                  pattern="tu0*",
                                  replacement = "TU "),
         order=str_replace(order, "unclassified_(.*)", "Unclassified \\1"),
         order=str_replace(order, "(.*)_unclassified", "Unclassified \\1"),
         order=str_replace_all(order, "_", " "))%>%
  select(otu, order) #can swap in any taxonomic level, make sure to do it above too

##
#class-level
##
taxonomy <- read_tsv(file="final.opti_mcc.0.03.cons.mar.24.bact.taxonomy") %>%
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
         class=str_replace(class, "unclassified_(.*)", "Unclassified \\1"),
         class=str_replace(class, "(.*)_unclassified", "Unclassified \\1"),
         class=str_replace_all(class, "_", " "))%>%
  select(otu, class)

metadata <- read_excel("mar.24.bact.metadata.xlsx") %>%
  rename(group = sample) %>%
  mutate(combo = glue("{location} {treatment}"))

composite <- inner_join(shared, taxonomy, by="otu") %>%
  group_by(group, genus) %>%
  summarize(count = sum(count), .groups="drop") %>%
  group_by(group) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  select(-count) %>%
  inner_join(., metadata, by="group") %>%
  #filter(location == "BF")
  filter(combo == "BF live" | combo == "Stiles live") 

summary_composite <- composite %>%
  subset(rel_abund !="0") %>%
  mutate(rel_abund = 100 * (rel_abund + 1/20000)) %>%
  subset(combo == "BF live") %>%
  #subset(combo == "Stiles live") %>%
  group_by(genus) %>%
  summarize(mean=mean(rel_abund),
            se = sd(rel_abund) / sqrt(n()),
            max = mean+se,
            min = mean-se,
            N = n(),
            .groups = "drop")

sig_taxa <- composite %>%
  nest(data = -order) %>%
  mutate(test = map(.x=data, ~pairwise.wilcox.test(x=.x$rel_abund, g=.x$combo, p.adjust.method = "BH") %>% tidy())) %>%
  unnest(test) %>% 
  filter(p.value < 0.05) %>%
  #select(genus, p.value) %>%
  ungroup()

#for phyla
composite %>%
  inner_join(sig_taxa, by="phylum") %>%
  mutate(rel_abund = 100 * (rel_abund + 1/20000),
         location = factor(location, levels = c("BF", "Stiles"))) %>%
  ggplot(aes(x=rel_abund, y=phylum, color=location, fill=location)) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8,
                                              jitter.width = 0.5),
              shape=21,
              size=1.4) +
  stat_summary(fun.data = median_hilow, fun.args = list(conf.int=0.5),
               geom="pointrange",
               position = position_dodge(width=0.1),
               color="darkgrey", show.legend = FALSE,
               size = 0.1) +
  scale_x_log10() +
  scale_color_manual(NULL,
                     breaks = c("BF", "Stiles"),
                     values = c('orange2', 'midnightblue'),
                     labels = c("Bottom Farm", "Stiles Farm")) +
  scale_fill_manual(NULL,
                    breaks = c("BF", "Stiles"),
                    values = c('orange2', 'midnightblue'),
                    labels = c("Bottom Farm", "Stiles Farm"))  +
  labs(x= "Relative abundance", y=NULL,
       title="Significantly Different Phyla (16S rRNA Gene)") +
  theme_classic() +
  theme(axis.text.y = element_markdown(size=9),
        legend.text = element_markdown(size=10),
        legend.key.size = unit(15, "pt"),
        plot.title = element_text(hjust=0.5))

ggsave("sig.diff.phyla.16S.location.png", width=11, height=8)

#for class
composite %>%
  inner_join(sig_taxa, by="class") %>%
  mutate(rel_abund = 100 * (rel_abund + 1/20000),
         location = factor(location, levels = c("BF", "Stiles"))) %>%
  ggplot(aes(x=rel_abund, y=class, color=location, fill=location)) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8,
                                              jitter.width = 0.5),
              shape=21,
              size=1.4) +
  stat_summary(fun.data = median_hilow, fun.args = list(conf.int=0.5),
               geom="pointrange",
               position = position_dodge(width=0.1),
               color="darkgrey", show.legend = FALSE,
               size = 0.1) +
  scale_x_log10() +
  scale_color_manual(NULL,
                     breaks = c("BF", "Stiles"),
                     values = c('orange2', 'midnightblue'),
                     labels = c("Bottom Farm", "Stiles Farm")) +
  scale_fill_manual(NULL,
                    breaks = c("BF", "Stiles"),
                    values = c('orange2', 'midnightblue'),
                    labels = c("Bottom Farm", "Stiles Farm"))  +
  labs(x= "Relative abundance", y=NULL,
       title="Significantly Different Classes (16S rRNA Gene)") +
  theme_classic() +
  theme(axis.text.y = element_markdown(size=9),
        legend.text = element_markdown(size=10),
        legend.key.size = unit(15, "pt"),
        plot.title = element_text(hjust=0.5))

ggsave("sig.diff.classes.16S.location.png", width=11, height=8)

#for order
composite %>%
  inner_join(sig_taxa, by="order") %>%
  mutate(rel_abund = 100 * (rel_abund + 1/20000),
         location = factor(location, levels = c("BF", "Stiles"))) %>%
  ggplot(aes(x=rel_abund, y=order, color=location, fill=location)) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8,
                                              jitter.width = 0.5),
              shape=21,
              size=1.4) +
  stat_summary(fun.data = median_hilow, fun.args = list(conf.int=0.5),
               geom="pointrange",
               position = position_dodge(width=0.1),
               color="darkgrey", show.legend = FALSE,
               size = 0.1) +
  scale_x_log10() +
  scale_color_manual(NULL,
                     breaks = c("BF", "Stiles"),
                     values = c('orange2', 'midnightblue'),
                     labels = c("Bottom Farm", "Stiles Farm")) +
  scale_fill_manual(NULL,
                    breaks = c("BF", "Stiles"),
                    values = c('orange2', 'midnightblue'),
                    labels = c("Bottom Farm", "Stiles Farm"))  +
  labs(x= "Relative abundance", y=NULL,
       title="Significantly Different Orders (Bacterial)<br>Bottom Farm vs. Stiles Farm") +
  theme_light() +
  theme(axis.text.y = element_markdown(size=9),
        legend.text = element_markdown(size=10),
        legend.key.size = unit(15, "pt"),
        legend.position = "bottom",
        plot.title = element_markdown(hjust=0.5))

ggsave("sig.diff.orders.16S.location.png", width=12, height=12)

#for BF live vs. Stiles live
composite %>%
  inner_join(sig_taxa, by="order") %>%
  mutate(rel_abund = 100 * (rel_abund + 1/20000),
         combo = factor(combo, levels = c("BF live", "Stiles live"))) %>%
         #location = factor(location, levels = c("BF", "Stiles"))) %>%
  ggplot(aes(x=rel_abund, y=order, color=combo, fill=combo)) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8,
                                              jitter.width = 0.5),
              shape=21,
              size=1.4) +
  stat_summary(fun.data = median_hilow, fun.args = list(conf.int=0.5),
               geom="pointrange",
               position = position_dodge(width=0.1),
               color="darkgrey", show.legend = FALSE,
               size = 0.1) +
  scale_x_log10() +
  scale_color_manual(NULL,
                     breaks = c("BF live", "Stiles live"),
                     values = c("#E97451", "#7393B3"),
                     labels = c("Bottom Farm Live", "Stiles Farm Live")) +
  scale_fill_manual(NULL,
                    breaks = c("BF live", "Stiles live"),
                    values = c("#E97451", "#7393B3"),
                    labels = c("Bottom Farm Live", "Stiles Farm Live"))  +
  labs(x= "Relative abundance", y=NULL,
       title="Significantly Different Orders (Bacterial)<br>Bottom Farm Live vs. Stiles Farm Live") +
  theme_light() +
  theme(axis.text.y = element_markdown(size=9),
        legend.text = element_markdown(size=10),
        legend.key.size = unit(15, "pt"),
        legend.position = "bottom",
        plot.title = element_markdown(hjust=0.5))

ggsave("sig.diff.orders.16S.blive.slive.png", width=12, height=12)
