library(tidyverse)
library(broom)
library(ggtext)
library(readxl)
library(glue)

set.seed(19980609)

shared <- read_tsv("final.opti_mcc.0.03.subsample.mar.24.fungal.shared",
                   col_types = cols(Group = col_character(),
                                    .default = col_double())) %>%
  rename_all(tolower) %>%
  select(group, starts_with("otu")) %>%
  pivot_longer(-group, names_to="otu", values_to="count") %>%
  mutate(otu = str_replace(string=otu,
                           pattern="otu",
                           replacement = "Otu"))

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
         genus = str_replace(string=genus,
                             pattern="\\*unclassified_(.*)\\*",
                             replacement="Unclassified \\1"),
         genus = str_replace_all(genus, "_", " "),
         taxon = glue("{genus}<br>({pretty_otu})")) %>%
  select(otu, genus) #can swap in any taxonomic level

metadata <- read_excel("mar.24.fungal.metadata.xlsx") %>%
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
  filter(combo == "BF live" | combo == "Stiles live")

summary_composite <- composite %>%
  subset(rel_abund !="0") %>%
  mutate(rel_abund = 100 * (rel_abund + 1/20000)) %>%
  #subset(combo == "BF live") %>%
  subset(combo == "Stiles live") %>%
  group_by(genus) %>%
  summarize(mean=mean(rel_abund),
            se = sd(rel_abund) / sqrt(n()),
            max = mean+se,
            min = mean-se,
            N = n(),
            .groups = "drop")

sig_taxa <- composite %>%
  nest(data = -genus) %>%
  mutate(test = map(.x=data, ~pairwise.wilcox.test(x=.x$rel_abund, g=.x$treatment, p.adjust.method = "BH") %>% tidy())) %>%
  unnest(test) %>% 
  filter(p.value < 0.05) %>%
  #select(genus, p.value) %>%
  ungroup()

#for order
composite %>%
  inner_join(sig_taxa, by="order") %>%
  mutate(rel_abund = 100 * (rel_abund + 1/20000),
         #combo = factor(location, levels = c("BF live", "BF dead", "BF none", "Stiles live", "Stiles dead", "Stiles none"))) %>%
         location = factor(location, levels = c("BF", "Stiles"))) %>%
  ggplot(aes(x=rel_abund, y=order, color=treatment, fill=treatment)) +
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
                     breaks = c("live", "dead", "none"),
                     values = c('dodgerblue', 'red', "black"),
                     labels = c("Live", "Heat-killed", "Bulk soil")) +
  scale_fill_manual(NULL,
                    breaks = c("live", "dead", "none"),
                    values = c('dodgerblue', 'red', "black"),
                    labels = c("Live", "Heat-killed", "Bulk soil")) +
  labs(x= "Relative abundance", y=NULL,
       title="Significantly Different Orders (16S rRNA Gene)") +
  theme_classic() +
  theme(axis.text.y = element_markdown(size=9),
        legend.text = element_markdown(size=10),
        legend.key.size = unit(15, "pt"),
        plot.title = element_text(hjust=0.5))

ggsave("sig.diff.orders.ITS.trt.png", width=11, height=8)

#for family
composite %>%
  inner_join(sig_taxa, by="family") %>%
  mutate(rel_abund = 100 * (rel_abund + 1/20000),
         #combo = factor(location, levels = c("BF live", "BF dead", "BF none", "Stiles live", "Stiles dead", "Stiles none"))) %>%
         location = factor(location, levels = c("BF", "Stiles"))) %>%
  ggplot(aes(x=rel_abund, y=family, color=treatment, fill=treatment)) +
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
                     breaks = c("live", "dead", "none"),
                     values = c('dodgerblue', 'red', "black"),
                     labels = c("Live", "Heat-killed", "Bulk soil")) +
  scale_fill_manual(NULL,
                    breaks = c("live", "dead", "none"),
                    values = c('dodgerblue', 'red', "black"),
                    labels = c("Live", "Heat-killed", "Bulk soil")) +
  labs(x= "Relative abundance", y=NULL,
       title="Significantly Different Families (16S rRNA Gene)") +
  theme_classic() +
  theme(axis.text.y = element_markdown(size=9),
        legend.text = element_markdown(size=10),
        legend.key.size = unit(15, "pt"),
        plot.title = element_text(hjust=0.5))

ggsave("sig.diff.families.ITS.trt.BF.png", width=11, height=8)

#for genera
composite %>%
  inner_join(sig_taxa, by="genus") %>%
  mutate(rel_abund = 100 * (rel_abund + 1/20000),
         #combo = factor(location, levels = c("BF live", "BF dead", "BF none", "Stiles live", "Stiles dead", "Stiles none"))) %>%
         location = factor(location, levels = c("BF", "Stiles"))) %>%
  ggplot(aes(x=rel_abund, y=genus, color=treatment, fill=treatment)) +
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
                     breaks = c("live", "dead", "none"),
                     values = c('dodgerblue', 'red', "black"),
                     labels = c("Live", "Heat-killed", "Bulk soil")) +
  scale_fill_manual(NULL,
                    breaks = c("live", "dead", "none"),
                    values = c('dodgerblue', 'red', "black"),
                    labels = c("Live", "Heat-killed", "Bulk soil")) +
  labs(x= "Relative abundance", y=NULL,
       title="Significantly Different Genera (16S rRNA Gene)<br>Bottom Farm") +
  theme_classic() +
  theme(axis.text.y = element_markdown(size=9),
        legend.text = element_markdown(size=10),
        legend.key.size = unit(15, "pt"),
        plot.title = element_markdown(hjust=0.5))

ggsave("sig.diff.genera.ITS.trt.BF.png", width=11, height=8)
