library(tidyverse)
library(broom)
library(ggtext)
library(readxl)
library(glue)

set.seed(19980609)

shared <- read_tsv("final.opti_mcc.0.03.subsample.fungal.BF.shared",
                   col_types = cols(Group = col_character(),
                                    .default = col_double())) %>%
  rename_all(tolower) %>%
  select(group, starts_with("otu")) %>%
  pivot_longer(-group, names_to="otu", values_to="count") %>%
  mutate(otu = str_replace(string=otu,
                           pattern="otu",
                           replacement = "Otu"))

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
           sep=";") %>%
  mutate(pretty_otu = str_replace(string=otu,
                                  pattern="tu0*",
                                  replacement = "TU "),
         order=str_replace(order, "unclassified_(.*)", "Unclassified \\1"),
         order=str_replace(order, "(.*)_unclassified", "Unclassified \\1"),
         order=str_replace_all(order, "_", " "))%>%
  select(otu, order) #can swap in any taxonomic level, make sure to do it above too

metadata <- read_excel("metadata.ITS.in.off.BF.xlsx") %>%
  rename(group = sample) %>%
  mutate(combo = glue("{season} {treatment}"))

composite <- inner_join(shared, taxonomy, by="otu") %>%
  group_by(group, order) %>%
  summarize(count = sum(count), .groups="drop") %>%
  group_by(group) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  select(-count) %>%
  inner_join(., metadata, by="group") %>%
  filter(combo == "off live" | combo == "in live") 

sig_taxa <- composite %>%
  nest(data = -order) %>%
  mutate(test = map(.x=data, ~pairwise.wilcox.test(x=.x$rel_abund, g=.x$combo, p.adjust.method = "BH") %>% tidy())) %>%
  unnest(test) %>% 
  filter(p.value < 0.05) %>%
  #select(genus, p.value) %>%
  ungroup()

#for order
composite %>%
  inner_join(sig_taxa, by="order") %>%
  mutate(rel_abund = 100 * (rel_abund + 1/20000),
         season = factor(season, levels = c("off", "in"))) %>%
  ggplot(aes(x=rel_abund, y=order, color=season, fill=season)) +
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
  scale_color_manual(name=NULL,
                     breaks=c("in", "off"),
                     values=c('darkgreen', 'dodgerblue3'),
                     labels=c("In-season", "Off-season")) +
  scale_fill_manual(name=NULL,
                    breaks=c("in", "off"),
                    values=c('darkgreen', 'dodgerblue3'),
                    labels=c("In-season", "Off-season")) +
  labs(x= "Relative abundance", y=NULL,
       title="Significantly Different Orders (Fungal)<br>In- vs. Off-season Bottom Farm") +
  theme_light() +
  theme(axis.text.y = element_markdown(size=9),
        legend.text = element_markdown(size=10),
        legend.key.size = unit(15, "pt"),
        legend.position = "bottom",
        plot.title = element_markdown(hjust=0.5))

ggsave("sig.diff.orders.ITS.in.off.BF.png", width=12, height=10)

#for in live vs. off live
composite %>%
  inner_join(sig_taxa, by="order") %>%
  mutate(rel_abund = 100 * (rel_abund + 1/20000),
         combo = factor(combo, levels = c("in live", "off live"))) %>%
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
                     breaks = c("in live", "off live"),
                     values = c("#E97451", "#7393B3"),
                     labels = c("In-season Live", "Off-season Live")) +
  scale_fill_manual(NULL,
                    breaks = c("in live", "off live"),
                    values = c("#E97451", "#7393B3"),
                    labels = c("In-season Live", "Off-season Live"))  +
  labs(x= "Relative abundance", y=NULL,
       title="Significantly Different Orders (Fungal)<br>Bottom Farm In- vs. Off-season") +
  theme_light() +
  theme(axis.text.y = element_markdown(size=9),
        legend.text = element_markdown(size=10),
        legend.key.size = unit(15, "pt"),
        legend.position = "bottom",
        plot.title = element_markdown(hjust=0.5))

ggsave("sig.diff.orders.ITS.ilive.olive.png", width=12, height=10)

##
#looking at sig. diff. taxa between sclerotia treatments
##

metadata <- read_excel("metadata.in.off.24.16S.BF.xlsx") %>%
  rename(group = sample) %>%
  mutate(combo = glue("{season} {treatment}")) %>%
  subset(season !="off") #want to look at each location separately for sclerotia treatment

composite <- inner_join(shared, taxonomy, by="otu") %>%
  group_by(group, order) %>%
  summarize(count = sum(count), .groups="drop") %>%
  group_by(group) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  select(-count) %>%
  inner_join(., metadata, by="group")

sig_taxa <- composite %>%
  nest(data = -order) %>%
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
         season = factor(season, levels = c("in", "off"))) %>%
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

ggsave("sig.diff.orders.16S.in.off.BF.trt.png", width=11, height=8)
