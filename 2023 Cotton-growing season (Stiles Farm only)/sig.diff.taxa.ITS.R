library(tidyverse)
library(broom)
library(ggtext)
library(readxl)

set.seed(19980609)

shared <- read_tsv("final.opti_mcc.0.03.subsample.shared",
                   col_types = cols(Group = col_character(),
                                    .default = col_double())) %>%
  rename_all(tolower) %>%
  select(group, starts_with("otu")) %>%
  pivot_longer(-group, names_to="otu", values_to="count")

taxonomy <- read_tsv(file="final.opti_mcc.0.03.cons.taxonomy") %>%
  rename_all(tolower) %>%
  select (otu, taxonomy) %>%
  mutate(otu = tolower(otu)) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="\\(\\d*\\)", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern=";$", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="k__", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="p__", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="c__", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="o__", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="f__", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="g__", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="s__", replacement="")) %>%
  separate(taxonomy, into=c("kingdom", "phylum", "class", "order", "family", "genus"), sep=";") %>%
  pivot_longer(c("kingdom", "phylum", "class", "order", "family", "genus"),
               names_to="level",
               values_to="taxon") %>%
  subset(level == "class") %>% #genus can be swapped for classification of choice
  select(otu, taxon)

metadata <- read_excel("ITS.2023.metadata.xlsx") %>%
  rename(group = sample)

composite <- inner_join(shared, taxonomy, by="otu") %>%
  mutate(taxon=str_replace(taxon, "unclassified_(.*)", "Unclassified \\1"),
         taxon=str_replace(taxon, "(.*)_unclassified", "Unclassified \\1")) %>%
  group_by(group, taxon) %>%
  summarize(count = sum(count), .groups="drop") %>%
  group_by(group) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  select(-count) %>%
  inner_join(., metadata, by="group") 
  #filter(trt == "live" | trt == "dead")

sig_taxa <- composite %>%
  nest(data = -taxon) %>%
  mutate(test = map(.x=data, ~pairwise.wilcox.test(x=.x$rel_abund, g=.x$trt, p.adjust.method = "BH") %>% tidy())) %>%
  unnest(test) %>% 
  filter(p.value < 0.05) %>%
  ungroup()

composite %>%
  inner_join(sig_taxa, by="taxon") %>%
  mutate(rel_abund = 100 * (rel_abund + 1/20000),
         #need to italicize for family and genus levels
         #taxon = str_replace(taxon, "Unclassified (.*)", "Unclassified *\\1*"),
         #taxon = str_replace(string=taxon,
                              #pattern="^(\\S*)$",
                              #replacement="*\\1*"),
         taxon=str_replace_all(taxon, "_", " "),
         trt = factor(trt, levels = c("live", "dead", "none"))) %>%
  ggplot(aes(x=rel_abund, y=taxon, color=trt, fill=trt)) +
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
                     values = c('#2c7bb6', '#d7191c', "black"),
                     labels = c("Live", "Heat<br>Killed", "Bulk soil")) +
  scale_fill_manual(NULL,
                    breaks = c("live", "dead", "none"),
                    values = c('#2c7bb6', '#d7191c', "black"),
                    labels = c("Live", "Heat<br>Killed", "Bulk soil"))  +
  labs(x= "Relative abundance", y=NULL,
       title="Significantly Different Class (ITS Region)") +
  theme_classic() +
  theme(axis.text.y = element_markdown(),
        legend.text = element_markdown(size=10),
        legend.key.size = unit(15, "pt"),
        plot.title = element_text(hjust=0.5))

ggsave("sig.diff.classes.ITS.png", width=10, height=5)
