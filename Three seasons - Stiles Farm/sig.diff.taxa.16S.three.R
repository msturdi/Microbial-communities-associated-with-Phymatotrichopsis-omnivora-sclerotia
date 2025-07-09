library(tidyverse)
library(broom)
library(ggtext)
library(readxl)
library(glue)
library(MoMAColors)

set.seed(19980609)

shared <- read_tsv("final.opti_mcc.0.03.subsample.three.bacterial.shared",
                   col_types = cols(Group = col_character(),
                                    .default = col_double())) %>%
  rename_all(tolower) %>%
  select(group, starts_with("otu")) %>%
  pivot_longer(-group, names_to="otu", values_to="count") %>%
  mutate(otu = str_replace(string=otu,
                           pattern="otu",
                           replacement = "Otu"))

taxonomy <- read_tsv(file="final.opti_mcc.0.03.cons.three.bacterial.taxonomy") %>%
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
         order=str_replace(order, "unclassified_(.*)", "Unclassified \\1"),
         order=str_replace(order, "(.*)_unclassified", "Unclassified \\1"),
         order=str_replace_all(order, "_", " "))%>%
  select(otu, order) #can swap in any taxonomic level, make sure to do it above too

metadata <- read_excel("metadata.three.bacterial.xlsx") %>%
  rename(group = sample) %>%
  mutate(combo = glue("{season} {treatment}"))

composite <- inner_join(shared, taxonomy, by="otu") %>%
  group_by(group, order) %>%
  summarize(count = sum(count), .groups="drop") %>%
  group_by(group) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  select(-count) %>%
  inner_join(., metadata, by="group")
  #filter(combo == "off live" | combo == "in live") 

composite_ave <- composite %>%
  subset(rel_abund !="0") %>%
  group_by(order, season) %>%
  summarize(mean=mean(rel_abund),
            se = sd(rel_abund) / sqrt(n()),
            max = mean+se,
            min = mean-se,
            N = n(),
            .groups = "drop")

##
#for three seasons tables comparing rel abund across the seasons
##
composite_ave_scl_trts <- composite %>%
  subset(rel_abund !="0") %>%
  mutate(rel_abund = 100 * (rel_abund + 1/20000)) %>%
  #subset(treatment == "live") %>%
  #subset(treatment == "dead") %>%
  subset(treatment == "none") %>%
  group_by(order, combo) %>%
  summarize(mean=mean(rel_abund),
            se = sd(rel_abund) / sqrt(n()),
            max = mean+se,
            min = mean-se,
            N = n(),
            .groups = "drop")

##
#stat analysis by sclerotia treatment
##
composite_live <- composite %>%
  subset(treatment == "live")

sig_taxa_live <- composite_live %>%
  nest(data = -order) %>%
  mutate(test = map(.x=data, ~pairwise.wilcox.test(x=.x$rel_abund, g=.x$season, p.adjust.method = "BH") %>% tidy())) %>%
  unnest(test) %>% 
  filter(p.value < 0.05) %>%
  #select(genus, p.value) %>%
  ungroup()

composite_dead <- composite %>%
  subset(treatment == "dead")

sig_taxa_dead <- composite_dead %>%
  nest(data = -order) %>%
  mutate(test = map(.x=data, ~pairwise.wilcox.test(x=.x$rel_abund, g=.x$season, p.adjust.method = "BH") %>% tidy())) %>%
  unnest(test) %>% 
  filter(p.value < 0.05) %>%
  #select(genus, p.value) %>%
  ungroup()

composite_bulk <- composite %>%
  subset(treatment == "none")

sig_taxa_bulk <- composite_bulk %>%
  nest(data = -order) %>%
  mutate(test = map(.x=data, ~pairwise.wilcox.test(x=.x$rel_abund, g=.x$season, p.adjust.method = "BH") %>% tidy())) %>%
  unnest(test) %>% 
  filter(p.value < 0.05) %>%
  #select(genus, p.value) %>%
  ungroup()
##

sig_taxa <- composite %>%
  nest(data = -order) %>%
  mutate(test = map(.x=data, ~pairwise.wilcox.test(x=.x$rel_abund, g=.x$season, p.adjust.method = "BH") %>% tidy())) %>%
  unnest(test) %>% 
  filter(p.value < 0.05) %>%
  #select(genus, p.value) %>%
  ungroup()

#for order
composite %>%
  inner_join(sig_taxa, by="order") %>%
  mutate(rel_abund = 100 * (rel_abund + 1/20000),
         season = factor(season, levels = c("in23", "off24", "in24"))) %>%
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
                     breaks=c("in23", "off24", "in24"),
                     values=c(moma.colors("Palermo", 3, direction=1, type="discrete")),
                     labels=c("In-season '23", "Off-season '24", "In-season '24")) +
  scale_fill_manual(name=NULL,
                    breaks=c("in23", "off24", "in24"),
                    values=c(moma.colors("Palermo", 3, direction=1, type="discrete")),
                    labels=c("In-season '23", "Off-season '24", "In-season '24")) +
  labs(x= "Relative abundance", y=NULL,
       title="Significantly Different Orders (Bacterial)<br>Three Seasons at the Stiles Farm") +
  theme_light() +
  theme(axis.text.y = element_markdown(size=9),
        legend.text = element_markdown(size=10),
        legend.key.size = unit(15, "pt"),
        legend.position = "bottom",
        plot.title = element_markdown(hjust=0.5))

ggsave("sig.diff.orders.16S.three.png", width=13, height=16)
