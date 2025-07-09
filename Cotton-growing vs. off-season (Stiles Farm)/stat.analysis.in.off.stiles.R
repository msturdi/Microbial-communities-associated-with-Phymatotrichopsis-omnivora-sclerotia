library(tidyverse)
library(readxl)
library(purrr)
library(broom)
library(dbplyr)
library(ggtext)
library(glue)
library(MoMAColors)

display.all.moma(6, colorblind_only=T)

metadata <- read_excel(path="metadata.16S.in.off.stiles.xlsx") %>%
  mutate(combo = glue("{season} {treatment}"))

alpha <- read_tsv(file="final.opti_mcc.groups.ave-std.bacterial.stiles.summary", 
                  col_types=cols(group = col_character())) %>%
  filter(method=='ave') %>%
  select(group, sobs, invsimpson, shannon, shannoneven, coverage)

meta_alpha <- inner_join(metadata, alpha, by=c('sample'='group'))

meta_alpha %>%
  nest(data = -combo) %>%
  mutate(summary_data=map(data, ~summary(.x$shannon) %>% tidy)) %>%
  #mutate(summary_data=map(data, ~summary(.x$invsimpson) %>% tidy)) %>%
  unnest(cols=summary_data) %>%
  select(-data)

meta_alpha %>% pull(shannon) %>% shapiro.test()
meta_alpha %>% pull(sobs) %>% shapiro.test() 
meta_alpha %>% pull(shannoneven) %>% shapiro.test() 

#tests for homoscedasticity
bartlett.test(shannon~season, data=meta_alpha)
bartlett.test(sobs~season, data=meta_alpha)
bartlett.test(shannoneven~season, data=meta_alpha)

#using a two-sample t-test b/c I have one measurement variable, and one nominal variable that only has 2 groups
t.test(shannon~season, data=meta_alpha, #not sig.
       var.equal=TRUE,
       conf.level=0.95)

t.test(sobs~season, data=meta_alpha, #significant
       var.equal=TRUE,
       conf.level=0.95)

t.test(invsimpson~season, data=meta_alpha, #not sig.
       var.equal=TRUE,
       conf.level=0.95)

t.test(shannoneven~season, data=meta_alpha, #not sig.
       var.equal=TRUE,
       conf.level=0.95)

#assessing sig. of  Shannon
shannon_aov <- aov(shannon~season, data=meta_alpha)
summary(shannon_aov)
TukeyHSD(shannon_aov)

#assessing sig. of raw Shannon values
kruskal.test(shannon~combo, data=meta_alpha)
pairwise.wilcox.test(g=meta_alpha$combo, x=meta_alpha$shannon, p.adjust.method="BH")

shannon_summary_season <- meta_alpha %>%
  group_by(season) %>%
  summarize(mean = mean(shannon),
            se = sd(shannon) / sqrt(n()),
            max = mean+se,
            min = mean-se,
            N = n(),
            .groups = "drop")

#season only
meta_alpha %>%
  ggplot(aes(x=season, y=shannon, color=season)) +
  geom_boxplot(width=0.3, outlier.shape = NA) +
  geom_jitter(aes(fill=season), position = position_jitterdodge(dodge.width = 0.2, jitter.width = 0.2), shape=21, size=1.5) +
  stat_summary(fun.y = mean, geom="point", shape=23, size=2) + 
  scale_color_manual(name=NULL,
                     breaks=c("in", "off"),
                     values=c('darkgreen', 'dodgerblue3'),
                     labels=c("Cotton-growing<br>season", "Off-season")) +
  scale_fill_manual(name=NULL,
                    breaks=c("in", "off"),
                    values=c('darkgreen', 'dodgerblue3'),
                    labels=c("In-season", "Off-season")) +
  scale_x_discrete(limits=c("in", "off"),
                   labels=c("Cotton-growing<br>season", "Off-season")) +
  scale_y_continuous(limits=c(6.35,6.75)) +
  labs(x=NULL,
       y="Bacterial Shannon Diversity",
       title="Cotton-growing vs. Off Season Bacterial Diversity<br>**Stiles Farm (Conducive)**") +
  theme_classic() +
  theme(axis.text.x = element_markdown(size=16),
        axis.text.y = element_markdown(size=12),
        axis.title.y = element_markdown(size=16),
        legend.position = "none",
        plot.title = element_markdown(hjust=0.5, size=18, vjust=2))

ggsave("shannon.16S.in.off.stiles.png", width=9, height=7)

#just season sobs
meta_alpha %>%
  ggplot(aes(x=season, y=sobs, color=season)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(fill=season), position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.3), shape=21, size=1.5) +
  stat_summary(fun.y = mean, geom="point", shape=23, size=2) + 
  scale_color_manual(name=NULL,
                     breaks=c("in", "off"),
                     values=c('darkgreen', 'dodgerblue3'),
                     labels=c("In-season", "Off-season")) +
  scale_fill_manual(name=NULL,
                    breaks=c("in", "off"),
                    values=c('darkgreen', 'dodgerblue3'),
                    labels=c("In-season", "Off-season")) +
  scale_x_discrete(limits=c("in", "off"),
                   labels=c("In-season", "Off-season")) +
  labs(x=NULL,
       y="Observed Diversity",
       title="Comparison of In- vs. Off-season Bacterial Observed Diversity<br>Stiles Farm") +
  theme_classic() +
  theme(axis.text.x = element_markdown(size=10),
        axis.text.y = element_markdown(size=10),
        legend.position = "none",
        plot.title = element_markdown(hjust=0.5, size=12, vjust=2))

ggsave("sobs.16S.in.off.stiles.png", width=8, height=6)

#just season evenness 
meta_alpha %>%
  ggplot(aes(x=season, y=shannoneven, color=season)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(fill=season), position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.3), shape=21, size=1.5) +
  stat_summary(fun.y = mean, geom="point", shape=23, size=2) + 
  scale_color_manual(name=NULL,
                     breaks=c("in", "off"),
                     values=c('darkgreen', 'dodgerblue3'),
                     labels=c("In-season", "Off-season")) +
  scale_fill_manual(name=NULL,
                    breaks=c("in", "off"),
                    values=c('darkgreen', 'dodgerblue3'),
                    labels=c("In-season", "Off-season")) +
  scale_x_discrete(limits=c("in", "off"),
                   labels=c("In-season", "Off-season")) +
  labs(x=NULL,
       y="Shannon Evenness",
       title="Comparison of In- vs. Off-season Bacterial Evenness<br>Stiles Farm") +
  theme_classic() +
  theme(axis.text.x = element_markdown(size=10),
        axis.text.y = element_markdown(size=10),
        legend.position = "none",
        plot.title = element_markdown(hjust=0.5, size=12, vjust=2))

ggsave("even.16S.in.off.stiles.png", width=8, height=6)

#sclerotia trts
meta_alpha %>%
  ggplot(aes(x=combo, y=shannon, color=combo)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(fill=combo), position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.3), shape=21, size=1.5) +
  stat_summary(fun.y = mean, geom="point", shape=23, size=2) + 
  scale_color_manual(name=NULL,
                     breaks=c("in live", "in dead", "in none", "off live", "off dead", "off none"),
                     values=c(moma.colors("OKeeffe", 6, direction=1, type="discrete")),
                     labels=c("In-season<br>live", "In-season<br>heat-killed", "In-season<br>bulk", "Off-season<br>live", "Off-season<br>heat-killed", "Off-season<br>bulk")) +
  scale_fill_manual(name=NULL,
                    breaks=c("in live", "in dead", "in none", "off live", "off dead", "off none"),
                    values=c(moma.colors("OKeeffe", 6, direction=1, type="discrete")),
                    labels=c("In-season<br>live", "In-season<br>heat-killed", "In-season<br>bulk", "Off-season<br>live", "Off-season<br>heat-killed", "Off-season<br>bulk"))  +
  scale_x_discrete(limits=c("in live", "in dead", "in none", "off live", "off dead", "off none"),
                   labels=c("In-season<br>live", "In-season<br>heat-killed", "In-season<br>bulk", "Off-season<br>live", "Off-season<br>heat-killed", "Off-season<br>bulk")) +
  labs(x=NULL,
       y="Shannon",
       title="Comparison of In- vs. Off-season Bacterial Diversity<br>Stiles Farm") +
  theme_classic() +
  theme(axis.text.x = element_markdown(size=10),
        axis.text.y = element_markdown(size=10),
        legend.position = "none",
        plot.title = element_markdown(hjust=0.5, size=12, vjust=2))

ggsave("shannon.16S.in.off.stiles.trt.png", width=8, height=6)

##
#location faceted sobs, invsimpson, evenness
##

#combining the other metics into one column
meta_alpha_long <- meta_alpha %>%
  pivot_longer(cols = c(sobs, invsimpson, shannoneven),  
               names_to = "other.metric",                 
               values_to = "value")

#pulling maximum value for significance stars
meta_alpha_long %>%
  group_by(other.metric) %>%
  summarize(max_value = max(value, na.rm = TRUE))

meta_alpha_long$other.metric <- factor(meta_alpha_long$other.metric, levels = c("sobs", "invsimpson", "shannoneven"))
location_facet <- meta_alpha_long %>%
  ggplot(aes(x=season, y=value, color=season)) +
  geom_boxplot(width = 0.3, outlier.shape = NA) +
  geom_jitter(aes(fill=season), position = position_jitterdodge(dodge.width = 0.2, jitter.width = 0.2), shape=21, size=1.5) +
  stat_summary(fun.y = mean, geom="point", shape=23, size=2) + 
  facet_wrap(~ other.metric,
             scales = "free_y",
             nrow = 1,
             labeller = labeller(other.metric = c("sobs" = "Observed",
                                                  "invsimpson" = "Inverse Simpson",
                                                  "shannoneven" = "Evenness"))) +
  scale_color_manual(name=NULL,
                     breaks=c("in", "off"),
                     values=c('darkgreen', 'dodgerblue3'),
                     labels=c("Cotton-growing<br>season", "Off-season")) +
  scale_fill_manual(name=NULL,
                    breaks=c("in", "off"),
                    values=c('darkgreen', 'dodgerblue3'),
                    labels=c("Cotton-growing<br>season", "Off-season")) +
  scale_x_discrete(limits=c("in", "off"),
                   labels=c("Cotton-growing<br>season", "Off-season")) +
  labs(x=NULL,
       y=NULL,
       title="In- vs. Off-season Additional Bacterial Diversity Indices<br>Stiles Farm") +
  theme_light() +
  theme(axis.text.x = element_markdown(size=16),
        axis.text.y = element_markdown(size=12),
        strip.text = element_markdown(size=16),
        legend.position = "none",
        plot.title = element_markdown(hjust=0.5, size=18, vjust=2))

text_sobs <- tibble(other.metric = factor("sobs", levels = c("sobs", "invsimpson", "shannoneven")),
                    x = 1,
                    y = 4120,
                    label = "*")

location_facet +
  geom_text(data=text_sobs, 
            aes(x=x, y=y, label=label), size=9, 
            inherit.aes = FALSE) 

ggsave("other.metrics.16S.in.off.stiles.location.png", width=14, height=8)
