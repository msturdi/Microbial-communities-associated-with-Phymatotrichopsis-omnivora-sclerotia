library(tidyverse)
library(readxl)
library(purrr)
library(broom)
library(dbplyr)
library(ggtext)
library(glue)
library(MoMAColors)

display.all.moma(12, colorblind_only=F)

metadata <- read_excel(path="metadata.three.bacterial.xlsx") %>%
  mutate(combo = glue("{season} {treatment}"))

alpha <- read_tsv(file="final.opti_mcc.groups.ave-std.three.bacterial.summary", 
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
t.test(shannon~season, data=meta_alpha,
       var.equal=TRUE,
       conf.level=0.95)

t.test(sobs~season, data=meta_alpha,
       var.equal=TRUE,
       conf.level=0.95)

t.test(shannoneven~season, data=meta_alpha,
       var.equal=TRUE,
       conf.level=0.95)

#assessing sig. of  Shannon
shannon_aov <- aov(shannon~season, data=meta_alpha)
summary(shannon_aov)
TukeyHSD(shannon_aov)

#assessing sig. of raw Shannon values
kruskal.test(shannon~season, data=meta_alpha)
pairwise.wilcox.test(g=meta_alpha$season, x=meta_alpha$shannon, p.adjust.method="BH")

kruskal.test(sobs~season, data=meta_alpha)
pairwise.wilcox.test(g=meta_alpha$season, x=meta_alpha$sobs, p.adjust.method="BH")

kruskal.test(shannoneven~season, data=meta_alpha)
pairwise.wilcox.test(g=meta_alpha$season, x=meta_alpha$shannoneven, p.adjust.method="BH")

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
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(fill=season), position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.3), shape=21, size=1.5) +
  stat_summary(fun.y = mean, geom="point", shape=23, size=2) + 
  scale_color_manual(name=NULL,
                     breaks=c("in23", "off24", "in24"),
                     values=c(moma.colors("Palermo", 3, direction=1, type="discrete")),
                     labels=c("In-season '23", "Off-season '24", "In-season '24")) +
  scale_fill_manual(name=NULL,
                    breaks=c("in23", "off24", "in24"),
                    values=c(moma.colors("Palermo", 3, direction=1, type="discrete")),
                    labels=c("In-season '23", "Off-season '24", "In-season '24")) +
  scale_x_discrete(limits=c("in23", "off24", "in24"),
                   labels=c("In-season '23", "Off-season '24", "In-season '24")) +
  labs(x=NULL,
       y="Shannon",
       title="Comparison of Three Seasons of Bacterial Shannon Diversity<br>Stiles Farm") +
  theme_classic() +
  theme(axis.text.x = element_markdown(size=10),
        axis.text.y = element_markdown(size=10),
        legend.position = "none",
        plot.title = element_markdown(hjust=0.5, size=12, vjust=2))

ggsave("shannon.16S.three.png", width=8, height=6)

#just season sobs
meta_alpha %>%
  ggplot(aes(x=season, y=sobs, color=season)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(fill=season), position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.3), shape=21, size=1.5) +
  stat_summary(fun.y = mean, geom="point", shape=23, size=2) + 
  scale_color_manual(name=NULL,
                     breaks=c("in23", "off24", "in24"),
                     values=c(moma.colors("Palermo", 3, direction=1, type="discrete")),
                     labels=c("In-season '23", "Off-season '24", "In-season '24")) +
  scale_fill_manual(name=NULL,
                    breaks=c("in23", "off24", "in24"),
                    values=c(moma.colors("Palermo", 3, direction=1, type="discrete")),
                    labels=c("In-season '23", "Off-season '24", "In-season '24")) +
  scale_x_discrete(limits=c("in23", "off24", "in24"),
                   labels=c("In-season '23", "Off-season '24", "In-season '24")) +
  labs(x=NULL,
       y="Observed Diversity",
       title="Comparison of Three Seasons of Bacterial Observed Diversity<br>Stiles Farm") +
  theme_classic() +
  theme(axis.text.x = element_markdown(size=10),
        axis.text.y = element_markdown(size=10),
        legend.position = "none",
        plot.title = element_markdown(hjust=0.5, size=12, vjust=2))

ggsave("sobs.16S.three.png", width=8, height=6)

#just season evenness 
meta_alpha %>%
  ggplot(aes(x=season, y=shannoneven, color=season)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(fill=season), position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.3), shape=21, size=1.5) +
  stat_summary(fun.y = mean, geom="point", shape=23, size=2) + 
  scale_color_manual(name=NULL,
                     breaks=c("in23", "off24", "in24"),
                     values=c(moma.colors("Palermo", 3, direction=1, type="discrete")),
                     labels=c("In-season '23", "Off-season '24", "In-season '24")) +
  scale_fill_manual(name=NULL,
                    breaks=c("in23", "off24", "in24"),
                    values=c(moma.colors("Palermo", 3, direction=1, type="discrete")),
                    labels=c("In-season '23", "Off-season '24", "In-season '24")) +
  scale_x_discrete(limits=c("in23", "off24", "in24"),
                   labels=c("In-season '23", "Off-season '24", "In-season '24"))  +
  labs(x=NULL,
       y="Shannon Evenness",
       title="Comparison of Three Seasons of Bacterial Evenness<br>Stiles Farm") +
  theme_classic() +
  theme(axis.text.x = element_markdown(size=10),
        axis.text.y = element_markdown(size=10),
        legend.position = "none",
        plot.title = element_markdown(hjust=0.5, size=12, vjust=2))

ggsave("even.16S.three.png", width=8, height=6)

#sclerotia trts
meta_alpha %>%
  ggplot(aes(x=combo, y=shannon, color=combo)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(fill=combo), position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.3), shape=21, size=1.5) +
  stat_summary(fun.y = mean, geom="point", shape=23, size=2) + 
  scale_color_manual(name=NULL,
                     breaks=c("in23 live", "in23 dead", "in23 none", "off24 live", "off24 dead", "off24 none", "in24 live", "in24 dead", "in24 none"),
                     values=c(moma.colors("OKeeffe", 9, direction=1, type="continuous")),
                     labels=c("In-season '23<br>live", "In-season '23<br>heat-killed", "In-season '23<br>bulk", "Off-season '24<br>live", "Off-season '24<br>heat-killed", 
                              "Off-season '24<br>bulk", "In-season '24<br>live", "In-season '24<br>heat-killed", "In-season '24<br>bulk")) +
  scale_fill_manual(name=NULL,
                    breaks=c("in23 live", "in23 dead", "in23 none", "off24 live", "off24 dead", "off24 none", "in24 live", "in24 dead", "in24 none"),
                    values=c(moma.colors("OKeeffe", 9, direction=1, type="continuous")),
                    labels=c("In-season '23<br>live", "In-season '23<br>heat-killed", "In-season '23<br>bulk", "Off-season '24<br>live", "Off-season '24<br>heat-killed", 
                             "Off-season '24<br>bulk", "In-season '24<br>live", "In-season '24<br>heat-killed", "In-season '24<br>bulk"))  +
  scale_x_discrete(limits=c("in23 live", "in23 dead", "in23 none", "off24 live", "off24 dead", "off24 none", "in24 live", "in24 dead", "in24 none"),
                   labels=c("In-season '23<br>live", "In-season '23<br>heat-killed", "In-season '23<br>bulk", "Off-season '24<br>live", "Off-season '24<br>heat-killed", 
                            "Off-season '24<br>bulk", "In-season '24<br>live", "In-season '24<br>heat-killed", "In-season '24<br>bulk")) +
  labs(x=NULL,
       y="Shannon",
       title="Comparison of In- vs. Off-season Bacterial Diversity<br>Stiles Farm") +
  theme_classic() +
  theme(axis.text.x = element_markdown(size=10),
        axis.text.y = element_markdown(size=10),
        legend.position = "none",
        plot.title = element_markdown(hjust=0.5, size=12, vjust=2))

ggsave("shannon.16S.three.trt.png", width=10, height=6)
