library(tidyverse)
library(readxl)
library(purrr)
library(broom)
library(dbplyr)
library(ggtext)
library(glue)
library(MoMAColors)

metadata <- read_excel(path="mar.24.fungal.metadata.xlsx") %>%
  mutate(combo = glue("{location} {treatment}"))

alpha <- read_tsv(file="final.opti_mcc.groups.ave-std.mar.24.fungal.summary", 
                  col_types=cols(group = col_character())) %>%
  filter(method=='ave') %>%
  select(group, sobs, invsimpson, shannon, shannoneven, coverage)

meta_alpha <- inner_join(metadata, alpha, by=c('sample'='group'))

meta_alpha %>%
  nest(data = -location) %>%
  mutate(summary_data=map(data, ~summary(.x$shannon) %>% tidy)) %>%
  #mutate(summary_data=map(data, ~summary(.x$invsimpson) %>% tidy)) %>%
  unnest(cols=summary_data) %>%
  select(-data)

#a way to test if the data are normally distributed, if it is, the points should all fall on a line
ggplot(meta_alpha, aes(
  sample=shannon,
  #sample=invsimpson,
  group=location, color=location)) + geom_qq() + stat_qq_line()

#scales the data if it is not normally distributed 
meta_alpha <- mutate(meta_alpha, 
                     scaled_shannon=shannon^3
                     #scaled_invsimpson=invsimpson^3
)

meta_alpha %>% pull(shannon) %>% shapiro.test() 
meta_alpha %>% pull(scaled_shannon) %>% shapiro.test() #normal

#tests for homoscedasticity
bartlett.test(shannon~location, data=meta_alpha)

#using a two-sample t-test b/c I have one measurement variable, and one nominal variable that only has 2 groups
t.test(shannon~location, data=meta_alpha, #not sig.
       var.equal=TRUE,
       conf.level=0.95)

t.test(sobs~location, data=meta_alpha, #not sig.
       var.equal=TRUE,
       conf.level=0.95)

t.test(invsimpson~location, data=meta_alpha, #significant
       var.equal=TRUE,
       conf.level=0.95)

t.test(shannoneven~location, data=meta_alpha, #not sig.
       var.equal=TRUE,
       conf.level=0.95)

kruskal.test(shannon~location, data=meta_alpha)
pairwise.wilcox.test(g=meta_alpha$combo, x=meta_alpha$scaled_shannon, p.adjust.method="BH")

shannon_summary_location <- meta_alpha %>%
  group_by(location) %>%
  summarize(mean = mean(scaled_shannon),
            se = sd(scaled_shannon) / sqrt(n()),
            max = mean+se,
            min = mean-se,
            N = n(),
            .groups = "drop")

#a way to make a bar plot with error bars without needing to have preceeding summary steps
#stat_summary(fun.data = mean_se, geom="errorbar", width=0.5) +
#stat_summary(fun.data = mean_se, geom="bar", show.legend = FALSE)

#location only
meta_alpha %>%
  ggplot(aes(x=location, y=shannon, color=location)) +
  geom_boxplot(width = 0.3, outlier.shape = NA) +
  geom_jitter(aes(fill=location), position = position_jitterdodge(dodge.width = 0.2, jitter.width = 0.2), shape=21, size=1.5) +
  stat_summary(fun.y = mean, geom="point", shape=23, size=2) + 
  scale_color_manual(name=NULL,
                     breaks=c("BF", "Stiles"),
                     values=c('orange2', 'midnightblue'),
                     labels=c("Bottom Farm", "Stiles Farm")) +
  scale_fill_manual(name=NULL,
                    breaks=c("BF", "Stiles"),
                    values=c('orange2', 'midnightblue'),
                    labels=c("Bottom Farm", "Stiles Farm")) +
  scale_x_discrete(limits=c("BF", "Stiles"),
                   labels=c("Bottom Farm<br>(Non-conducive)", "Stiles Farm<br>(Conducive)")) +
  labs(x=NULL,
       y="Fungal Shannon Diversity",
       title="Off-season Fungal Shannon Diversity") +
  theme_classic() +
  theme(axis.text.x = element_markdown(size=16),
        axis.text.y = element_markdown(size=12),
        axis.title.y = element_markdown(size=16),
        legend.position = "none",
        plot.title = element_text(hjust=0.5, size=18, vjust=2))

ggsave("shannon.ITS.mar.24.location.png", width=9, height=7)

#for sclerotia trts
meta_alpha %>%
  ggplot(aes(x=combo, y=sobs, color=combo)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(fill=combo), position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.3), shape=21, size=1.5) +
  stat_summary(fun.y = mean, geom="point", shape=23, size=2) +
  scale_color_manual(name=NULL,
                     breaks=c("BF dead", "BF live", "BF none", "Stiles dead", "Stiles live", "Stiles none"),
                     values=c(moma.colors("Klein", 6, direction=1, type="discrete")),
                     labels=c("BF dead", "BF live", "BF bulk", "Stiles dead", "Stiles live", "Stiles bulk")) +
  scale_fill_manual(name=NULL,
                    breaks=c("BF dead", "BF live", "BF none", "Stiles dead", "Stiles live", "Stiles none"),
                    values=c(moma.colors("Klein", 6, direction=1, type="discrete")),
                    labels=c("BF dead", "BF live", "BF bulk", "Stiles dead", "Stiles live", "Stiles bulk")) +
  scale_x_discrete(limits=c("BF dead", "BF live", "BF none", "Stiles dead", "Stiles live", "Stiles none"),
                   labels=c("BF<br>heat-killed", "BF<br>live", "BF<br>bulk", "Stiles<br>heat-killed", "Stiles<br>live", "Stiles<br>bulk")) +
  labs(x=NULL,
       y="Shannon",
       title="Comparison of Off-season Fungal Shannon Diversity") +
  theme_classic() +
  theme(axis.text.x = element_markdown(size=10),
        axis.text.y = element_markdown(size=10),
        legend.position = "none",
        plot.title = element_text(hjust=0.5, size=12, vjust=2))

ggsave("shannon.ITS.mar.24.trt.png", width=8, height=6)

##
#location faceted sobs, invsimpson, evenness
##

#combining the other metrics into one column
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
  ggplot(aes(x=location, y=value, color=location)) +
  geom_boxplot(width = 0.3, outlier.shape = NA) +
  geom_jitter(aes(fill=location), position = position_jitterdodge(dodge.width = 0.2, jitter.width = 0.2), shape=21, size=1.5) +
  stat_summary(fun.y = mean, geom="point", shape=23, size=2) + 
  facet_wrap(~ other.metric,
             scales = "free_y",
             nrow = 1,
             labeller = labeller(other.metric = c("sobs" = "Observed",
                                                  "invsimpson" = "Inverse Simpson",
                                                  "shannoneven" = "Evenness"))) +
  scale_color_manual(name=NULL,
                     breaks=c("BF", "Stiles"),
                     values=c('orange2', 'midnightblue'),
                     labels=c("Bottom Farm", "Stiles Farm")) +
  scale_fill_manual(name=NULL,
                    breaks=c("BF", "Stiles"),
                    values=c('orange2', 'midnightblue'),
                    labels=c("Bottom Farm", "Stiles Farm")) +
  scale_x_discrete(limits=c("BF", "Stiles"),
                   labels=c("Bottom Farm<br>(Non-conducive)", "Stiles Farm<br>(Conducive)")) +
  labs(x=NULL,
       y=NULL,
       title="Off-season Additional Fungal Diversity Indices") +
  theme_light() +
  theme(axis.text.x = element_markdown(size=16),
        axis.text.y = element_markdown(size=12),
        strip.text = element_markdown(size=16),
        legend.position = "none",
        plot.title = element_text(hjust=0.5, size=18, vjust=2))

text_invsimpson <- tibble(other.metric = factor("invsimpson", levels = c("sobs", "invsimpson", "shannoneven")),
                          x = 2,
                          y = 20.5,
                          label = "*")

location_facet +
  geom_text(data=text_invsimpson, 
            aes(x=x, y=y, label=label), size=9, 
            inherit.aes = FALSE)

ggsave("other.metrics.fungal.mar.24.location.png", width=14, height=8)
