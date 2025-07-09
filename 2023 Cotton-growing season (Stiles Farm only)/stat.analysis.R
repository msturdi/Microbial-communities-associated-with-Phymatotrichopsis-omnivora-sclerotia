library(tidyverse)
library(readxl)
library(purrr)
library(broom)
library(dbplyr)
library(ggtext)

metadata <- read_excel(path="ITS.2023.metadata.xlsx")

alpha <- read_tsv(file="final.opti_mcc.groups.ave-std.summary", 
                  col_types=cols(group = col_character())) %>%
  filter(method=='ave') %>%
  select(group, sobs, invsimpson, shannon, coverage)

meta_alpha <- inner_join(metadata, alpha, by=c('sample'='group'))

meta_alpha %>%
  nest(data = -trt) %>%
  mutate(summary_data=map(data, ~summary(.x$shannon) %>% tidy)) %>%
  #mutate(summary_data=map(data, ~summary(.x$invsimpson) %>% tidy)) %>%
  unnest(cols=summary_data) %>%
  select(-data)

#a way to test if the data are normally distributed, if it is, the points should all fall on a line
ggplot(meta_alpha, aes(
                      sample=shannon,
                      #sample=invsimpson,
                      group=trt, color=trt)) + geom_qq() + stat_qq_line()

#scales the data if it is not normally distributed 
meta_alpha <- mutate(meta_alpha, 
                     scaled_shannon=shannon^3
                     #scaled_invsimpson=invsimpson^3
                     )

meta_alpha %>% pull(shannon) %>% shapiro.test() #just barely normal
meta_alpha %>% pull(scaled_shannon) %>% shapiro.test() #comfortable normal

meta_alpha %>% pull(invsimpson) %>% shapiro.test() 

#assessing sig. of scaled Shannon
trt_shannon_aov <- aov(scaled_shannon~trt, data=meta_alpha)
summary(trt_shannon_aov)
TukeyHSD(trt_shannon_aov) #run if experimental-wide p-value<0.05

#assessing sig. of raw Shannon values
kruskal.test(shannon~trt, data=meta_alpha)
pairwise.wilcox.test(g=meta_alpha$trt, x=meta_alpha$shannon, p.adjust.method="BH")

#assessing sig. of invsimpson values
trt_invsimpson_aov <- aov(invsimpson~trt, data=meta_alpha)
summary(trt_invsimpson_aov)
TukeyHSD(trt_invsimpson_aov) #run if experimental-wide p-value<0.05

shannon_summary <- meta_alpha %>%
  group_by(trt) %>%
  summarize(mean = mean(shannon),
            se = sd(shannon) / sqrt(n()),
            max = mean+se,
            min = mean-se,
            N = n(),
            .groups = "drop")

invsimpson_summary <- meta_alpha %>%
  group_by(trt) %>%
  summarize(mean = mean(invsimpson),
            se = sd(invsimpson) / sqrt(n()),
            max = mean+se,
            min = mean-se,
            N = n(),
            .groups = "drop")

#a way to make a bar plot with error bars without needing to have preceeding summary steps
#stat_summary(fun.data = mean_se, geom="errorbar", width=0.5) +
#stat_summary(fun.data = mean_se, geom="bar", show.legend = FALSE)

meta_alpha %>%
  ggplot(aes(x=trt, y=shannon, color=trt, fill=trt)) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.5), shape=21, size=1.8) +
  stat_summary(fun.data = median_hilow, fun.args = list(conf.int=0.5),
               geom="pointrange",
               position = position_dodge(width=0.1),
               color="darkgrey", show.legend = FALSE,
               size = 0.1) +  
  scale_color_manual(name=NULL,
                     breaks=c("live", "dead", "none"),
                     values=c('#2c7bb6', '#d7191c', 'black'),
                     labels=c("Live", "Heat killed", "Bulk soil")) +
  scale_fill_manual(name=NULL,
                    breaks=c("live", "dead", "none"),
                    values=c('#2c7bb6', '#d7191c', 'black'),
                    labels=c("Live", "Heat<br>killed", "Bulk soil")) +
  scale_x_discrete(limits=c("live", "dead", "none"),
                   labels=c("Live", "Heat killed", "Bulk soil")) +
  labs(x=NULL,
       y="Scaled Shannon",
       title="Comparison of Shannon Values between Treatments") +
  theme_classic() +
  theme(axis.text.x = element_markdown(),
        title = element_markdown(),
        legend.position = "none",
        plot.title = element_text(hjust=0.5, size=12, vjust=2))

ggsave("shannon.ITS.png", width=8, height=6)
