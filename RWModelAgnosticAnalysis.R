#Probabilistic Analysis

library(easystats)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggrepel)

source("Functions.R")
theme_set(cowplot::theme_cowplot())

#load data

probabilistic   <- read.csv("Data/ProbabilisticCleaned.csv") %>% dplyr::select(-X)

names(probabilistic)

# Plotting ----------------------------------------------------------------

probabilisticPlot <- probabilistic %>%
  mutate(
    PersecOrd = cut(
      Persec,
      breaks = c(-Inf, 4, 10, 17, 27, Inf),
      labels = c("1", "2", "3", "4", "5")
    ),
    SocRefOrd = cut(
      SocRef,
      breaks = c(-Inf, 9, 15, 20, 24, Inf),
      labels = c("1", "2", "3", "4", "5")
    ),
    PersecOrd = factor(
      PersecOrd,
      levels = c("1", "2", "3", "4", "5"),
      ordered = T
    ),
    PersecLevel = ifelse(Persec >= 3.88, "High", "Low"),
    SocRefLevel = ifelse(SocRef >= 7.66, "High", "Low"),
    ICARLevel   = ifelse(ICARTot >=5.00, "High", "Low"),
    SocRefOrd = factor(
      SocRefOrd,
      levels = c("1", "2", "3", "4", "5"),
      ordered = T
    ),
    ICAROrd = cut(
      ICARTot,
      breaks = (5),
      labels = c("1", "2", "3", "4", "5")
    ),
    ICAROrd = factor(
      ICAROrd,
      levels = c("1", "2", "3", "4", "5"),
      ordered = T
    )
  ) %>%
  pivot_longer(Urn1:Urn3, names_to = "Urns", values_to = "Value")



#proportion correct by paranoia
prop1 <- ggplot(probabilisticPlot) +

  stat_summary(aes(Trial, Correct, color = PersecLevel), geom = "line") +

  scale_color_manual(name = "Persecutory Ideation", values = c("#A50F15", "#FC9272")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0,  0.5,  1))+
  coord_cartesian(ylim = c(0,1))+

  labs(x = "Trial",
       y = "Proportion of highest earning card chosen") +

  facet_wrap( ~ Block) +

  bbplot::bbc_style() +

  theme(
    strip.text.x = element_text(hjust = 0.5),
    axis.title = element_text(size = 14),
    legend.position = c(0.1, 0.9),
    legend.title = element_text(),
    legend.direction = "horizontal"
  )

library(RColorBrewer)

#proportion correct by ICAR
prop2 <- ggplot(probabilisticPlot) +

  stat_summary(aes(Trial, Correct, color = ICARLevel), geom = "line") +

  scale_color_manual(name = "ICAR Total", values = c("#006D2C", "#A1D99B")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1))+
  coord_cartesian(ylim = c(0,1))+

  labs(x = "Trial",
       y = "Proportion of Correct Answers") +

  facet_wrap( ~ Block) +

  bbplot::bbc_style() +

  theme(
    strip.text.x = element_text(hjust = 0.5),
    axis.title = element_text(size = 14),
    legend.position = c(0.2, 0.9),
    legend.title = element_text(),
    legend.direction = "horizontal"
  )

sum1 <- probabilisticPlot %>%
  dplyr::select(Persec, Value, Urns, ID, Block) %>%
  mutate(Urns = recode(Urns,
                       Urn1 = 'Card1',
                       Urn2 = 'Card2',
                       Urn3 = 'Card3')) %>%
  distinct() %>%
  ggplot() +

  geom_smooth(aes(Persec,
                  Value,
                  color = Urns),
              method = "lm") +
  ggpubr::stat_cor(
    aes(
      Persec,
      Value,
      color = Urns),
      show.legend = F,
      label.y.npc = 0.8
    ) +

  scale_color_brewer(palette = 'Dark2') +
  scale_y_continuous(labels = c(0, 5, 10, 15, 20, 25), breaks = c(0, 5, 10, 15, 20, 25))+
  coord_cartesian(ylim = c(0, 25))+

  labs(y = "Sum of each card chosen in each block",
       x = "Persecutory Ideation") +

  facet_wrap(~ Block) +
  bbplot::bbc_style() +

  theme(
    strip.text.x = element_blank(),
    axis.title = element_text(size = 14),
    legend.position = c(0.7, 1),
    legend.title = element_text(),
    legend.direction = "vertical"
  )

sum2 <- probabilisticPlot %>%
    dplyr::select(ICARTot, Value, Urns, ID, Block) %>%
    mutate(ICARTot = as.numeric(ICARTot)) %>%
    distinct() %>%

  ggplot() +

  geom_smooth(
    aes(
    ICARTot, Value, color = Urns),
    method = "lm") +

  scale_color_brewer(palette = 'Dark2') +
  scale_y_continuous(labels = c(0, 5, 10, 15, 20, 25), breaks = c(0, 5, 10, 15, 20, 25))+
  coord_cartesian(ylim = c(0, 25))+

  ggpubr::stat_cor(aes(ICARTot, Value, color = Urns), label.y.npc = 0.8) +

  labs(
    y = "Sum of each card chosen in each block",
    x = "General Reasoning Ability (ICAR Total Score)"
    ) +

  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10), labels = c(0, 2, 4, 6, 8, 10))+

  facet_wrap( ~ Block) +
  bbplot::bbc_style() +

  theme(
    strip.text.x = element_blank(),
    axis.title = element_text(size = 14),
    legend.position = 'none',
    legend.title = element_text(),
    legend.direction = "vertical"
  )

library(patchwork)

#Create text annotation for plots
ann_text1 <- data.frame(Persec = 15, Value = 2, lab = "Card 1 highest earner",
                        Block = factor("Block 1", levels = c("Block 1", "Block 2")))
ann_text2 <- data.frame(Persec = 15, Value = 2, lab = "Card 3 highest earner",
                       Block = factor("Block 2", levels = c("Block 1", "Block 2")))
ann_text3 <- data.frame(ICARTot = 5, Value = 2, lab = "Card 1 highest earner",
                        Block = factor("Block 1", levels = c("Block 1", "Block 2")))
ann_text4 <- data.frame(ICARTot = 5, Value = 2, lab = "Card 3 highest earner",
                        Block = factor("Block 2", levels = c("Block 1", "Block 2")))

#Build plot window with plot objects
(prop1 +
    theme(
     legend.text = element_text(size = 14)) |

 prop2 +
    theme(
     axis.title.y = element_blank(),
     legend.text = element_text(size = 14),
     axis.text.y = element_blank()))/

(sum1 +
    theme(
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
      legend.direction = 'horizontal') +

   geom_text(
     data = ann_text1,
        aes(Persec, Value),
     label = "Card 1 highest earner", family = "Arial", size = 3) +
   geom_text(
     data = ann_text2,
        aes(Persec, Value),
     label = "Card 3 highest earner", family = "Arial", size = 3)|

 sum2 +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      legend.title = element_blank()) +

  geom_text(
    data = ann_text3,
        aes(ICARTot, Value),
    label = "Card 1 highest earner", family = "Arial", size = 3) +
  geom_text(
    data = ann_text4,
        aes(ICARTot, Value),
    label = "Card 3 highest earner", family = "Arial", size = 3)
 )

# Proportion of cards selected by trial

probabilisticPlot %>%
  mutate(Selection = recode(Selection,
                       urn1.png = 'Card1',
                       urn2.png = 'Card2',
                       urn3.png = 'Card3')) %>%
  distinct() %>%
  ggplot() +

  geom_bar(aes(Trial, fill = Selection),
           stat = 'count',
           position = 'fill') +

  scale_fill_brewer(name = "Selection", palette = "Dark2")+

  facet_wrap( ~ Block)

# Analysis ----------------------------------------------------------------



# learning by correct answers
# block 1

#unadjusted
summary(lme4::glmer(Correct ~ scale(Persec) + (1|ID), data = probabilistic %>% filter(Block == 'Block 1'), family = 'binomial'))
confint(lme4::glmer(Correct ~ scale(Persec) + (1|ID), data = probabilistic %>% filter(Block == 'Block 1'), family = 'binomial'))

#adjusted # Model P1
modelp1 <- probabilistic %>% filter(Block == 'Block 1') %>% mutate(Sex = ifelse(Sex == "Male", 1, 0))
globalModel1 <-
  lme4::glmer(
    Correct ~ scale(Persec) + scale(ICARTot) + Age + Sex + Control + (1|ID),
    data = modelp1,
    na.action = na.fail,
    family = 'binomial'
  )

model.compare(globalModel1)

#Model P1b
urnsassess1b <- probabilisticPlot %>%
  dplyr::select(ID, Urns, Value, Block, Persec, ICARTot, Age, Sex, Control) %>%
  filter(Urns == 'Urn2',
         Block == 'Block 1') %>%
  distinct()

model.compare(lm(Value ~ Persec + scale(ICARTot) + scale(Age) + Sex + Control,
                 data = urnsassess1b,
                 na.action = na.fail))

#Model P1c
urnsassess1c <- probabilisticPlot %>%
  dplyr::select(ID, Urns, Value, Block, Persec, ICARTot, Age, Sex, Control) %>%
  filter(Urns == 'Urn3',
         Block == 'Block 1') %>%
  distinct()

summary(lm(Value ~ Persec + scale(ICARTot) + scale(Age) + Sex + Control,
                 data = urnsassess1c,
                 na.action = na.fail))

# block 2

#unadjusted
summary(lme4::glmer(Correct ~ scale(Persec) + (1|ID), data = probabilistic %>% filter(Block == 'Block 2'), family = 'binomial'))
confint(lme4::glmer(Correct ~ scale(Persec) + (1|ID), data = probabilistic %>% filter(Block == 'Block 2'), family = 'binomial'))

#adjusted
#Model P2
modelp2<- probabilistic %>% filter(Block == 'Block 2') %>% mutate(Sex = ifelse(Sex == "Male", 1, 0))
globalModel2 <-
  lme4::glmer(
    Correct ~ scale(Persec) + scale(ICARTot) + Age + Sex + Control + (1|ID),
    data = modelp2,
    na.action = na.fail,
    family = 'binomial'
  )

model.compare(globalModel2)

#Model P2b
urnsassess2b <- probabilisticPlot %>%
  dplyr::select(ID, Urns, Value, Block, Persec, ICARTot, Age, Sex, Control) %>%
  filter(Urns == 'Urn2',
         Block == 'Block 2') %>%
  distinct()

model.compare(lm(Value ~ Persec + scale(ICARTot) + scale(Age) + Sex + Control,
                 data = urnsassess2b,
                 na.action = na.fail))

#Model P2c
urnsassess2c <- probabilisticPlot %>%
  dplyr::select(ID, Urns, Value, Block, Persec, ICARTot, Age, Sex, Control) %>%
  filter(Urns == 'Urn1',
         Block == 'Block 2') %>%
  distinct()

model.compare(lm(Value ~ Persec + scale(ICARTot) + scale(Age) + Sex + Control,
                 data = urnsassess2c,
                 na.action = na.fail))

#learning by explicit answers block 1
modelp3a <- probabilistic %>%
  filter(Block == "Block 1") %>%
  mutate(Sex = ifelse(Sex == "Male", 1, 0)) %>%
  dplyr::select(Final1, Persec, ICARTot, Age, Sex, ID, Control) %>%
  distinct()

#unadjusted
summary(glm(Final1 ~ Persec, modelp3a, family = 'binomial'))
confint(glm(Final1 ~ Persec, modelp3a, family = 'binomial'))

#adjusted

globalModel3a <-
  glm(
    Final1 ~ scale(Persec) + scale(ICARTot) + Age + Sex + Control,
    data = modelp3,
    na.action = na.fail,
    family = 'binomial'
  )

model.compare(globalModel3a)

#learning by explicit answers block 2
modelp3b <- probabilistic %>%
  filter(Block == "Block 2") %>%
  mutate(Sex = ifelse(Sex == "Male", 1, 0)) %>%
  dplyr::select(Final1, Persec, ICARTot, Age, Sex, ID, Control) %>%
  distinct()

#unadjusted
summary(glm(Final1 ~ Persec, modelp3b, family = 'binomial'))
confint(glm(Final1 ~ Persec, modelp3b, family = 'binomial'))

#adjusted

globalModel4 <-
  glm(
    factor(Final1) ~ scale(Persec) + scale(ICARTot) + Age + Sex + Control,
    data = modelp3b,
    na.action = na.fail,
    family = 'binomial'
  )

model.compare(globalModel3b)


#Model P4
TotalRewarddf <- probabilistic %>%
  group_by(ID, Block) %>%
  mutate(TotalReward = sum(Response)) %>%
  dplyr::select(Persec, ICARTot, ID, TotalReward, Age, Sex, Control, Block) %>%
  na.omit() %>%
  distinct()

TotalRewarddf_block1 <- TotalRewarddf %>% filter(Block == 'Block 1')
TotalRewarddf_block2 <- TotalRewarddf %>% filter(Block == 'Block 2')

confint(lm(scale(TotalReward) ~ scale(Persec),
              data = TotalRewarddf_block1,
              na.action = na.fail))
confint(lm(scale(TotalReward) ~ scale(Persec),
              data = TotalRewarddf_block2,
              na.action = na.fail))

model.compare(lm(scale(TotalReward) ~ scale(Persec) + scale(ICARTot) + scale(Age) + Sex + scale(Control),
              data = TotalRewarddf_block1,
              na.action = na.fail))
model.compare(lm(scale(TotalReward) ~ scale(Persec) + scale(ICARTot) + scale(Age) + Sex + scale(Control),
              data = TotalRewarddf_block2,
              na.action = na.fail))

ggplot(probabilistic %>%
         mutate(Selection = ifelse(Selection == 'urn1.png', 1, ifelse(Selection == 'urn2.png', 2, 3)),
                `Card 1` = ifelse(Selection == 1, 1, 0),
                `Card 2` = ifelse(Selection == 2, 1, 0),
                `Card 3` = ifelse(Selection == 3, 1, 0)) %>%
         pivot_longer(`Card 1`:`Card 3`, 'Urn_p', values_to = 'Value'),
       aes(Trial, Value, color = Persec>1))+
  stat_summary(geom = 'line')+
  ggpubr::stat_compare_means(label = 'p.signif', hide.ns = T, label.y = 0.4, size = 6, show.legend = F)+
  facet_wrap(Urn_p~Block, nrow = 3, ncol = 2)+
  scale_color_brewer(palette = 'Set1', name = 'Paranoia Score Over The Median')+
  labs(y = 'p(Card Chosen)')+
  theme_minimal()+
  theme(text = element_text(size = 14),
        legend.position = c(0.25, 0.22),
        legend.background = element_rect(colour = 'black'))
