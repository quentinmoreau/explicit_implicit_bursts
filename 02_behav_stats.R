library(lme4)
library(car)
library(emmeans)
library(dplyr)
library(doBy)

# Read data
data_path <- "data/derivatives/behav_df_cleaned_new.csv"
data <- read.csv(data_path)

data <- data[data$reach_vis_abs_err<60,]

# Factorize variables
data$coh_cat <- ordered(data$coh_cat, levels = c('zero', 'low', 'med', 'high'))
data$group <- as.factor(data$group)

# Create block_type variable
data$block_type <- ifelse(data$block == 0, 'baseline',
                          ifelse(data$block > 0 & data$block < 7, 'adaptation',
                                 ifelse(data$block == 7, 'washout', NA)))
data$block_type <- as.factor(data$block_type)

########################## Absolute error magnitude #############################
print('Absolute error magnitude')
model <- lmer(reach_vis_abs_err ~ group*block_type + (1| subject), data = data)
print(Anova(model, contrasts=list(group=contr.sum, block_type=contr.sum, type=3)))
print(emmeans(model, pairwise~block_type|group, type='response', pbkrtest.limit = 15987))
print(emmeans(model, pairwise~group|block_type, type='response', pbkrtest.limit = 15987))

########################## Signed error - Implicit Group ##########################
print('Signed error - Implicit Group')
data_implicit <- data[data$group == "Implicit", ]

model <- lmer(reach_vis_err ~ block_type + (1| subject), data = data_implicit)
print(Anova(model, contrasts=list(block_type=contr.sum, type=3)))
print(emmeans(model, pairwise~block_type, type='response', pbkrtest.limit = 15829))


########################## Signed error - Implicit Group within Adaptation ##########################
print('Signed error - Implicit Group within Adaptation')
adapt_blk_data_implicit=data_implicit[data_implicit$block>0 & data_implicit$block<7,]

model <- lmer(reach_vis_err ~ block + (1| subject), data = adapt_blk_data_implicit)
print(Anova(model))


########################## Analysis for Coherence and Error magnitude ##########################
print('Analysis for Coherence and Error magnitude - Implicit group')
adapt_blk_data_implicit=data_implicit[data_implicit$block>0 & data_implicit$block<7,]
model <- lmer(reach_vis_abs_err ~ coh_cat * block + (1| subject), data = adapt_blk_data_implicit)
print(Anova(model))

print('Analysis for Coherence and Error magnitude - Explicit group')
data_explicit <- data[data$group == "Explicit", ]
adapt_blk_data_explicit=data_explicit[data_explicit$block>0 & data_explicit$block<7,]
model <- lmer(reach_vis_abs_err ~ coh_cat * block + (1| subject), data = adapt_blk_data_explicit)
print(Anova(model))
em <- emtrends(model, pairwise ~ coh_cat, var = "block", pbkrtest.limit = 12010)
print(summary(em, infer = c(TRUE, TRUE)))






############################ Analysis Comparing Block Types ############################
print('RT - Analysis Comparing Block Types')
model <- lmer(reach_rt ~ block_type * group + (1| subject), data = data)
print(Anova(model, contrasts=list(group=contr.sum, block_type=contr.sum), type=3))
print(emmeans(model, pairwise ~ block_type|group, type='response', pbkrtest.limit = 15829))
print(emmeans(model, pairwise ~ group|block_type, type='response', pbkrtest.limit = 15829))

############################ Analysis Comparing Blocks Within Adaptation  ############################
print('RT - Analysis Comparing Blocks Within Adaptation - Implicit')
# subset the data only for adaptation blocks
data_implicit <- data[data$group == "Implicit", ]
adapt_blk_data_implicit=data_implicit[data_implicit$block>0 & data_implicit$block<7,]

model <- lmer(reach_rt ~ coh_cat * block + (1| subject), data = adapt_blk_data_implicit)
print(Anova(model, contrasts=list(coh_cat=contr.sum), type=3))

print('RT - Analysis Comparing Blocks Within Adaptation - Explicit')
data_explicit <- data[data$group == "Explicit", ]
adapt_blk_data_explicit=data_explicit[data_explicit$block>0 & data_explicit$block<7,]

model <- lmer(reach_rt ~ coh_cat * block + (1| subject), data = data_explicit)
print(Anova(model, contrasts=list(coh_cat=contr.sum), type=3))

em <- emtrends(model, pairwise ~ coh_cat, var = "block", pbkrtest.limit = 12010)
print(summary(em, infer = c(TRUE, TRUE)))
