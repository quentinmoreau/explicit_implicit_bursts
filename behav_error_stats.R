library(lme4)
library(car)
library(emmeans)
library(dplyr)
library(doBy)

# Read data
data_path <- "data/behavioral/behav_df_cleaned_new.csv"
data <- read.csv(data_path)

# Factorize variables
data$coh_cat <- ordered(data$coh_cat, levels = c('zero', 'low', 'med', 'high'))
data$group <- as.factor(data$group)

# Create block_type variable
data$block_type <- ifelse(data$block == 0, 'baseline',
                          ifelse(data$block > 0 & data$block < 7, 'adaptation',
                                 ifelse(data$block == 7, 'washout', NA)))
data$block_type <- as.factor(data$block_type)

########################## Absolute error magnitude #############################
model <- lmer(reach_vis_abs_err ~ group*block_type + (1| subject), data = data)
Anova(model, contrasts=list(group=contr.sum, block_type=contr.sum, type=3))
emmeans(model, pairwise~block_type|group, type='response', pbkrtest.limit = 15987)

########################## Signed error - Implicit Group ########################## 
data_implicit <- data[data$group == "Implicit", ]

model <- lmer(reach_vis_err ~ block_type + (1| subject), data = data_implicit)
Anova(model, contrasts=list(block_type=contr.sum, type=3))

# posthoc
emmeans(model, pairwise~block_type, type='response', pbkrtest.limit = 15829)


########################## Signed error - Implicit Group within Adaptation ########################## 
adapt_blk_data_implicit=data_implicit[data_implicit$block>0 & data_implicit$block<7,]

# Run the LMM for adaptation blocks
model <- lmer(reach_vis_err ~ block + (1| subject), data = adapt_blk_data_implicit)

Anova(model)


########################## Analysis for Coherence and Error magnitude ########################## 
adapt_blk_data_implicit=data_implicit[data_implicit$block>0 & data_implicit$block<7,]
model <- lmer(reach_vis_abs_err ~ coh_cat * block + (1| subject), data = adapt_blk_data_implicit)
Anova(model)

data_explicit <- data[data$group == "Explicit", ]
adapt_blk_data_explicit=data_explicit[data_explicit$block>0 & data_explicit$block<7,]
model <- lmer(reach_vis_abs_err ~ coh_cat * block + (1| subject), data = adapt_blk_data_explicit)
Anova(model)

em <- emtrends(model, pairwise ~ coh_cat * block, var = "block", pbkrtest.limit = 12010)
print(summary(em, infer = c(TRUE, TRUE)))

