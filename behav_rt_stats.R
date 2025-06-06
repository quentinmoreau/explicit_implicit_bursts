library(lme4)
library(car)
library(emmeans)
library(dplyr)

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

############################ Analysis Comparing Block Types ############################  
# Run the LMM
model <- lmer(reach_rt ~ block_type * group + (1| subject), data = data)
Anova(model, contrasts=list(group=contr.sum, block_type=contr.sum), type=3)
emmeans(model, pairwise ~ block_type|group, type='response', pbkrtest.limit = 15829)
emmeans(model, pairwise ~ group|block_type, type='response', pbkrtest.limit = 15829)

############################ Analysis Comparing Blocks Within Adaptation  ############################  
# subset the data only for adaptation blocks
data_implicit <- data[data$group == "Implicit", ]
adapt_blk_data_implicit=data_implicit[data_implicit$block>0 & data_implicit$block<7,]

# Run the LMM for adaptation blocks
model <- lmer(reach_rt ~ coh_cat * block + (1| subject), data = adapt_blk_data_implicit)
Anova(model, contrasts=list(coh_cat=contr.sum), type=3)

data_explicit <- data[data$group == "Explicit", ]
adapt_blk_data_explicit=data_explicit[data_explicit$block>0 & data_explicit$block<7,]

# Run the LMM for adaptation blocks
model <- lmer(reach_rt ~ coh_cat * block + (1| subject), data = data_explicit)
Anova(model, contrasts=list(coh_cat=contr.sum), type=3)

em <- emtrends(model, pairwise ~ coh_cat * block, var = "block", pbkrtest.limit = 12010)
print(summary(em, infer = c(TRUE, TRUE)))
