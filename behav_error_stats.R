library(lme4)
library(car)
library(emmeans)
library(dplyr)
library(doBy)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(viridisLite)

# Read data
data_path <- "C:/Users/qmoreau/Documents/Work/Beta_bursts/Behavioral/behav_df_cleaned_new.csv"
data <- read.csv(data_path)

# Factorize variables
data$coh_cat <- ordered(data$coh_cat, levels = c('zero', 'low', 'med', 'high'))
data$group <- as.factor(data$group)

# Create block_type variable
data$block_type <- ifelse(data$block == 0, 'baseline',
                          ifelse(data$block > 0 & data$block < 7, 'adaptation',
                                 ifelse(data$block == 7, 'washout', NA)))
data$block_type <- as.factor(data$block_type)


########################## Analysis for Implicit Group ########################## 
data_implicit <- data[data$group == "Implicit", ]
data_implicit$target_error<-data_implicit$reach_vis_err


# Remove outliers above 3 SDs
mean_val <- mean(data_implicit$target_error, na.rm = TRUE)
sd_val <- sd(data_implicit$target_error, na.rm = TRUE)
data_implicit <- data_implicit[abs(data_implicit$target_error - mean_val) <= 3 * sd_val | is.na(data_implicit$target_error), ]

#### Analysis Comparing Block Types for #### 

# Run the LMM
mixed_model_TargetError_block_type_imp <- lmer(target_error ~ block_type + (1| subject), data = data_implicit)
#mixed_model_TargetError_block_type_imp <- lmer(target_error ~ coh_cat + block_type + (1| subject), data = data_implicit)
Anova(mixed_model_TargetError_block_type_imp, contrasts=list(coh_cat=contr.sum, block_type=contr.sum, type=3))

# posthoc
emmeans(mixed_model_TargetError_block_type_imp, pairwise~block_type, type='response', pbkrtest.limit = 15829)
#emmeans(mixed_model_TargetError_block_type_imp, pairwise~coh_cat, type='response', pbkrtest.limit = 15829)


# plotting the block_type main effect
emm_blocktype_imp <- emmeans(mixed_model_TargetError_block_type_imp, ~ block_type, pbkrtest.limit = 15829)
emm_blocktype_imp_df <- as.data.frame(emm_blocktype_imp)
emm_blocktype_imp_df$block_type <- factor(emm_blocktype_imp_df$block_type, levels = c("baseline", "adaptation", "washout"))

# Define palette
palette1 <- magma(6)[3:5]

# Plot block_type main effect
p1 <- ggplot(emm_blocktype_imp_df, aes(x = block_type, y = emmean, ymin = emmean - SE, ymax = emmean + SE)) +
  geom_point(size = 6, aes(fill = block_type, color = block_type)) +
  geom_errorbar(width = 0.2, aes(color = block_type)) +
  labs(x = "Block", y = "Target Error", title = "Target Error Analysis per Block Type") +
  theme_minimal() +
  scale_fill_manual(values = palette1) +
  scale_color_manual(values = palette1)
print(p1 + theme(legend.position = "none"))

# plotting the coh_cat main effect
emm_coh_cat_imp <- emmeans(mixed_model_TargetError_block_type_imp, ~ coh_cat, pbkrtest.limit = 15829)
emm_coh_cat_imp_df <- as.data.frame(emm_coh_cat_imp)
emm_coh_cat_imp_df$coh_cat <- factor(emm_coh_cat_imp_df$coh_cat)

palette1 <- magma(6)[3:6]

p2 <- ggplot(emm_coh_cat_imp_df, aes(x = coh_cat, y = emmean, ymin = emmean - SE, ymax = emmean + SE)) +
  geom_point(size = 6, aes(fill = coh_cat, color = coh_cat)) +
  geom_errorbar(width = 0.2, aes(color = coh_cat)) +
  labs(x = "coh_cat", y = "Target Error", title = "Target Error Analysis per coh_cat") +
  theme_minimal() +
  scale_fill_manual(values = palette1) +
  scale_color_manual(values = palette1)
print(p2 + theme(legend.position = "none"))


#### Analysis Comparing Blocks Within Adaptation  ####
adapt_blk_data_implicit=data_implicit[data_implicit$block>0 & data_implicit$block<7,]

# Run the LMM for adaptation blocks
mixed_model_TargetError_adapt_block_imp <- lmer(target_error ~ coh_cat * block + (1| subject), data = adapt_blk_data_implicit)

Anova(mixed_model_TargetError_adapt_block_imp, contrasts=list(coh_cat=contr.sum, block=contr.sum, type=3))

# posthoc
emmeans(mixed_model_TargetError_adapt_block_imp, pairwise~block, type='response', pbkrtest.limit = 15829)
emmeans(mixed_model_TargetError_adapt_block_imp, pairwise~coh_cat, type='response', pbkrtest.limit = 15829)

# Add predicted values to adaptation data
adapt_blk_data_implicit$predicted_targeterror <- predict(mixed_model_TargetError_adapt_block_imp, newdata = adapt_blk_data_implicit, type = "response")


# Scatterplot of target_error vs. block
p2 <- ggplot(adapt_blk_data_implicit, aes(x = block, y = target_error)) +
  geom_point(alpha = 0.1) +
  geom_smooth(aes(y = predicted_targeterror), method = "lm", se = FALSE, linewidth = 2, color = 'black') +
  labs(x = "Block", y = "Target Error") +
  theme(panel.background = element_blank())  # Remove background
ggtitle("Scatterplot of Target Error vs. Block")
print(p2)


########################## Analysis for Explicit Group ########################## 
data_explicit <- data[data$group == "Explicit", ]
data_explicit$error_magnitude<-data_explicit$reach_vis_abs_err

# Compute the mean and standard deviation, ignoring NA values
mean_val <- mean(data_explicit$error_magnitude, na.rm = TRUE)
sd_val <- sd(data_explicit$error_magnitude, na.rm = TRUE)
# Remove outliers above 3 SDs
data_explicit <- data_explicit[abs(data_explicit$error_magnitude - mean_val) <= 3 * sd_val | is.na(data_explicit$error_magnitude), ]

#### Analysis Comparing Block Types for #### 

# Run the LMM
mixed_model_ErrorMagnitude_block_type_exp <- lmer(error_magnitude ~ block_type + (1| subject), data = data_explicit)
# mixed_model_ErrorMagnitude_block_type_exp <- lmer(error_magnitude ~ coh_cat + block_type + (1| subject), data = data_explicit)

Anova(mixed_model_ErrorMagnitude_block_type_exp, type=3)

# posthoc
emmeans(mixed_model_ErrorMagnitude_block_type_exp,pairwise~block_type, type='response', pbkrtest.limit = 15829)
# emmeans(mixed_model_ErrorMagnitude_block_type_exp,pairwise~coh_cat, type='response')

# plotting the block_type main effect
emm_blocktype_exp <- emmeans(mixed_model_ErrorMagnitude_block_type_exp, ~ block_type, pbkrtest.limit = 15829)
emm_blocktype_exp_df <- as.data.frame(emm_blocktype_exp)
emm_blocktype_exp_df$block_type <- factor(emm_blocktype_exp_df$block_type, levels = c("baseline", "adaptation", "washout"))

palette3 <- viridis(6)
palette3 <- palette3[3:5]

p3 <- ggplot(emm_blocktype_exp_df, aes(x = block_type, y = emmean, ymin = emmean - SE, ymax = emmean + SE)) +
  geom_point(size = 6, aes(fill = block_type, color = block_type)) +
  geom_errorbar(width = 0.2, aes(color = block_type)) +
  labs(x = "Block", y = "Error Magnitude", title = "Target Error Analysis per Block Type") +
  theme_minimal() +
  # facet_wrap(~ reorder(group, -as.numeric(group))) + # Reordering based on the numeric value of 'group'
  scale_fill_manual(values = palette3) +
  scale_color_manual(values = palette3)
print(p3 + theme(legend.position = "none"))

# plotting the coh_cat main effect
# emm_coh_cat_exp <- emmeans(mixed_model_ErrorMagnitude_block_type_exp, ~ coh_cat, pbkrtest.limit = 15829)
# emm_coh_cat_exp_df <- as.data.frame(emm_coh_cat_exp)
# emm_coh_cat_exp_df$coh_cat <- factor(emm_coh_cat_exp_df$coh_cat, levels = c("zero", "low", "med", "high"))

# posthoc
pairwise_comp_exp <- contrast(emm_coh_cat_exp, method = "pairwise", adjust = "tukey")
print(pairwise_comp_exp)

palette4 <- viridis(6)
palette4 <- palette4[3:6]
# 
# p4 <- ggplot(emm_coh_cat_exp_df, aes(x = coh_cat, y = emmean, ymin = emmean - SE, ymax = emmean + SE)) +
#   geom_point(size = 6, aes(fill = coh_cat, color = coh_cat)) +
#   geom_errorbar(width = 0.2, aes(color = coh_cat)) +
#   labs(x = "Coherence Category", y = "Error Magnitude", title = "Target Error Analysis per Block Type") +
#   theme_minimal() +
#   # facet_wrap(~ reorder(group, -as.numeric(group))) + # Reordering based on the numeric value of 'group'
#   scale_fill_manual(values = palette4) +
#   scale_color_manual(values = palette4)
# print(p4 + theme(legend.position = "none"))


#### Analysis Comparing Blocks Within Adaptation  ####
adapt_blk_data_explicit=data_explicit[data_explicit$block>0 & data_explicit$block<7,]

mixed_model_ErrorMagnitude_adapt_block_exp <- lmer(error_magnitude ~ coh_cat * block + (1| subject), data = adapt_blk_data_explicit)
Anova(mixed_model_ErrorMagnitude_adapt_block_exp)

interaction_trends <- emtrends(mixed_model_ErrorMagnitude_adapt_block_exp, pairwise ~ coh_cat * block, var = "block")
print(summary(interaction_trends, infer = c(TRUE, TRUE)))

adapt_blk_data_explicit$predicted_error_magnitude <- predict(mixed_model_ErrorMagnitude_adapt_block_exp, newdata = adapt_blk_data_explicit, type = "response")

# Scatterplot of reach_rt vs. block, colored by group
p5 <- ggplot(adapt_blk_data_explicit, aes(x = block, y = error_magnitude, color = coh_cat)) +
  geom_point(alpha = 0.1) +
  geom_smooth(aes(y = predicted_error_magnitude), method = "lm", se = FALSE, linewidth = 2) +
  labs(x = "Block", y = "Error Magnitude") +
  theme(panel.background = element_blank())  # Remove background
print(p5)
