library(lme4)
library(car)
library(emmeans)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(viridisLite)
source('C:/Users/qmoreau/Documents/Work/Beta_bursts/Bursts/plot_fit.R')

# Read data
data <-  read.csv('C:/Users/qmoreau/Documents/Work/Beta_bursts/Behavioral/behav_df_cleaned_new.csv')

# Factorize variables

data$coh_cat <- factor(data$coh_cat, levels = c('zero','low','med','high'), ordered = TRUE)
data$group<-as.factor(data$group)

############################ Analysis Comparing Block Types ############################  

# Create block type variable
data$block_type <- ifelse(data$block == 0, 'baseline',
                          ifelse(data$block > 0 & data$block < 7, 'adaptation',
                                 ifelse(data$block == 7, 'washout', NA)))
data$block_type <- as.factor(data$block_type)

# remove outliers
mean_val <- mean(data$reach_rt, na.rm = TRUE)
sd_val <- sd(data$reach_rt, na.rm = TRUE)
data <- data[abs(data$reach_rt- mean_val) <= 3 * sd_val | is.na(data$reach_rt), ]

# Run the LMM
mixed_model_RT_block_type <- lmer(reach_rt ~ block_type * group + (1| subject), data = data)
# mixed_model_RT_block_type <- lmer(reach_rt ~ coh_cat * group + block_type * group + (1| subject), data = data)
Anova(mixed_model_RT_block_type, contrasts=list(coh_cat=contr.sum, group=contr.sum, block_type=contr.sum), type=3)

# plotting the group:block_type interaction
emm_group_blocktype <- emmeans(mixed_model_RT_block_type, ~ group * block_type, pbkrtest.limit = 15829)
emm_group_blocktype_df <- as.data.frame(emm_group_blocktype)
emm_group_blocktype_df$block_type <- factor(emm_group_blocktype_df$block_type, levels = c("baseline", "adaptation", "washout"))

# posthoc
emmeans(mixed_model_RT_block_type,pairwise~group|block_type, type='response')
emmeans(mixed_model_RT_block_type,pairwise~block_type|group, type='response')

# Define palettes

palette1 <- magma(6)
palette1 <- palette1[3:5]
palette2 <- mako(7)
palette2 <- palette2[4:6]

# Create interaction variable
emm_group_blocktype_df$block_type_group <- interaction(emm_group_blocktype_df$block_type, emm_group_blocktype_df$group)


p1 <- ggplot(emm_group_blocktype_df, aes(x = block_type, y = emmean, ymin = emmean - SE, ymax = emmean + SE)) +
  geom_point(size = 6, aes(fill = block_type_group, color = block_type_group)) +
  geom_errorbar(width = 0.2, aes(color = block_type_group)) +
  labs(x = "Block", y = "RTs", title = "Reaction Times Analysis per Block Type") +
  theme_minimal() +
  facet_wrap(~ reorder(group, -as.numeric(group))) + # Reordering based on the numeric value of 'group'
  scale_fill_manual(values = c(palette2, palette1)) +
  scale_color_manual(values = c(palette2, palette1))
print(p1 + theme(legend.position = "none"))


############################ Analysis Comparing Blocks Within Adaptation  ############################  
# subset the data only for adaptation blocks
adapt_blk_data=data[data$block>0 & data$block<7,]

# Run the LMM for adaptation blocks
mixed_model_RT_adapt_block <- lmer(reach_rt ~ coh_cat * group * block + (1| subject), data = adapt_blk_data)
Anova(mixed_model_RT_adapt_block, contrasts=list(coh_cat=contr.sum, group=contr.sum), type=3)

# posthoc
interaction_trends <- emtrends(mixed_model_RT_adapt_block, pairwise ~ group, var = "block")
print(summary(interaction_trends, infer = c(TRUE, TRUE)))

# Add predicted values to adaptation data
adapt_blk_data$predicted_reach_rt <- predict(mixed_model_RT_adapt_block, newdata = adapt_blk_data, type = "response")

# Scatterplot of reach_rt vs. block, colored by group
p2 <- ggplot(adapt_blk_data, aes(x = block, y = reach_rt, color = group)) +
  geom_point(alpha = 0.1) +
  geom_smooth(aes(y = predicted_reach_rt), method = "lm", se = FALSE, linewidth = 2) +
  labs(x = "Block", y = "Reaction Time (ms)") +
  theme(panel.background = element_blank())  # Remove background
ggtitle("Reach_rt within Adaptation blocks")
print(p2)
 
