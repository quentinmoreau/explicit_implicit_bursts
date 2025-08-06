library('lme4')
library('car')
library('emmeans')

for(epo in c('vis','mot')) {
    print(paste0('===== Epoch: ', epo, '====='))
    df <- read.csv(paste0('./output/burst_features_',epo,'.csv'))

    df$group <- as.factor(df$group)
    contrasts(df$group) <- contr.sum
      
    m <- lmer(fwhm_time ~ group + (1 | subject), data = df)
    print(Anova(m, contrasts=list(group=contr.sum)))

    m <- lmer(peak_amp_base ~ group + (1 | subject), data = df)
    print(Anova(m, contrasts=list(group=contr.sum)))

    m <- lmer(fwhm_freq ~ group + (1 | subject), data = df)
    print(Anova(m, contrasts=list(group=contr.sum)))

    m <- lmer(peak_freq ~ group + (1 | subject), data = df)
    print(Anova(m, contrasts=list(group=contr.sum)))
}