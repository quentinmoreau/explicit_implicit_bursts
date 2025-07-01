run_trial_permutation <- function(varname, use_next_trial = TRUE) {
  library(permutes)
  library(tidyverse)
  library(doParallel)
  
  source("utils.R")
  options(contrasts = c("contr.sum", "contr.poly"))
  
  cl <- makeCluster(8, outfile = '')
  registerDoParallel(cl)
  
  for (epoch in c("vis", "mot")) {
    all_results <- list()
    
    for (pc in 7:7) {
      df <- read_csv(paste0("output/PC_", pc, "_", epoch, "_trial_burst_counts.csv")) %>%
        filter(block > 0 & block < 7)
      
      window_cols <- grep("^time_", names(df), value = TRUE)
      window_cols <- window_cols[3:(length(window_cols) - 2)]
      
      for (q in 1:4) {
        df_q <- df %>% 
          filter(quartile == q)
        
        if (use_next_trial) {
          df_q <- df_q %>% add_next_trial_value(varname)
          value_col <- paste0("next_", varname)
        } else {
          value_col <- varname
        }
        
        df_long <- df_q %>%
          select(subject, group, block, all_of(value_col), all_of(window_cols)) %>%
          pivot_longer(cols = all_of(window_cols), names_to = "time", values_to = "burst_count") %>%
          mutate(time = as.numeric(gsub("time_", "", time)))
        
        formula <- as.formula(paste0(value_col, " ~ group * burst_count + (1 | subject)"))
        
        perms <- clusterperm.lmer(
          formula,
          data = df_long,
          series.var = ~time,
          parallel = TRUE,
          progress = '',
          type = 'anova'
        )
        
        perms$pc <- pc
        perms$quartile <- q
        all_results[[paste0("PC", pc, "_Q", q)]] <- perms
      }
    }
    
    filename_suffix <- if (use_next_trial) paste0(varname, "_next_trial") else varname
    results <- bind_rows(all_results)
    write_csv(results, paste0("output/PC_", epoch, "_", filename_suffix, "_perm.csv"))
  }
}

run_burst_rate_trend_analysis <- function(varname, use_next_trial = TRUE) {
  library(lme4)
  library(lmerTest)
  library(tidyverse)
  library(scales)
  source("utils.R")
  
  for (epoch in c("vis", "mot")) {
    all_trends <- list()
    
    for (pc in 7:7) {
      df <- read_csv(paste0("output/PC_", pc, "_", epoch, "_trial_burst_counts.csv")) %>%
        filter(block > 0 & block < 7)
      
      window_cols <- grep("^time_", names(df), value = TRUE)
      window_cols <- window_cols[3:(length(window_cols) - 2)]
      
      for (q in 1:4) {
        df_q <- df %>% filter(quartile == q)
        
        if (use_next_trial) {
          df_q <- df_q %>% add_next_trial_value(varname)
          value_col <- paste0("next_", varname)
        } else {
          value_col <- varname
        }
        
        trends <- list()
        
        for (col in window_cols) {
          df_q$burst_count <- df_q[[col]]
          time_val <- as.numeric(gsub("time_", "", col))
          
          model <- lmer(
            formula = as.formula(paste0(value_col, " ~ group * burst_count + (1 | subject)")),
            data = df_q
          )
          
          tr_group <- emtrends(model, ~group, var = "burst_count", infer = TRUE)
          tr_df <- summary(tr_group) %>%
            as_tibble() %>%
            mutate(
              time = time_val,
              quartile = q,
              pc = pc
            )
          
          trends[[col]] <- tr_df
        }
        
        trends_df <- bind_rows(trends) %>%
          rename(estimate = burst_count.trend)
        
        all_trends[[paste0("PC", pc, "_Q", q)]] <- trends_df
      }
    }
    
    results <- bind_rows(all_trends)
    
    suffix <- if (use_next_trial) paste0("burst_rate_", varname, "_next_trial") else paste0("burst_rate_", varname)
    write_csv(results, paste0("output/PC_", epoch, "_", suffix, ".csv"))
  }
}


vars<-c('aim_vis_abs_err', 'reach_vis_abs_err', 'reach_dur', 'reach_rt')

for(var in vars) {
  run_trial_permutation(var, use_next_trial = TRUE)
  run_burst_rate_trend_analysis(var, use_next_trial = TRUE)
  run_trial_permutation(var, use_next_trial = FALSE)
  run_burst_rate_trend_analysis(var, use_next_trial = FALSE)
}

