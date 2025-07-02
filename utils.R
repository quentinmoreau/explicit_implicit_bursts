assign_cluster_ids <- function(sig_vec, min_size = 4) {
  r <- rle(sig_vec)
  cluster_id <- rep(NA, length(sig_vec))
  cluster_counter <- 0
  idx <- 1
  for (i in seq_along(r$lengths)) {
    len <- r$lengths[i]
    val <- r$values[i]
    if (val && len >= min_size) {
      cluster_counter <- cluster_counter + 1
      cluster_id[idx:(idx + len - 1)] <- cluster_counter
    }
    idx <- idx + len
  }
  return(cluster_id)
}

add_next_trial_value <- function(df, colname) {
  df <- df %>%
    arrange(subject, block, trial) %>%
    mutate(next_trial = trial + 1)
  
  next_df <- df %>%
    select(subject, block, trial, !!sym(colname)) %>%
    rename(
      next_trial = trial,
      next_value = !!sym(colname)
    )
  
  df_joined <- left_join(df, next_df, by = c("subject", "block", "next_trial"))
  
  # Rename the joined column to reflect the input
  new_colname <- paste0("next_", colname)
  df_joined <- df_joined %>%
    rename(!!new_colname := next_value)
  
  return(df_joined)
}