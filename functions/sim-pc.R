#' Simulate point count data with removal and distance sampling
#' 
#' @param n_obs How many sampling events in this point count data set?
#' @param phi Singing/cue rate
#' @param tau Effective detection radius
#' @param den Density for bSims
#' @param n_protocols Number of point count protocols to simulate
#' @param n_time_bins Number of time bins per protocol (vector)
#' @param n_dis_bins Number of distance bins per protocol (vector)
#' @param max_times Maximum time for each protocol
#' @param max_dist Maximum distance for each protocol
#' @param n_cores Number of cores for replicates
#' 

sim_pc <- function(n_obs = 1000,
                   phi = NULL,
                   tau = NULL,
                   den = NULL,
                   n_protocols = NULL,
                   n_time_bins = NULL,
                   n_dist_bins = NULL,
                   max_times = NULL,
                   max_dist = NULL,
                   n_cores = NULL)
{
  # Randomly allocate sample size per protocol
  protocol_prop <- round(n_obs * dirmult::rdirichlet(n = 1, alpha = rep(3, n_protocols)))
  protocol_num <- NULL
  for (i in 1:n_protocols)
  {
    protocol_num <- c(protocol_num, rep(i, protocol_prop[i]))
  }

  # Due to rounding, ensure that the total sample size remains == n_obs
  if (length(protocol_num) > n_obs)
  {
    protocol_num <- protocol_num[1:n_obs]
  }else if (length(protocol_num) < n_obs)
  {
    remainder <- n_obs - length(protocol_num)
    to_add <- rep(4, remainder)
    protocol_num <- c(protocol_num, to_add)
  }
  
  rem_names <- c("N", "Protocol",
                 paste0(rep("Int", times = max(n_time_bins)), 1:max(n_time_bins)))
  rem_df_list <- vector(mode = "list", length = n_protocols)
  
  dis_names <- c("N", "Protocol",
                 paste0(rep("Int", times = max(n_dist_bins)), 1:max(n_dist_bins)))
  dis_df_list <- vector(mode = "list", length = n_protocols)
  
  for (p in 1:n_protocols)
  {
    n_obs_p <- unname(table(protocol_num)[p])
    
    rem_df_list[[p]] <- data.frame(as.data.frame(matrix(data = NA,
                                                        nrow = n_obs_p,
                                                        ncol = length(rem_names))))
    names(rem_df_list[[p]]) <- rem_names
    
    dis_df_list[[p]] <- as.data.frame(matrix(data = NA,
                                             nrow = n_obs_p,
                                             ncol = length(dis_names)))
    names(dis_df_list[[p]]) <- dis_names
    
    rem_df_list[[p]]$Protocol <- dis_df_list[[p]]$Protocol <- p
  }

  sim_landscape <- bsims_all(density = den, vocal_rate = phi, tau = tau/100)
  
  #cluster_protocols <- makeCluster(n_cores[1], type = "PSOCK")
  #registerDoParallel(cluster_protocols)
  
  for (p in 1:n_protocols)#each(p = 1:n_protocols, .packages = c('bSims', 'doParallel')) %dopar%
  {
    n_replicates <- nrow(rem_df_list[[p]])
    
    cluster_reps <- makeCluster(n_cores, type = "PSOCK")
    sim_reps <- sim_landscape$replicate(n_replicates, cl = cluster_reps)
    stopCluster(cluster_reps)
    
    tint <- max_times[p, 1:n_time_bins[p]]
    rint <- max_dist[p, 1:n_dist_bins[p]] / 100 
    
    for (n in 1:n_replicates)
    {
      tr <- bsims_transcribe(sim_reps[[n]], tint=tint, rint=rint)
      tally <- tr$removal
      N <- sum(tr$abundance)
      dis_df_list[[p]][n, 2 + c(1:n_dist_bins[p])] <- unname(rowSums(tally))
      rem_df_list[[p]][n, 2 + c(1:n_time_bins[p])] <- unname(colSums(tally))
      
      dis_df_list[[p]][n, "N"] <- rem_df_list[[p]][n, "N"] <- N
    }
  }
  #stopCluster(cluster_protocols)
  
  rem_df <- do.call(rbind, rem_df_list)
  dis_df <- do.call(rbind, dis_df_list)
  
  time_design <- matrix(data = NA,
                        nrow = n_obs,
                        ncol = ncol(max_times))
  time_design <- data.frame(max_times[protocol_num, ])
  names(time_design) <- paste0(rep("Int", times = ncol(max_times)), 1:ncol(max_times))
  
  dist_design <- matrix(data = NA,
                        nrow = n_obs,
                        ncol = ncol(max_dist))
  dist_design <- data.frame(max_dist[protocol_num, ])
  names(dist_design) <- paste0(rep("Int", times = ncol(max_dist)), 1:ncol(max_dist))
  
  return(list(rem = rem_df,
              dis = dis_df,
              dist_design = dist_design,
              time_design = time_design))
}
