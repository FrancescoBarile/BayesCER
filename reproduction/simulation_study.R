library(doParallel)
library(foreach)
library(BayesCER)

# DGP
n_vec = c( 10, 100, 1000)
N_vec = c( 10, 100, 1000, 2000)
alpha_vec = c(0.1, 0.2, 0.3)
scenario_mat = expand.grid(N_vec, alpha_vec)

# hyperparameters
W_0 = 1/2
n_0 = 1
a = b = 1
K = 2

# mcmcm options
n_iter=1000
burnin=1000
thinning=10
mcmc_options=list( n_iter=n_iter, burnin=burnin, thinning=thinning  )

N_rep = 100


set.seed(1)
all_seeds <- sample.int(
  .Machine$integer.max,
  nrow(scenario_mat)
)

set.seed(2)
seed_array <- array(
  sample.int(.Machine$integer.max,
             nrow(scenario_mat) * length(n_vec) * N_rep),
  dim = c(nrow(scenario_mat), length(n_vec), N_rep)
)

n_cores = max( (1:detectCores())[ N_rep %% (1:detectCores())==0 ] )


project_path = getwd()

out_path = paste0(project_path, "/data")
if (!dir.exists( out_path )) {
  dir.create(out_path)
}

folder_out <- paste0(out_path, "/output")
if (!dir.exists( folder_out )) {
  dir.create(folder_out)
}

myCluster <- makeCluster(n_cores, type = "PSOCK")

st_all <- Sys.time()
registerDoParallel(myCluster)
for( s in 1:nrow(scenario_mat)){

  set.seed(all_seeds[s])

  N = scenario_mat[s,1]
  M = choose(N,2)

  W=rep(W_0, M)

  alpha = scenario_mat[s,2]

  prior_list=list(W=W, n_0=n_0, a=a, b=b)

  # DGP centroid
  C_mat = BayesCER::generate_sbm_fast(N, K = 2,
                                      intra_prob = 0.8, inter_prob = 0.1,
                            z = sample(seq_len(K), N, TRUE))
  C = BayesCER:::from_full_adj_to_low_adj(C_mat)

  # ECC DGP distribution
  summary_C = BayesCER::exact_graph_stats(C)

  ecc_true = clustering(alpha, edges = summary_C[1], d2 = summary_C[2],  summary_C[3], N)
  density_true = (M * alpha + sum(C) * ( 1 - 2*alpha ) )/M

  for( n in n_vec ){

    folder <- paste0(folder_out,
                     "/n_", n,
                     "_N_",  N,
                     "_alpha_", alpha)


    if (!dir.exists( folder )) {
      dir.create(folder)
    }

    foreach(r = 1:N_rep,
            .packages = c("Rcpp","BayesCER")
    ) %dopar% {

      seed_snr = seed_array[s, which(n_vec == n), r]
      set.seed(seed_snr)

      # generate synthetic data
      st <- Sys.time()
      G_sum <- as.integer(rbinom(M,n,alpha+(1-2*alpha)*C))
      G_bar <- G_sum/n
      en <- Sys.time()
      time_start_rgen = st
      time_end_rgen   = en
      time_rgen = c(time_start_rgen, time_end_rgen)

      # Inference
      # EB inference
      st <- Sys.time()
      eb_infer = gCER_inference(G_sum, n, prior_list=prior_list,
                                mcmc=F,
                                parallel=F)
      en <- Sys.time()

      time_start_EB = st
      time_end_EB   = en
      time_EB = c(time_start_EB, time_end_EB)

      C_hat = eb_infer$C_hat
      alpha_hat = eb_infer$alpha_hat

      discrepancy_C_hat_EB     = sum(abs(C_hat-C))/M
      discrepancy_alpha_hat_EB = abs(alpha_hat-alpha)

      ecc_hat_EB      = eb_infer$ecc
      density_EB = eb_infer$density

      #Gibbs
      st <- Sys.time()
      gibbs_infer = gCER_inference(G_sum, n, prior_list=prior_list, mcmc=T,
                                   mcmc_options=mcmc_options,
                                   parallel=F,
                                   seed = seed_snr)
      en <- Sys.time()

      time_start_gibbs = st
      time_end_gibbs   = en
      time_gibbs       = c(time_start_gibbs, time_end_gibbs)

      C_gibbs = as.integer(gibbs_infer$C_mean>=.5)
      alpha_gibbs = mean(gibbs_infer$alpha_chain)
      discrepancy_C_hat_gibbs     = sum(abs(C_gibbs-C))/M
      discrepancy_alpha_hat_gibbs = abs(alpha_gibbs-alpha)

      ecc_hat_gibbs = mean(gibbs_infer$ecc_dist)
      density_gibbs = mean(gibbs_infer$density_dist)

      temp = list( n = n, N=N, alpha=alpha, rep = r,

                   ecc_true=ecc_true,
                   density_true=density_true,

                   time_EB=time_EB, time_gibbs=time_gibbs, time_rgen=time_rgen,

                   discrepancy_C_hat_EB=discrepancy_C_hat_EB,
                   discrepancy_alpha_hat_EB=discrepancy_alpha_hat_EB,

                   discrepancy_C_hat_gibbs=discrepancy_C_hat_gibbs,
                   discrepancy_alpha_hat_gibbs=discrepancy_alpha_hat_gibbs,

                   ecc_hat_EB=ecc_hat_EB,
                   ecc_hat_gibbs=ecc_hat_gibbs,

                   density_EB=density_EB,
                   density_gibbs=density_gibbs
      )

      saveRDS(temp, paste0(folder, "/n_", n,
                           "_N_",  N, "_alpha_", alpha, "_rep_", r, ".rds" ) )
      rm(temp)
      gc()

      print(paste0("n_", n,
                   "_N_",  N,
                   "_alpha_", alpha,
                   "_rep_", r ))

    }
  }
}
stopCluster(myCluster)
en_all <- Sys.time()
en_all - st_all

output_path = folder_out

output_folders = list.files(output_path)

output_df  = data.frame()
for( fo in 1:length(output_folders) ){
  folder_fo = list.files(paste0(output_path, "/", output_folders[fo]))
  for (fi in 1:length(folder_fo) ){
    output_fi = readRDS(paste0(output_path, "/",
                               output_folders[fo], "/",
                               folder_fo[fi]  ))
    temp_df = data.frame(n                   = output_fi$n,
                         N                   = output_fi$N,
                         alpha               = output_fi$alpha,
                         rep                 = output_fi$rep,
                         starting_time_EB    = output_fi$time_EB[1],
                         ending_time_EB      = output_fi$time_EB[2],
                         starting_time_gibbs = output_fi$time_gibbs[1],
                         ending_time_gibbs   = output_fi$time_gibbs[2],
                         starting_time_rgen  = output_fi$time_rgen[1],
                         ending_time_rgen    = output_fi$time_rgen[2],
                         ecc_true            = output_fi$ecc_true,
                         density_true                 = output_fi$density_true,
                         discrepancy_C_hat_EB         = output_fi$discrepancy_C_hat_EB,
                         discrepancy_C_hat_gibbs      = output_fi$discrepancy_C_hat_gibbs,
                         discrepancy_alpha_hat_EB     = output_fi$discrepancy_alpha_hat_EB,
                         discrepancy_alpha_hat_gibbs  = output_fi$discrepancy_alpha_hat_gibbs,
                         ecc_hat_EB                   = output_fi$ecc_hat_EB,
                         ecc_hat_gibbs                = output_fi$ecc_hat_gibbs,
                         density_EB                   = output_fi$density_EB,
                         density_gibbs                = output_fi$density_gibbs
    )

    output_df = rbind(output_df,  temp_df)

  }
}


library(tidyverse)
output_df_long <- output_df |>
  pivot_longer(
    cols = -c(n , N, alpha, rep, ecc_true, density_true,
              starting_time_rgen, ending_time_rgen
    ),
    names_to = c(".value", "method_temp"),
    names_pattern = "(.*)_(EB|gibbs)"
  ) %>% mutate(method=case_when(method_temp == "gibbs" ~ "Gibbs",
                                TRUE ~ method_temp)) %>%
  select(-method_temp) %>%
  relocate(method, .after = ending_time)





