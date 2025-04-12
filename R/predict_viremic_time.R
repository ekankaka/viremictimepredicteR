
#' Thus Function uses posterior samples from a markov chain monte carlo (mcmc) 
#' estimation to predict a viremic time and credible interval 
#' given a single diversity metric and specification of weights 
#' based on the count of eligible unique sequences
#' 
#' @param mcmc_samples. mcmc list
#' 
#' @param x_new. Numeric. New diversity estimate to use for predicting viremic time.
#' Currently supported diversity metrics include:
#' - mean pairwise distance from the "raw" model (rawMPD)
#' - mean pairwise distance from the "TN93" model (tn93MPD)
#' - nucleotide diversity from the "raw" model (rawPI)
#' - nucleotide diversity from the "TN93" model (tn93PI)
#' - weighted fraction of polymorphic sites (WFPS)
#' - WFPS at third codon positions only (WFPS_codons)
#' @param weight_new. Numeric. Weight to use during estimation.
#' Options include:
#' 1 (None)
#' uniqueseqs (raw count of unique sequences)
#' log(uniqueseqs) (log-transformed counts of unique sequences)
#'
#' @return Retuns list of median viremic time, 
#' and lower and upper bounds of credible interval.
#' 
#' @examples
#' res <- predict_y(mcmc_samples = mymcmcobject, x_new = mynewx, weight_new = mynewweight)
#'
#' @export
predict_y <- function(mcmc_samples, x_new, weight_new){
  
  a_samples <- mcmc_samples[,"a"]
  b_samples <- mcmc_samples[,"b"]
  sig_samples = mcmc_samples[,"sig"]
  prec_sig_samples <- 1/sig_samples
  
  # Compute predictive mean for new x
  mu_pred <- a_samples + b_samples * x_new
  
  # Compute scaled precision for the new point
  prec_pred <- weight_new * prec_sig_samples
  
  # Simulate from predictive distribution (e.g., posterior predictive draws)
  y_pred_draws <- rnorm(length(mu_pred), mean = mu_pred, sd = sqrt(1 / prec_pred))
  
  # Compute mean and credible interval (95%)
  summary_stats <- c(
    mean = mean(y_pred_draws),
    sd = sd(y_pred_draws),
    quantile(y_pred_draws, probs = c(0.025, 0.5, 0.975))
  )
  
  return(list(median = summary_stats[4], 
           lwr = summary_stats[3], 
           upr = summary_stats[5]))
}
#'
#'
#'
#'
#'#' This Function predicts viremic time (and credible interval)
#' for multiple weighting methods and diversity metrics.
#'
#' Based on imput, the function selects markov chain monte carlo posterior samples to use,
#' for the specified HIV region, weighting method, and diversity metric.
#' 
#' @param distances. List. Result from the function calculate_distance.
#'
#' @param hiv_region. Character. Currently supported HIV hiv_regions include:
#' - Matrix (especially p17 portion of Matrix in gag)
#' - RT (reverse transcriptase region of pol)
#' - gp41 (gp41 region of env)
#'
#' @param weights_type. Character. Choose which weights to use for estimation. 
#' Currently supported weights include: 
#' - None
#' - UniqueseqsAsis
#' - UniqueseqsLogTransformed (Recommended)
#' With default settings, result contain estimates for each weighting method.
#'
#' @param diversity_metrics. Character. 
#' Choose which diversity metrics to use for estimation
#' Currently supported diversity metrics include:
#' - mean pairwise distance from the "raw" model (rawMPD)
#' - mean pairwise distance from the "TN93" model (tn93MPD)
#' - nucleotide diversity from the "raw" model (rawPI)
#' - nucleotide diversity from the "TN93" model (tn93PI)
#' - weighted fraction of polymorphic sites (WFPS)
#' - WFPS at third codon positions only (WFPScodons)
#' By default, result contain estimates for each metric.
#'
#' @return Returns dataframe of input specifications 
#' and viremic time estimates for each specification
#'
#' @examples
#' viremic_time <- predict_viremic_time(distances = dist_list, hiv_region = "gp41")
#'
#' @export
predict_viremic_time <- function(
    distances,
    hiv_region,
    weights_type = c("None", "UniqueseqsAsis", "UniqueseqsLogTransformed"),
    diversity_metrics = c("rawMPD", "tn93MPD", "rawPI", "tn93PI", "WFPS", "WFPScodons")){
  
  # sequence type
  sequence_type = distances$sequence_type
  
  # unique sequences
  uniqueseqs = distances$uniqueseqs
  
  # alignment width
  alignmentwidth = distances$alignmentwidth
  
  # check if hiv region is eligible for prediction
  hiv_region_formatted <- case_when(
    grepl("matrix", tolower(hiv_region)) ~ "Matrix",
    grepl("p17", tolower(hiv_region)) ~ "Matrix",
    grepl("rt", tolower(hiv_region)) ~ "RT",
    grepl("gp41", tolower(hiv_region)) ~ "gp41",
    TRUE ~ NA)
  
  eligible_regions = c("Matrix", "RT", "gp41")
  if (!(tolower(hiv_region_formatted) %in% tolower(eligible_regions))) stop("hiv_region must be one of: Matrix, RT, or gp41")
  
  # read mcmc samples
  samples = readRDS("data/samples.rds")
  

  
  # select weights to use for prediction
  predictions = data.frame(sequence_type = NA, hiv_region = NA, uniqueseqs = NA,
                           alignmentwidth = NA, weights_type = NA, 
                           diversity_metric = NA, median = NA, lwr = NA, upr = NA)
  
  for (i in 1:length(weights_type)){
    for (j in 1:length(diversity_metrics)){
      print("########################################")
      
      # prediction model to use for prediction
      model_name = paste(sequence_type, 
                         hiv_region_formatted, 
                         weights_type[i],
                         diversity_metrics[j], 
                         sep = "_")
      
      print(model_name)
      
      # mcmc samples for this model
      mod_sim = samples[[model_name]]
      mod_csim = as.mcmc(do.call(rbind, mod_sim))
      
      # new x value (diversity) to use for prediction
      x_new = distances[[diversity_metrics[j]]]
      
      # weight to use in prediction
      weight = case_when(
        weights_type == "None" ~ 1,
        weights_type == "UniqueseqsAsis" ~ uniqueseqs,
        weights_type == "UniqueseqsLogTransformed" ~ log(uniqueseqs),
        TRUE ~ NA)
      
      # predict
      res <- predict_y(mcmc_samples = mod_csim, x_new = x_new, weight_new = weight)
      out <- c(sequence_type = sequence_type, hiv_region = hiv_region_formatted, 
               uniqueseqs = uniqueseqs, alignmentwidth = alignmentwidth,
               weights_type = weights_type[i], 
               diversity_metric = diversity_metrics[j], res)
      
      predictions <- data.frame(rbind(predictions, out)) %>% filter(!is.na(hiv_region))
      print(predictions[nrow(predictions),])
      
      }
      
  }
  
  return(predictions)
  
  }
