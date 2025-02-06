
library(ape)
library(phangorn)
library(pegas)
library(Biostrings)
library(strex) # for str_elem
library(rjags)

#' Function to predict time since infection given a distance
#'
#' This function selects markov chain monte carlo posterior samples,
#' for the specified HIV region and distance metric.
#' @param hiv_region. Character. Currently supported HIV hiv_regions (selected due to better performance) include:
#' - Envgp120V5
#' - Polp51RT
#' - Gagp17Matrix
#' - Envgp41
#' @param distance_metric. Numeric. Currently supported distance metrics include:
#' - mean pairwise distance from the "raw" model (rawMPD)
#' - mean pairwise distance from the "TN93" model (tn93MPD)
#' - nucleotide diversity from the "raw" model (rawPI)
#' - nucleotide diversity from the "TN93" model (tn93PI)
#' - weighted fraction of polymorphic sites (WFPS)
#' - WFPS at third codon positions only (WFPS_codons)
#' @param x_new. Numeric. distance to be used to estimate time since infection.
#'
#' @return predicted, lower bound of 95% credible interval, upper bound of 95% credible interval.
#' @examples
#' TSI_rawMPD = predict_TSI("Gagp17Matrix", "rawMPD", 0.003)
#' TSI_WFPS = predict_TSI("Gagp17Matrix", "WFPS", 0.003)
#'
#' @export
predict_TSI <- function(hiv_region, distance_metric, x_new){

  # load mcmc_samples from bayesian linear regression
  all_mcmc_samples <- readRDS("data/samples_model1.rds")

  # select mcmc sample for this gene hiv_region and distance metric
  model_name <- grep(paste(hiv_region,distance_metric, sep="_"), names(all_mcmc_samples), value = T)

  if (distance_metric == "WFPS_codons"){
    model_name <- grep(paste(hiv_region,"WFPS_codons", sep="_"), names(all_mcmc_samples), value = T)
  } else if (distance_metric == "WFPS"){
    model_name <- grep(paste(hiv_region,"WFPS$", sep="_"), names(all_mcmc_samples), value = T)
  } else {
    model_name <- grep(paste(hiv_region,distance_metric, sep="_"), names(all_mcmc_samples), value = T)
  }

  mod_sim <- all_mcmc_samples[[model_name]]
  mod_csim = as.mcmc(do.call(rbind, mod_sim))

  # Compute posterior predictive distribution
  y_pred_samples <- mod_csim[,"a"] + mod_csim[,"b"] * x_new

  # Compute mean and credible interval (95%)
  y_pred_mean <- mean(y_pred_samples)
  y_pred_lwr <- quantile(y_pred_samples, probs = 0.025)
  y_pred_upr <- quantile(y_pred_samples, probs = 0.975)
  res <- c(y_pred_mean, y_pred_lwr, y_pred_upr)
  names(res) = c("pred",  "lwr", "upr")
  print(res)
  return(res)

}
