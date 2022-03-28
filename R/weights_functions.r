# Functions for model stacking,  

# Read in a bunch of loo files (in parallel)
read_loo_matrix = function(file_names, cores = 48){
  loo_lpd_list = parallel::mclapply(file_names, 
                                    function(.x) {
                                      read_rds(.x)$pointwise[, "elpd_loo"]
                                    }, mc.cores = cores)
  loo_lpd_matrix = do.call(cbind, loo_lpd_list)
  loo_lpd_matrix
}

# More efficient version of loo::stacking_weeights()
# seed that function for doc
stacking_weights_fast <-
  function(lpd_point,
           optim_method = "BFGS",
           optim_control = list()) {
    
    stopifnot(is.matrix(lpd_point))
    N <- nrow(lpd_point)
    K <- ncol(lpd_point)
    if (K < 2) {
      stop("At least two models are required for stacking weights.")
    }
    
    exp_lpd_point <- exp(lpd_point)
    simplex_exp_lpd =  exp_lpd_point[,-K] - exp_lpd_point[, K]
    
    negative_log_score_loo <- function(w) {
      # objective function: log score
      # stopifnot(length(w) == K - 1)
      w_full <- c(w, 1 - sum(w))
      # cat("iter\n")
      -sum(log(exp_lpd_point %*% w_full))
    }
    
    # Numerator for original gradient block
    
    gradient <- function(w) {
      # gradient of the objective function
      stopifnot(length(w) == K - 1)
      w_full <- c(w, 1 - sum(w))
      # Denominator of the original code (inverted for multiplication)
      point_weights = 1 / t(exp_lpd_point  %*% w_full)
      # Gradient:
      -as.numeric(point_weights %*% simplex_exp_lpd)
    }
    
    ui <- rbind(rep(-1, K - 1), diag(K - 1))  # K-1 simplex constraint matrix
    ci <- c(-1, rep(0, K - 1))
    w <- constrOptim(
      theta = rep(1 / K, K - 1),
      f = negative_log_score_loo,
      grad = gradient,
      ui = ui,
      ci = ci,
      method = optim_method,
      control = optim_control
    )$par
    
    wts <- structure(
      c(w, 1 - sum(w)),
      names = paste0("model", 1:K),
      class = c("stacking_weights")
    )
    
    return(wts)
  }