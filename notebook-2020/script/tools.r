
#########################################################################
## Tools to analyze results from the cluster

select_lambda <- function(parm, quant, n_select) {
  p <- ncol(parm)
  n <- nrow(parm)
  quantile_parm <- rep(NA, p)
  for (i in 1:p) quantile_parm[i] <- sort(parm[, i])[quant * n]
  selected <- sort(quantile_parm, decreasing = T)[1:n_select]
  covariate <- rep(NA, n_select)
  for (i in 1:n_select) covariate[i] <- which(quantile_parm %in% selected[i])
  covariate
  # which(quantile_parm %in% selected)
}

construct_plot_data <- function(parm, nIter, nChains, names) {
  iteration <- rep(1:nIter, nChains)
  chain <- rep(1:nChains, each = nIter)
  posterior.sample <- data.frame(parm, iteration, as.factor(chain))
  # names(posterior.sample) <- c(paste0("log_lambda[", index, "]"), "iteration", "chain")
  names(posterior.sample) <- names
  posterior.sample <- posterior.sample %>%
    gather(key = parameter, value = value, -chain, -iteration)
}

trace_plot <- function(posterior.sample) {
  trace.plot <- ggplot(data = posterior.sample,
                       aes(x = iteration, y = value, color = chain)) +
    geom_line() + theme_bw() + facet_wrap(~ parameter)
  print(trace.plot)
}

density_hist <- function(posterior.sample, bins = 30) {
  density.plot <- ggplot(data = posterior.sample,
                         aes(x = value, color = chain, fill = chain)) +
    geom_histogram(alpha = 0.25, position = "identity", bins = bins) + 
    theme_bw() + facet_wrap(~ parameter) +
    theme(text = element_text(size = 15))
  print(density.plot)
}

quant_select_plot <- function(parm, quant, threshold = 3.5) {
  index <- 1:ncol(parm)
  parm_quant <- apply(parm, 2, quantile, quant)
  ggplot(data = data.frame(index = index, parm_quant = parm_quant),
         aes(x = index, y = parm_quant)) + geom_point(size = 0.25) +
    geom_text(aes(label = ifelse(parm_quant > threshold, index, '')),
              hjust = 0, vjust = 0) +
    theme_bw() + theme(text = element_text(size = 18))
}

quant_select_plot2 <- function(parm1, parm2, quant, threshold = 3.5, alpha = 0.05,
                               x = 0.95, y = 0.8, index_offset = 0) {
  index <- c(1:ncol(parm1), 1:ncol(parm2)) + index_offset
  parm_quant <- c(apply(parm1, 2, quantile, quant),
                  apply(parm2, 2, quantile, quant))
  method <- c(rep("(full) HMC", ncol(parm1)), rep("HMC + Laplace", ncol(parm2)))

  ggplot(data = data.frame(index = index, parm_quant = parm_quant, method = method),
         aes(x = index, y = parm_quant, color = method)) + 
    geom_text(aes(label = ifelse(parm_quant > threshold, index, '')),
              hjust = 0, vjust = 0) +
    geom_point(size = 0.25, alpha = alpha) +
    theme_bw() +
    theme(
      legend.position = c(x, y),
      legend.justification = c("right", "top"),
      legend.box.just = "right",
      legend.margin = margin(6, 6, 6, 6),
      text = element_text(size = 15)
    ) + xlab("covariate index") + ylab("90th quantile")
}

summary_table <- function(log_lambda, tau, caux, index) {
  standard_post <- data.frame(log_lambda, tau, caux)
  names(standard_post) <- c(paste0("log_lambda[", index, "]"), 
                            "tau", "caux")
  standard_post <- as_draws(standard_post)
  draw_summary <- summarise_draws(standard_post)
  time <- time <- sum(colSums(get_elapsed_time(stanfit))) + 
    sum(colSums(get_elapsed_time(stanfit2)))
  draw_summary$eff_bulk <- draw_summary$ess_bulk / time
  draw_summary$eff_tail <- draw_summary$ess_tail / time
  draw_summary
}


sample_comparison_plot <- function(plot_data, x = 1, y = 0.79) {
  ggplot(data = plot_data) +
    geom_histogram(aes(x = value, fill = method), alpha = 0.5, color = "black",
                   bins = 30, position = "identity") + theme_bw() +
    facet_wrap(~key, scale = "free", nrow = 1, labeller = "label_parsed") +
    theme(
      legend.title = element_blank(),
      legend.position = c(x, y),
      legend.justification = c("right", "top"),
      legend.box.just = "right",
      legend.margin = margin(6, 6, 6, 6),
      text = element_text(size = 25)
    )
}

eff_comparison_plot <- function(plot_data, x = 0.95, y = 0.98) {
  ggplot(data = plot_data,
         aes(x = parameter, y = eff, fill = method)) +
    geom_bar(stat = "identity", width = 0.3, alpha = 0.8, position = "dodge") + 
    # facet_wrap(~ parameter, scale = "free", nrow = 1) +
    theme_bw() + theme(text = element_text(size = 10)) + coord_flip() +
    ylab("ESS / s") + xlab(" ") +
    theme(
      legend.position = c(x, y),
      legend.justification = c("right", "top"),
      legend.box.just = "right",
      legend.margin = margin(6, 6, 6, 6)
    ) 
}
