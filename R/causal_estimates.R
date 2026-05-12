# Functions to take Brms models and derive Causal Estimates from them ####
library(posterior)
library(ggdist)
library(marginaleffects)

# Function to Extract Draws and calculate the ATE
extract_effect <- function(model, data) {
  # Transform
  data_for_pred <- data %>%
    mutate(treatment = as.numeric(treatment == "intervention"))

  # Extract Mean Predictions
  pred_draws <- avg_predictions(
    model = model,
    by = "treatment",
    type = "response",
    newdata = data_for_pred
  ) %>%
    posterior_draws()

  # Get the predictions for T 1 and T 0
  ate_draws <- pred_draws %>%
    pivot_wider(
      id_cols = drawid,
      names_from = treatment,
      values_from = draw,
      names_prefix = "mu_"
    ) %>%
    mutate(ATE = mu_1 - mu_0)

  # Combine with draws
  combine_return <- list(
    prediction_draws = pred_draws,
    ate_draws        = ate_draws
  )

  # Return the ATE
  return(combine_return)
}

### Theme Function ###
theme_bayes <- function(base_size = 12) {
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.15)),
      plot.subtitle = element_text(colour = "grey30", size = rel(0.95)),
      plot.caption = element_text(
        colour = "grey40", size = rel(0.8),
        hjust = 0
      ),
      axis.title = element_text(face = "plain"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(colour = "grey92"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

#### Function to plot the ATE posterior Distributions ####
plot_ate_posterior <- function(
    ate_draws,
    outcome_label = "outcome",
    effect_color = "red4"
    ) {
  
  # Color Control 
  pal_bayes <- c(
    control = "#8E9AAF", treatment = "#264653",
    effect = effect_color, null = "#E76F51"
  )
  
  # ATE Plot 
  ate_posterior <- ggplot(data = ate_draws, aes(x = ATE)) +
    stat_halfeye(
      .width         = c(0.66, 0.95),
      point_interval = "median_qi",
      fill  = pal_bayes["effect"],
      alpha = 0.85
    )+
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30") +
    labs(
      x = "ATE",
      y = NULL,
      title    = "Posterior of the average treatment effect",
      subtitle = "Median, 66% and 95% credible intervals",
      caption  = "Dashed line = no effect."
    ) +
    theme_bayes()
  
  # Return the plot 
  return(ate_posterior)
}

#### Function to plot the P of effect exceeding a threshold 
plot_threshold_curve <- function(
  ate_draws, 
  favours_neg = TRUE,
  outcome_label = "outcome"
  ) {
  
  # Find Treashholds
  thr <- seq(min(ate_draws$ATE), max(ate_draws$ATE), length.out = 250)
  
  # Switch between effects of interests
  prob_fun <- if (favours_neg) {
    function(c) mean(ate_draws$ATE < c)
  } else {
    function(c) mean(ate_draws$ATE > c)
  }
  
  # Combine into dataframe
  df <- data.frame(c = thr, prob = vapply(thr, prob_fun, numeric(1)))
  
  # Color Control 
  pal_bayes <- c(
    control = "#8E9AAF", treatment = "#264653",
    effect = effect_color, null = "#E76F51"
  )
  
  # Plot the Decision Curve
  decision_curve <- ggplot(data = df ,aes(x = c, y = prob)) +
    geom_area(fill = pal_bayes["effect"], alpha = 0.25) +
    geom_line(colour = pal_bayes["effect"], linewidth = 1) +
    geom_hline(
      yintercept = c(0.5, 0.95),
      linetype = "dotted",
      colour = "grey40"
      ) +
    annotate(
      "text", x = max(thr), y = 0.95, label = "95% certainty",
       vjust = -0.4, hjust = 1, size = 3, colour = "grey30"
      ) +
    annotate(
      "text", x = max(thr), y = 0.50, label = "50% certainty",
      vjust = -0.4, hjust = 1, size = 3, colour = "grey30"
      ) +
    labs(
      x = paste0("Threshold c (", outcome_label, " units)"),
      y = "",
      title    = "Probability the effect exceeds a threshold",
      subtitle = "Replaces the binary question of significance with a continuous one") +
    theme_bayes()
  
  # Return the plot 
  return(decision_curve)
}

#### Function to plot Counterfactual marginal means on the outcome scale ####
plot_counterfactual_means <- function(
    pred_draws,
    outcome_label = "outcome",
    effect_color = "red4"
    ){
  
  # Color Control 
  pal_bayes <- c(
    control = "#8E9AAF", treatment = "#264653",
    effect = effect_color, null = "#E76F51"
  )
  
  # Plot Counterfactual Means
  pc_means <- ggplot(data = pred_draws, aes(
    x = factor(treatment, labels = c("Untreated", "Treated")),
    y = draw,
    fill = factor(treatment))) +
    stat_halfeye(
      .width         = c(0.66, 0.95),
      point_interval = "median_qi") +
    scale_fill_manual(
      values = unname(pal_bayes[c("control", "treatment")])) +
    labs(
      x = "Comparison Between Treated and Control" , 
      y = outcome_label,
      title    = "Counterfactual marginal means",
      subtitle = "Average outcome under each treatment arm") +
    theme_bayes()
  
  # Return the plot 
  return(pc_means)
}

#### Plot the Posterior Distribution ####

  
  
  
  
  
  











