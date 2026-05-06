#### Functions to help in the Descriptive Statistics #####
library(ggstatsplot)
library(patchwork)

#### Function to compute delta measures Pre Treatment - Post Treatment and return a dataframe #####
compule_delta_anthropometrics <- function(data_clean) {
  # Delta Anthropometrics
  data_delta <- data_clean %>%
    mutate(
      delta_height = height_pre - height_post,
      delta_weight = weight_pre - weight_post,
      delta_bmi    = bmi_pre - bmi_post,
      delta_waist  = waist_pre - waist_post,
      delta_hip    = hip_pre - hip_post
    )
  return(data_delta)
}

#### Function to compute delta measures Pre Treatment Exercise Numerical  - Post Treatment Exercise Numerical
compule_delta_exercise_metrics <- function(data_clean) {
  # Delta Exercimetrics
  data_delta <- data_clean %>%
    mutate(
      delta_exercise_type_1_freq = exercise_type_1_freq_pre - exercise_type_1_freq_post,
      delta_exercise_type_2_freq = exercise_type_2_freq_pre - exercise_type_2_freq_post
    )
  return(data_delta)
}

#### Function to compute delta measures Pre Treatment Diet Quantities - Post Treatment Diet Quantities
compute_delta_dietometrics <- function(data_clean) {
  # Delta Deitometrics
  data_delta <- data_clean %>%
    mutate(
      delta_carbs = carb_pre - carb_post,
      delta_energy = energy_pre - energy_post,
      delta_diet_score = diet_score_pre - diet_score_post
    )
  return(data_delta)
}

### Marginal distributions of post-outcomes facet by treatment arm ####
marginal_distribution <- function(data_clean, x, title) {
  # Histogram
  hist <- grouped_gghistostats(
    data = data_clean,
    x = {{ x }},
    grouping.var = treatment,
    results.subtitle = FALSE,
    centrality.plotting = FALSE,
    normal.curve = TRUE,
    ggtheme = ggthemes::theme_tufte()
  )
  # Boxplot
  boxplot <- ggbetweenstats(
    data = data_clean,
    x = treatment,
    y = {{ x }},
    xlab = "Condition",
    ylab = rlang::as_label(rlang::enquo(x))
  )
  # Combine plots
  combined <- (hist / boxplot) +
    plot_annotation(title = title)
  return(combined)
}

### Pairwise correlations of baselines ####
pairwise_corrleation_map <- function(data_clean, vars) {
  # Subset the data to include only the specifyed vars and plot them
  cor_map_full <- ggcorrmat(
    data = data_clean,
    cor.vars = vars,
    colors = c("#B2182B", "white", "#4D4D4D"),
    type = "nonparametric"
  )
  return(cor_map_full)
}

### Pre-vs-post scatter for each outcome facet by treatment ###
scatter_pre_post <- function(data_clean, x, y) {
  scatter_fig <-
    grouped_ggscatterstats(
      data = data_clean,
      x = {{ x }},
      y = {{ y }},
      grouping.var = treatment,
      annotation.args = list(
        title = paste0(
          "Relationship between ",
          rlang::as_name(rlang::enquo(x)),
          " and ",
          rlang::as_name(rlang::enquo(y)),
          " by treatment"
        )
      )
    )
  return(scatter_fig)
}
