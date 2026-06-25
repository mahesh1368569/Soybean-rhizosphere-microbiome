library(readxl)
library(tidyverse)
library(janitor)
library(car)
library(emmeans)
library(multcomp)
library(multcompView)
library(broom)
library(writexl)
library(lme4)
library(lmerTest)
library(performance)

dat <- read_excel("/Users/durgasmacmini/Documents/R-Git/Soybean-microbiome/Data files/Soybean_env.xlsx", sheet = 1) %>%
  clean_names()

names(dat)

dat <- dat %>%
  rename(
    rep = rep,
    n = nitrogen,
    s = sulphur
  ) %>%
  mutate(
    rep = factor(rep),
    n = factor(n, levels = c(0, 100, 150, 200)),
    s = factor(s, levels = c(0, 15, 25, 35))
  )

dat %>%
  count(n, s)

dat %>%
  count(rep, n, s)

str(dat)

id_vars <- c("potnumber", "sample_id", "exp_id", "rep", "n", "s")

response_vars <- setdiff(names(dat), id_vars)

response_vars

options(contrasts = c("contr.sum", "contr.poly"))

fit_one_lm <- function(response) {
  
  df <- dat %>%
    dplyr::select(dplyr::all_of(c(response, "rep", "n", "s"))) %>%
    tidyr::drop_na()
  
  df <- df %>%
    dplyr::mutate(
      rep = factor(rep),
      n = factor(n),
      s = factor(s)
    )
  
  formula_lm <- as.formula(paste(response, "~ rep + n * s"))
  
  model <- lm(formula_lm, data = df)
  
  anova_table <- car::Anova(model, type = 2) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term") %>%
    dplyr::filter(term %in% c("rep", "n", "s", "n:s")) %>%
    dplyr::mutate(response = response) %>%
    dplyr::select(response, term, dplyr::everything())
  
  return(anova_table)
}

all_anova_results <- purrr::map_dfr(response_vars, fit_one_lm)

all_anova_results

write_csv(all_anova_results, "Output_results/all_anova_results.csv")

### Means
raw_ns_summary <- dat %>%
  pivot_longer(
    cols = all_of(response_vars),
    names_to = "response",
    values_to = "value"
  ) %>%
  group_by(response, n, s) %>%
  summarise(
    n_obs = sum(!is.na(value)),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    se = sd / sqrt(n_obs),
    .groups = "drop"
  )

raw_ns_summary

write_csv(raw_ns_summary, "Output_results/raw_N_by_S_mean_SD_SE.csv")

## pairwise mean differencss

get_pairwise_ns <- function(response) {
  
  df <- dat %>%
    dplyr::select(dplyr::all_of(c(response, "rep", "n", "s"))) %>%
    tidyr::drop_na()
  
  df <- df %>%
    mutate(
      rep = factor(rep),
      n = factor(n),
      s = factor(s)
    )
  
  model <- lm(as.formula(paste(response, "~ rep + n * s")), data = df)
  
  emm <- emmeans(model, ~ n * s)
  
  pairs_out <- pairs(emm, adjust = "tukey") %>%
    as.data.frame() %>%
    mutate(response = response) %>%
    dplyr::select(response, contrast, estimate, SE, df, t.ratio, p.value)
  
  return(pairs_out)
}

pairwise_ns_differences <- purrr::map_dfr(response_vars, get_pairwise_ns)

pairwise_ns_differences

write_csv(pairwise_ns_differences, "Output_results/pairwise_N_by_S_mean_differences.csv")

dir.create("leaf_plots", showWarnings = FALSE)

for (resp in response_vars) {
  
  plot_data <- dat %>%
    dplyr::select(all_of(c(resp, "n", "s"))) %>%
    drop_na() %>%
    group_by(n, s) %>%
    summarise(
      mean = mean(.data[[resp]], na.rm = TRUE),
      sd   = sd(.data[[resp]], na.rm = TRUE),
      n_obs = sum(!is.na(.data[[resp]])),
      se   = sd / sqrt(n_obs),
      .groups = "drop"
    )
  
  p <- ggplot(plot_data, aes(x = n, y = mean, group = s, color = s)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.15) +
    labs(
      title = paste("Effect of Nitrogen and Sulphur on", resp),
      x = "Nitrogen",
      y = resp,
      color = "Sulphur"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11)
    )
  
  ggsave(
    filename = paste0("leaf_plots/", resp, "_interaction_plot.png"),
    plot = p,
    width = 7,
    height = 5,
    dpi = 300
  )
}

dir.create("violin_NS_raw", showWarnings = FALSE)
dir.create("violin_N_main", showWarnings = FALSE)
dir.create("violin_S_main", showWarnings = FALSE)

plot_violin_ns <- function(response) {
  
  df <- dat %>%
    dplyr::select(all_of(c(response, "n", "s"))) %>%
    tidyr::drop_na() %>%
    mutate(ns = interaction(n, s, sep = "_"))
  
  p <- ggplot(df, aes(x = ns, y = .data[[response]])) +
    geom_violin(trim = FALSE, fill = "grey85", color = "black") +
    geom_jitter(width = 0.12, size = 2, alpha = 0.8) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 3) +
    labs(
      title = paste("Raw distribution of", response, "across N × S treatments"),
      x = "Nitrogen × Sulphur",
      y = response
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  return(p)
}

for (resp in response_vars) {
  p <- plot_violin_ns(resp)
  
  ggsave(
    filename = paste0("violin_NS_raw/", resp, "_violin_NS_raw.png"),
    plot = p,
    width = 8,
    height = 5,
    dpi = 300
  )
}

## Main effects violin plots for N

plot_violin_n <- function(response) {
  
  df <- dat %>%
    dplyr::select(all_of(c(response, "n"))) %>%
    tidyr::drop_na()
  
  p <- ggplot(df, aes(x = n, y = .data[[response]])) +
    geom_violin(trim = FALSE, fill = "grey85", color = "black") +
    geom_jitter(width = 0.12, size = 2, alpha = 0.8) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 3) +
    stat_summary(
      fun.data = mean_se,
      geom = "errorbar",
      width = 0.15
    ) +
    labs(
      title = paste("Main effect of Nitrogen on", response),
      x = "Nitrogen",
      y = response
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  return(p)
}

for (resp in response_vars) {
  p <- plot_violin_n(resp)
  
  ggsave(
    filename = paste0("violin_N_main/", resp, "_violin_N_main.png"),
    plot = p,
    width = 6,
    height = 5,
    dpi = 300
  )
}

## Main effect violin plots for Sulphur

plot_violin_s <- function(response) {
  
  df <- dat %>%
    dplyr::select(all_of(c(response, "s"))) %>%
    tidyr::drop_na()
  
  p <- ggplot(df, aes(x = s, y = .data[[response]])) +
    geom_violin(trim = FALSE, fill = "grey85", color = "black") +
    geom_jitter(width = 0.12, size = 2, alpha = 0.8) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 3) +
    stat_summary(
      fun.data = mean_se,
      geom = "errorbar",
      width = 0.15
    ) +
    labs(
      title = paste("Main effect of Sulphur on", response),
      x = "Sulphur",
      y = response
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  return(p)
}

for (resp in response_vars) {
  p <- plot_violin_s(resp)
  
  ggsave(
    filename = paste0("violin_S_main/", resp, "_violin_S_main.png"),
    plot = p,
    width = 6,
    height = 5,
    dpi = 300
  )
}
