# Install/Load required packages -----------------------------------------------
#' Check packages, install if not available, then load.
#'
#' @param A vector of packages
#'
#' @return Attach packages in the global environment or R session
#'
#' @example
#' install_load_packages(c("dplyr", "emmeans"))
install_load_packages <- function(packages) {
  # Check and install packages if not already installed.
  # Requires a vector of package names.

  for (package in packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
      install.packages(package)
    }
  }

  # Load all the packages
  sapply(packages, function(pkg) {
    suppressPackageStartupMessages({
      library(pkg, character.only = TRUE, warn.conflicts = FALSE)
    })
  })
}


# Plot data by block -----------------------------------------------------------
#' Visualize data by block
#'
#' @param dt A dataframe containing the data to be plotted.
#' Expected columns:
#' - `measure`: Numeric. The measure or response under investigation
#' - `treatment`: Factor/Character. Treatment groups
#' - `block_var`: Factor/Character. Blocking variable
#' @param y_label A character string specifying the label for the y-axis.
#' @param cbb_palette A character vector of colors for the boxplots.
#' @return A ggplot object representing the boxplot.
#'
#' @examples
# create_boxplot(df, "Growth rate")
create_boxplot <- function(df, y_label) {
  ggplot(data = df, aes(y = measure, x = treatment, fill = treatment)) +
    geom_boxplot(alpha = 0.4) +
    geom_jitter() +
    facet_grid(. ~ block_var) +
    scale_colour_manual(values = cbb_palette, aesthetics = "fill") +
    labs(y = y_label, x = "Treatment") +
    theme_bw() +
    theme(
      axis.title = element_text(size = 16),
      axis.text.x = element_text(angle = 340, hjust = 0),
      legend.title = element_blank(),
      legend.text = element_text(size = 12)
    )
}


# Generate summary data for a variable of interest ----------------------------
#' Title: Five-numbers summaries by treatment and blocking, in wide format
#'
#' @details For use in random effects meta-analysis.
#'
#' @param `df` dataframe
#' Expected columns:
#' - `measure`: Numeric. The measure or response under investigation
#' - `treatment`: Factor/Character. Treatment variable (grouping variable)
#' - `block_var`: Factor/Character. Blocking variable (grouping variable)
#' @return wide-formated dataframe of a reference and a treatment summaries:
#' Columns:
#' - `n`: Integer. Number of observations
#' - `mean`: Numeric. Mean of observations
#' - `median`: Numeric. Median of observations
#' - `sd`: Numeric. Standard deviation of observations
#'
#' @examples
#' summariseMB_wide(df)
summariseMB_wide <- function(df) {
  df %>%
    group_by(block_var, treatment) %>%
    dplyr::summarise(
      mean = mean(measure, na.rm = TRUE),
      median = median(measure, na.rm = TRUE),
      SD = sd(measure, na.rm = TRUE),
      no_Dp = n(),
      .groups = "drop"
    ) %>%
    arrange(treatment) %>%
    pivot_wider(
      names_from = treatment,
      values_from = c(no_Dp, mean, median, SD)
    ) %>%
    set_names(c(
      "block_var",
      "no_DP_V", "no_DP_T",
      "mean_V", "mean_T",
      "median_V", "median_T",
      "SD_V", "SD_T"
    )) %>%
    relocate("block_var", ends_with("_V"), ends_with("_T"))
}


# Linear model --------------------------------------------------------------
#' Fit linear model and calculate pairwise contrast between treatments.
#'
#' @param `dt_split_trt` List of dataframes
#' Expected format of dataframes:
#' - `measure`: Numeric. The measure or response under investigation
#' - `block_var`: Factor/Character. Blocking variable (grouping variable)
#' - `treatment`: Factor/Character. Two levels treatment variable
#'    - `ref_trt`: Reference
#'    - `trt_interest`: Treatment of interest
#' @param `trt_interest` Treatments of interest
#' @param `ref_trt` Reference treatment
#'
#' @return List of fitted model, and pairwise treatment contrast objects
#'
#' @examples
#' fit_lm_contrast_trt(dt_split_trt, trt_interest, ref_trt)
fit_lm_contrast_trt <- function(dt_split_trt, trt_interest, ref_trt) {
  # Use map2 to iterate over dt_split_trt and trt_interest
  results <- map2(dt_split_trt, trt_interest, ~ {
    trt_str <- .y
    veh_str <- ref_trt

    # Create the linear model
    m <- lm(measure ~ treatment + block_var + block_var:treatment,
      data = .x %>% droplevels(),
      na.action = "na.exclude",
      contrasts = list(treatment = "contr.sum", block_var = "contr.sum")
    )

    # Save the original data in the model object
    m$call$data <- .x %>% droplevels()

    # Calculate the estimated marginal means
    m_emm <- emmeans(m, tukey ~ treatment)
    
    m_est <- tibble(
        estimates <- summary(m_emm)$emmean,
        standard_errors <- summary(m_emm)$SE,
        lower_conf <- summary(m_emm)$lower.CL,
        upper_conf <- summary(m_emm)$upper.CL
    )
    
    m_contrasts <- contrast(m_emm, method = "pairwise")

    # Return a list of results for each iteration
    list(model = m, emmeans_estimate = m_est, emmeans_contrast = m_contrasts)
  })

  return(results)
}


#' Print tables of fitted linear model pairwise contrast between treatments.
#' 
#' @param lm_results List of fitted LM `model` and pairwise `emmeans_estimates`
#' 
#' @return Html tables of model and emmeans results
#' 
#' @example 
#' process_results(lm_results)
process_results <- function(lm_results) {
  lm_results %>% imap(~{
    emmeans_estimate_df <- .x %>% 
      pluck("emmeans_estimate") 
    
    emmeans_contrast_df <- .x %>% 
      pluck("emmeans_contrast") 
    
    # model_coef_df <- .x %>% pluck("model")
    
    # if(inherits(model_coef_df, "lm")){
    #     model_coef_df <- model_coef_df %>% summary() %>% broom::tidy()
    # }
    # if(inherits(model_coef_df, "lme")){
    #     model_coef_df <- model_coef_df %>% broom.mixed::tidy()
    # }

    # cat(str_glue("Model Coefficients: {.y}"))
    # model_coef_df %>%
    #     kable(format = "pipe") %>%
    #     kable_styling(bootstrap_options = c("striped", "hover")) %>%
    #     print()

    cat(str_glue("EM Means Estimate:{.y}"))
    emmeans_estimate_df %>%
        kable(format = "pipe") %>%
        kable_styling(bootstrap_options = c("striped", "hover")) %>%
        print()

    cat(str_glue("EM Means Contrast:{.y}"))
    emmeans_contrast_df %>%
        kable(format = "pipe") %>%
        kable_styling(bootstrap_options = c("striped", "hover")) %>%
        print()
  })
  invisible()
}


# Linear mixed model -----------------------------------------------------------
#' Fit LMM and perform pairwise contrasts on the fitted model.
#'
#' @param dt_split_trt List of dataframes
#' Expected format of dataframes:
#' - `measure`: Numeric. The measure or response under investigation
#' - `block_var`: Factor/Character. Blocking variable (grouping variable)
#' - `treatment`: Factor/Character. Two levels treatment variable
#'    - `ref_trt`: Reference
#'    - `trt_interest`: Treatment of interest
#' @return List of fitted and pairwise treatment contrasts objects.
#'
#' @examples
#' fit_lmm_contrast_trt(dt_split_trt)
fit_lmm_contrast_trt <- function(dt_split_trt) {
  # Use map to iterate over dt_split_trt
  results <- map(dt_split_trt, ~ {
    # Create the mixed-effects model
    m <- lme(measure ~ treatment, 
      random = ~ 1 | block_var / treatment,
      data = .x %>% droplevels(),
      method = "REML", na.action = "na.omit"
    )

    # Save the original data in the model object
    m$call$data <- .x %>% droplevels()

    # Calculate the estimated marginal means
    ## m_estimate <- pairs(emmeans(m, ~treatment, adjust = "none"))
    m_emm <- emmeans(m, tukey ~ treatment)
    
    m_est <- tibble(
        estimates <- summary(m_emm)$emmean,
        standard_errors <- summary(m_emm)$SE,
        lower_conf <- summary(m_emm)$lower.CL,
        upper_conf <- summary(m_emm)$upper.CL
    )
    
    m_contrasts <- contrast(m_emm, method = "pairwise")
    
    # Return a list of results for each iteration
    list(model = m, emmeans_estimate = m_est, emmeans_contrast = m_contrasts)
  })

  return(results)
}


# Generate output of regression model parameters -------------------------------
#' Output of fitted LM, LMM, and their pairwise treatment contrasts
#'
#' @param `results_lm` A list of fitted linear models, and their pairwise
#'                   treatment contrasts.
#' @param `results_lmm` A list of fitted linear mixed models, and their pairwise
#'                    treatment contrast.
#'
#' @return A list of data frames:
#' Columns:
#' - `model params`: Character. Name of model parameters
#' - `value`: Numeric. Estimates of model parameters
#'
#' @examples
#' create_output(lm_results, lmm_results)
# pmap(
#     list(lm = results_lm, lmm = results_lmm, trt_name = trt_interest),
#     function(lm, lmm, trt_name) {
#         # Get estimated model parameters
#         lm_fit <- summary(lm$model) %>%
#             tidy() %>%
#             add_column(model = "lm", .before = "term")
#         lmm_fit <- lmm$model %>%
#             broom::tidy() %>%
#             add_column(model = "lmm", .before = "effect")
#         shared_cols <- intersect(names(lm_fit), names(lmm_fit))
#         # Combine estimated model parameters
#         bind_rows(
#             lm_fit %>% select(all_of(shared_cols)), 
#             lmm_fit %>% select(all_of(shared_cols))
#         ) %>%
#             mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
#             filter(str_detect(term, "treatment1$")) %>%
#             # Create and style table of estimated model parameters
#             kable(
#                 .,
#                 caption = str_glue("Treatment: {trt_name}"),
#                 booktabs = TRUE,
#             ) %>%
#             kable_styling(bootstrap_options = c("striped", "hover"))
#     }
# )
#

# Generate output of meta-analytic random effects model parameters ----------
#' @param `r_ma_model_lst`  A list of fitted meta-analytic mixed-effects models.
#' @return A list of data frames:
#' Columns:
#' - `model params`: Character. Name of model parameters
#' - `value`: Numeric. Estimates of model parameters
#' @examples
#' create_ma_output(r_ma_model_results)
create_ma_output <- function(r_ma_model_lst) {
    n <- length(r_ma_model_lst)
    # Iterate over the index
    output <- map(seq_len(n), ~ {
        index <- .x
        
        
        # Extracting meta-analysis results if available
        ma_info <- vector()
        if (index <= length(r_ma_model_lst)) {
            ma_result <- r_ma_model_lst[[index]]
            ma_info <- c(
                estimate = round(ma_result$b[1], 4),
                se = round(ma_result$se[1], 4),
                pvalue = signif(ma_result$pval, digits = 4),
                I2 = round(ma_result$I2, 2),
                testOfHet = signif(ma_result$QEp, digits = 4)
            )
        }
        ma_info %>% enframe()
    }) %>%
        set_names(names(r_ma_model_lst))
    return(output)
}


# Model diagnostic plots to check the validity of models -----------------------

#' Generate Pearson residual plot w.r.t. blocking and treatment terms of models.
#'
#' @param `df` A dataframe that has columns of model terms.
#' Expected format of dataframe:
#' - `measure`: Numeric. The measure or response under investigation
#' - `block_var`: Factor/Character. Blocking variable (grouping variable)
#' - `treatment`: Factor/Character. Two levels treatment variable
#'    - `ref_trt`: Reference
#'    - `trt_interest`: Treatment of interest
#' @param `model` A fitted LM or LMM regression model
#' @param `groupingVariable` A blocking variable; in the data frame, `df`
#' @param `approach` The method used to calculate the model residuals.
#'
#' @return A Pearson residual plot.
#'
#' @examples
#' ResVpred(dt, model = model_obj, groupingVariable = "block_var")
ResVpred <- function(df,
                     model,
                     groupingVariable,
                     approach = "pearson",
                     treatment_string) {
  pred <- predict(model)
  res <- resid(model, type = approach)
  data_all <- data.frame(df, res, pred)
  plot <- ggplot(
    data = data_all,
    aes(
      y = {{ res }}, x = {{ pred }},
      colour = !!sym(groupingVariable)
    )
  ) +
    geom_point() +
    scale_colour_manual(values = cbb_palette, aesthetics = "colour") +
    labs(y = "Residual", x = "Predicted") +
    theme_bw() +
    ggtitle(treatment_string)

  print(plot)
}


#' Diagnostic plot to check normality of linear mixed model residuals.
#'
#' @param `df` A dataframe that has columns of model terms.
#' Expected format of dataframe:
#' - `measure`: Numeric. The measure or response under investigation
#' - `block_var`: Factor/Character. Blocking variable (grouping variable)
#' - `treatment`: Factor/Character. Two levels treatment variable
#'    - `ref_trt`: Reference
#'    - `trt_interest`: Treatment of interest
#' @param `model` A fitted LMM regression model.
#' @param `approach` The method used to calculate the model residuals.
#'
#' @return A Q-Q plot.
#'
#' @examples
#' ResidualQQplot(model = model_obj, df = dt)
ResidualQQplot <- function(df, model, approach = "pearson", title) {
  res <- resid(model)
  data_all <- data.frame(df, res)
  qqnorm(data_all$res, main = title)
  qqline(data_all$res)
}
