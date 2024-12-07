#' Plate Analysis Package
#' 
#' @description Functions for analyzing high-throughput screening plate data with NRFE metric
#' @name plateQC
#' @import ggplot2
#' @import parallel
#' @import drc
#' @import minpack.lm
#' @import caTools
#' @importFrom dplyr filter group_by mutate arrange left_join %>% n_distinct
#' @importFrom gridExtra arrangeGrob
#' @importFrom stats IQR mad median na.omit quantile var coef fitted nls predict
#' @importFrom utils head tail
NULL

utils::globalVariables(c(
  "CONC", "Column", "DRUG_NAME", "Row", "WELL", 
  "inhibition_percent", "logconc", "text2", "x2", "y2"
))

#' Remove outliers using IQR method
#' 
#' @description Removes outliers from a numeric vector using the Interquartile Range method
#' @param x Numeric vector
#' @return Numeric vector with outliers replaced by NA
#' @export
remove_outliers <- function(x) {
  qq <- unname(quantile(x, probs = c(.25, .75), na.rm = TRUE))
  outlier_detector <- 1.5 * IQR(x, na.rm = TRUE)
  x[x < (qq[1] - outlier_detector) | x > (qq[2] + outlier_detector)] <- NA
  x
}

#' Calculate standard error
#' @param x Numeric vector
#' @return Standard error value
#' @export
calc_stderr <- compiler::cmpfun(function(x) {
  sqrt(var(x, na.rm = TRUE) / length(na.omit(x)))
})

#' Calculate population standard deviation
#' @param x Numeric vector
#' @return Population standard deviation
#' @export
calc_pop_sd <- compiler::cmpfun(function(x) {
  sqrt(var(x) * (length(x) - 1) / length(x))
})

#' Calculate SSMD (Strictly Standardized Mean Difference)
#' @param x,y Numeric vectors to compare
#' @return SSMD value
#' @export
calc_ssmd <- compiler::cmpfun(function(x, y) {
  round((mean(x) - mean(y)) / sqrt(var(x) + var(y)))
})

#' Calculate Z-factor
#' @param x,y Numeric vectors to compare
#' @return Z-factor value
#' @export
calc_zfactor <- compiler::cmpfun(function(x, y) {
  round((1 - (3 * (calc_pop_sd(x) + calc_pop_sd(y))) / (abs(mean(x) - mean(y)))), 2)
})

#' Calculate robust Z-factor
#' @param x,y Numeric vectors to compare
#' @return Robust Z-factor value
#' @export
calc_robust_zfactor <- compiler::cmpfun(function(x, y) {
  round((1 - (3 * (mad(x) + mad(y))) / (abs(median(x) - median(y)))), 2)
})

#' Create plate row labels
#' @param n Number of rows
#' @return Character vector of row labels
#' @export
create_plate_labels <- function(n) {
  single <- LETTERS
  if (n <= 26) return(single[1:n])
  double <- paste0("A", LETTERS)
  all_labels <- c(single, double)
  return(all_labels[1:n])
}

#' Process plate data
#' @param plate_data Data frame containing plate data
#' @return List containing processed data and quality metrics
#' @export
process_plate <- function(plate_data) {
  # Input validation
  if (!is.data.frame(plate_data) || nrow(plate_data) == 0) {
    stop("Invalid plate data provided")
  }
  
  # Get plate metadata
  plate_id <- unique(plate_data$BARCODE)[1]
  
  # Check if inhibition_percent already exists
  if ("inhibition_percent" %in% names(plate_data)) {
    # Process plate data without recalculating inhibition
    plate_data <- plate_data %>% 
      mutate(
        Row = gsub("([A-Z]+).*", "\\1", WELL),
        Column = as.numeric(gsub("[A-Z]+", "", WELL))
      )
    
    # Create complete grid
    complete_wells <- expand.grid(
      Row = unique(plate_data$Row),
      Column = 1:max(plate_data$Column, na.rm = TRUE)
    ) %>%
      mutate(WELL = paste0(Row, Column))
    
    # Final processing
    processed_data <- complete_wells %>%
      left_join(plate_data, by = c("Row", "Column", "WELL")) %>%
      arrange(Row, Column) %>%
      mutate(
        Column = as.numeric(as.character(Column)),
        Row = factor(Row, levels = unique(Row))
      )
    
    # Set all quality metrics to NA
    quality_metrics <- list(
      zfactor = NA,
      ssmd = NA,
      robust_z_prime = NA,
      signal_vs_bg = NA
    )
    
    return(list(
      processed_data = processed_data,
      quality_metrics = quality_metrics
    ))
  }
  
  # If inhibition_percent doesn't exist, continue with original processing
  # Validate controls
  controls <- list(
    POS = plate_data$INTENSITY[plate_data$DRUG_NAME == "POS_CTRL"],
    NEG = plate_data$INTENSITY[plate_data$DRUG_NAME == "NEG_CTRL"]
  )
  
  # Check controls presence and validity
  for (ctrl_type in c("POS", "NEG")) {
    if (length(controls[[ctrl_type]]) == 0) {
      warning(sprintf("Missing %s_CTRL in plate %s", ctrl_type, plate_id))
      return(NULL)
    }
    
    if (all(is.na(controls[[ctrl_type]]))) {
      warning(sprintf("All %s_CTRL values are NA in plate %s", ctrl_type, plate_id))
      return(NULL)
    }
  }
  
  # Process controls
  tryCatch({
    pos_ctrl <- na.omit(remove_outliers(controls$POS))
    neg_ctrl <- na.omit(remove_outliers(controls$NEG))
    
    if (length(pos_ctrl) == 0 || length(neg_ctrl) == 0) {
      warning(sprintf("No valid controls after outlier removal in plate %s", plate_id))
      return(NULL)
    }
    
    avg_low <- mean(pos_ctrl)
    avg_high <- mean(neg_ctrl)
    
    if (avg_low >= avg_high) {
      warning(sprintf("Invalid control values: POS (%g) >= NEG (%g) in plate %s", 
                      avg_low, avg_high, plate_id))
      return(NULL)
    }
    
    # Calculate inhibition
    plate_data$inhibition_percent <- ((avg_high - plate_data$INTENSITY) / 
                                        (avg_high - avg_low)) * 100
    
    # Process plate data
    plate_data <- plate_data %>% 
      mutate(
        inhibition_percent = pmax(0, pmin(100, inhibition_percent)),
        Row = gsub("([A-Z]+).*", "\\1", WELL),
        Column = as.numeric(gsub("[A-Z]+", "", WELL))
      )
    
    # Create complete grid
    complete_wells <- expand.grid(
      Row = unique(plate_data$Row),
      Column = 1:max(plate_data$Column, na.rm = TRUE)
    ) %>%
      mutate(WELL = paste0(Row, Column))
    
    # Final processing
    processed_data <- complete_wells %>%
      left_join(plate_data, by = c("Row", "Column", "WELL")) %>%
      arrange(Row, Column) %>%
      mutate(
        Column = as.numeric(as.character(Column)),
        Row = factor(Row, levels = unique(Row))
      )
    
    # Calculate quality metrics
    quality_metrics <- list(
      zfactor = calc_zfactor(neg_ctrl, pos_ctrl),
      ssmd = calc_ssmd(neg_ctrl, pos_ctrl),
      robust_z_prime = calc_robust_zfactor(neg_ctrl, pos_ctrl),
      signal_vs_bg = round(mean(neg_ctrl) / mean(pos_ctrl), 1)
    )
    
    return(list(
      processed_data = processed_data,
      quality_metrics = quality_metrics
    ))
    
  }, error = function(e) {
    warning(sprintf("Error processing plate %s: %s", plate_id, e$message))
    return(NULL)
  })
}

#' Generate plate heatmap
#' @param plate_data Processed plate data
#' @param data_type String; type of data to visualize ("inhibition_percent" or "residuals_norm")
#' @return ggplot object
#' @export
create_plate_heatmap <- function(plate_data, data_type) {
  
  if(data_type == "inhibition_percent"){ 
    
    # Determine plate dimensions from data
    max_row <- max(which(create_plate_labels(32) %in% unique(plate_data$Row)))
    max_col <- max(plate_data$Column, na.rm = TRUE)
    
    ggplot(plate_data, aes(x = Column, y = Row, fill = .data[[data_type]])) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(data = subset(plate_data, DRUG_NAME == "NEG_CTRL"),
                aes(label = "x"), size = 3, color = "black") +
      geom_text(data = subset(plate_data, DRUG_NAME == "POS_CTRL"),
                aes(label = "*"), size = 5, color = "black") +
      scale_fill_gradient(
        low = "white",
        high = "#bb5788",
        name = "Inhibition\n(%)",
        limits = c(0, 100),
        na.value = "grey95"
      ) +
      scale_x_continuous(
        breaks = seq(1, max_col, by = ifelse(max_col > 24, 2, 1)),  # Adjust break intervals
        expand = c(0, 0)
      ) +
      scale_y_discrete(
        limits = rev(create_plate_labels(max_row)),
        expand = c(0, 0)
      ) +
      coord_fixed() +
      labs(
        title = sprintf("%dx%d Well Plate Inhibition Heatmap", max_row, max_col),
        subtitle = "Plate layout with controls",
        x = "Column",
        y = "Row"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        axis.text = element_text(size = ifelse(max_col > 24, 7, 8), color = "black"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        panel.grid = element_blank()
      )
  } else if(data_type == "residuals_norm"){
    
    # Determine plate dimensions from data
    max_row <- max(which(create_plate_labels(32) %in% unique(plate_data$Row)))
    max_col <- max(plate_data$Column, na.rm = TRUE)
    
    ggplot(plate_data, aes(x = Column, y = Row, fill = .data[[data_type]])) +
      geom_tile(color = "white", linewidth = 0.5) +
      scale_fill_gradient(
        low = "white",
        high = "#bb5788",
        name = "Inhibition\n(%)",
        limits = c(0, 100),
        na.value = "grey95"
      ) +
      scale_x_continuous(
        breaks = seq(1, max_col, by = ifelse(max_col > 24, 2, 1)),  # Adjust break intervals
        expand = c(0, 0)
      ) +
      scale_y_discrete(
        limits = rev(create_plate_labels(max_row)),
        expand = c(0, 0)
      ) +
      coord_fixed() +
      labs(
        title = sprintf("Cropped Plate Error Heatmap", max_row, max_col),
        subtitle = "Plate layout with controls",
        x = "Column",
        y = "Row"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        axis.text = element_text(size = ifelse(max_col > 24, 7, 8), color = "black"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        panel.grid = element_blank()
      )
  }
}


#' Create scatter plot of inhibition by row
#' @param plate_data Processed plate data
#' @return ggplot object
#' @export
create_row_scatter <- function(plate_data) {
  
  # Determine number of rows from data
  max_row <- max(which(create_plate_labels(32) %in% unique(plate_data$Row)))
  
  # Adjust point size based on number of wells
  point_size <- ifelse(max_row <= 8, 3,        # 96-well
                       ifelse(max_row <= 16, 2.5,      # 384-well
                              2))                       # 1536-well
  
  # Adjust jitter height based on plate size
  jitter_height <- ifelse(max_row <= 8, 0.3,    # 96-well
                          ifelse(max_row <= 16, 0.25,    # 384-well
                                 0.2))                    # 1536-well
  
  # Adjust text sizes based on number of rows
  text_size <- ifelse(max_row <= 8, 9,          # 96-well
                      ifelse(max_row <= 16, 8,         # 384-well
                             7))                       # 1536-well
  
  ggplot(plate_data, aes(x = inhibition_percent, y = Row)) +
    geom_point(aes(color = DRUG_NAME),
               size = point_size,
               alpha = 0.6,
               position = position_jitter(width = 0, height = jitter_height)) +
    scale_color_manual(
      values = c("POS_CTRL" = "#1B9E77", "NEG_CTRL" = "#D95F02", "TEST" = "#7570B3"),
      name = "Compound Type",
      labels = c("Positive Control", "Negative Control", "Test Compounds")
    ) +
    scale_x_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, by = 20)
    ) +
    scale_y_discrete(limits = rev(create_plate_labels(max_row))) +
    labs(
      x = "Inhibition (%)",
      y = "Row",
      title = sprintf("Distribution of Inhibition by Row (%d rows)", max_row)
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text = element_text(size = text_size),
      axis.title = element_text(size = 12, face = "bold"),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
}


#' Process dose response curves
#' @param plate_data Processed plate data
#' @param plot_dose_response Logical. If TRUE, generates dose-response curve plots
#' @param directory_curves String. Directory path for saving dose-response curve plots
#' @param cores_n Integer. Number of cores for parallel processing
#' @return Data frame with curve fitting results
#' @export
analyze_dose_response <- function(plate_data, plot_dose_response, directory_curves, cores_n) {
  # Filter data
  plate_data <- plate_data %>%
    filter(!is.na(DRUG_NAME), 
           !DRUG_NAME %in% c("NEG_CTRL", "POS_CTRL")) %>%
    group_by(DRUG_NAME) %>%
    filter(n_distinct(CONC) >= 5)
  
  # Analyze curves using parallel processing
  results <- mclapply(unique(plate_data$DRUG_NAME), function(drug) {
    tmp <- plate_data[plate_data$DRUG_NAME == drug, ]
    
    suppressWarnings(
      curve_results <- calculate_all(tmp$CONC, tmp$inhibition_percent,
                                     plate_data$BARCODE[[1]], drug, WELL = tmp$WELL,
                                     plot_dose_response = plot_dose_response,
                                     directory_curves = directory_curves)
    )
    
    # Merge results
    tmp <- merge(
      curve_results[, c("dose", "residuals", "WELL", "fitted")],
      tmp,
      by = "WELL"
    )
    
    # Calculate normalized residuals
    tmp$fitted_proportion <- tmp$fitted / 100
    tmp$normalization_factor <- 1 + (tmp$fitted_proportion * 
                                       (1 - tmp$fitted_proportion) / 0.25)
    tmp$residuals <- abs(tmp$residuals)
    tmp$residuals_norm <- tmp$residuals * tmp$normalization_factor
    
    return(tmp)
  }, mc.cores = cores_n)
  
  do.call(rbind, results)
}



#' Fit dose-response curves and get residuals
#' @param dose Numeric vector of drug concentrations
#' @param inhibition_percent Numeric vector of inhibition percentages
#' @param dataset Character vector identifying different datasets
#' @param drug_name Character string of drug name
#' @param WELL Character vector of well identifiers
#' @param plot_dose_response Logical; if TRUE, generates dose-response curve plots
#' @param directory_curves String; directory path for saving dose-response curve plots
#' @return Data frame with fitted values and analysis metrics
#' @export
calculate_all <- function(dose, inhibition_percent, dataset, drug_name, WELL = WELL, plot_dose_response = plot_dose_response, directory_curves = directory_curves) {
  if(all(inhibition_percent <= 0)) inhibition_percent <- rep(0, length(inhibition_percent))
  if(any(duplicated(inhibition_percent))) inhibition_percent <- seq(from = 0, length.out = length(inhibition_percent), by = 0.01) + inhibition_percent
  
  mat_tbl <- data.frame(inhibition_percent, dose, logconc = log10(dose), 
                        100-inhibition_percent, dataset = dataset, WELL = WELL)
  mat_tbl <- mat_tbl[order(mat_tbl[,"dose"]),]
  
  estimate_param <- tryCatch({
    drm(inhibition_percent ~ logconc, data = mat_tbl, 
        fct = LL.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50")),
        logDose = 10, control = drmc(errorm = FALSE))
  }, warning = function(w) {
    drm(inhibition_percent ~ logconc, data = mat_tbl, 
        fct = L.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50")),
        logDose = 10)
  }, error = function(e) {
    drm(inhibition_percent ~ logconc, data = mat_tbl, 
        fct = L.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50")),
        logDose = 10)
  })
  
  coef_estim <- coef(estimate_param)
  names(coef_estim) <- c("SLOPE","MIN","MAX","IC50")
  coef_estim["SLOPE"] <- coef_estim["SLOPE"] * -1
  
  coef_estim["IC50"] <- ifelse(coef_estim["MAX"] <= coef_estim["MIN"] | 
                                 coef_estim["IC50"] > max(mat_tbl$dose, na.rm = TRUE), 
                               max(mat_tbl$dose, na.rm = TRUE), coef_estim["IC50"])
  coef_estim["IC50"] <- ifelse(coef_estim["IC50"] < 0, 
                               min(mat_tbl$dose, na.rm = TRUE), coef_estim["IC50"])
  coef_estim["IC50"] <- ifelse(coef_estim["IC50"] < 0, 
                               mean(mat_tbl$dose, na.rm = TRUE), coef_estim["IC50"])
  coef_estim["IC50"] <- log10(coef_estim["IC50"])
  coef_estim["IC50"] <- ifelse(coef_estim["IC50"] < min(mat_tbl$logconc), 
                               max(mat_tbl$logconc), coef_estim["IC50"])
  coef_estim["IC50"] <- ifelse(all(mat_tbl$inhibition_percent < 0), 
                               max(mat_tbl$logconc, na.rm = TRUE), coef_estim["IC50"])
  
  coef_estim["MIN"] <- 0
  coef_estim["MAX"] <- max(mat_tbl$inhibition_percent, na.rm = TRUE)
  
  min_lower <- ifelse(min(mat_tbl$inhibition_percent, na.rm = TRUE) > 0,
                      min(mat_tbl$inhibition_percent, na.rm = TRUE), 0)
  min_lower <- ifelse(min_lower >= 100, 99, min_lower)
  
  coef_estim["MAX"] <- ifelse(coef_estim["MAX"] > 100, 100, coef_estim["MAX"])
  coef_estim["MAX"] <- ifelse(coef_estim["MAX"] < 0, 100, coef_estim["MAX"])
  
  max_lower <- ifelse(max(mat_tbl$inhibition_percent, na.rm = TRUE) > 100,
                      coef_estim["MAX"], max(mat_tbl$inhibition_percent, na.rm = TRUE))
  max_lower <- ifelse(max_lower < 0, coef_estim["MAX"], max(mat_tbl$inhibition_percent, na.rm = TRUE))
  max_lower <- ifelse(max_lower < 0, 0, max_lower)
  max_lower <- ifelse(max_lower > 100, 100, max_lower)
  
  run_avg <- caTools::runmean(mat_tbl$inhibition_percent, 10)
  max_upper <- ifelse(any(run_avg[-nrow(mat_tbl)] > run_avg[nrow(mat_tbl)]),
                      max(mat_tbl$inhibition_percent[run_avg > run_avg[nrow(mat_tbl)]]),
                      coef_estim["MAX"])
  max_upper <- ifelse(any(mat_tbl$inhibition_percent > max_upper),
                      mean(mat_tbl$inhibition_percent[mat_tbl$inhibition_percent > max_upper]) + 5,
                      max_upper)
  max_upper <- ifelse(max_upper < 0, coef_estim["MAX"], max_upper)
  max_upper <- ifelse(max_upper > 100, 100, max_upper)
  max_upper <- ifelse(max_lower > max_upper, coef_estim["MAX"], max_upper)
  
  mean_inh_last = mean(tail(mat_tbl$inhibition_percent, 2), na.rm = TRUE)
  if(mean_inh_last < 60) {
    if(mean_inh_last > 25) coef_estim["IC50"] <- mean(mat_tbl$logconc, na.rm = TRUE)
    else if(mean_inh_last < 25) coef_estim["IC50"] <- max(mat_tbl$logconc, na.rm = TRUE)
  }
  if(mean(mat_tbl$inhibition_percent[1:3], na.rm = TRUE) < 5) {
    coef_estim["IC50"] <- max(mat_tbl$logconc, na.rm = TRUE)
  }
  if(unname(coef_estim["MIN"]) == unname(coef_estim["MAX"])) {
    coef_estim["MAX"] <- coef_estim["MAX"] + 0.001
  }
  
  nls_result_ic50_old <- function() {
    tryCatch({
      nls(inhibition_percent ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))),
          data = mat_tbl, algorithm = "port",
          start = list(SLOPE = 1, MIN = coef_estim["MIN"][[1]], 
                       MAX = coef_estim["MAX"][[1]], IC50 = coef_estim["IC50"][[1]]),
          lower = list(SLOPE = 0.1, MIN = 0, MAX = max_lower, 
                       IC50 = min(mat_tbl$logconc)),
          upper = list(SLOPE = 5, MIN = 0, MAX = max_upper, 
                       IC50 = max(mat_tbl$logconc)),
          control = list(warnOnly = TRUE, minFactor = 1/2048))
    }, error = function(e) {
      minpack.lm::nlsLM(inhibition_percent ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))),
                        data = mat_tbl,
                        start = list(SLOPE = 1, MIN = coef_estim["MIN"][[1]],
                                     MAX = coef_estim["MAX"][[1]], IC50 = coef_estim["IC50"][[1]]),
                        lower = c(SLOPE = 0.1, MIN = 0, MAX = max_lower,
                                  IC50 = min(mat_tbl$logconc)),
                        upper = c(SLOPE = 5, MIN = 0, MAX = max_upper,
                                  IC50 = max(mat_tbl$logconc)))
    })
  }
  
  nls_result_ic50 <- nls_result_ic50_old()
  
  nls_result_ic50_2 <- tryCatch({
    nls(inhibition_percent ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))),
        data = mat_tbl, algorithm = "port",
        start = list(SLOPE = 1, MIN = coef_estim["MIN"][[1]],
                     MAX = coef_estim["MAX"][[1]], IC50 = median(mat_tbl$logconc)),
        lower = list(SLOPE = 0.1, MIN = 0, MAX = max_lower,
                     IC50 = min(mat_tbl$logconc)),
        upper = list(SLOPE = 5, MIN = 0, MAX = max_upper,
                     IC50 = max(mat_tbl$logconc)),
        control = list(warnOnly = TRUE, minFactor = 1/2048))
  }, warning = function(w) {
    nls_result_ic50_old()
  }, error = function(e) {
    nls_result_ic50_old()
  })
  
  nls_result_ic50 = tryCatch({
    summary(nls_result_ic50)
    nls_result_ic50
  }, error = function(e) {
    nls_result_ic50_2
  })
  
  sumIC50 = list(summary(nls_result_ic50), summary(nls_result_ic50_2))
  
  ic50std_resid <- round(sqrt(sum((sumIC50[[1]]$residuals)^2)/
                                (length(sumIC50[[1]]$residuals)-1)), 1)
  ic50std_resid2 <- round(sqrt(sum((sumIC50[[2]]$residuals)^2)/
                                 (length(sumIC50[[2]]$residuals)-1)), 1)
  
  switch_ = which.min(c(ic50std_resid, ic50std_resid2))
  nls_result_ic50 = list(nls_result_ic50, nls_result_ic50_2)[[switch_]]
  
  sumIC50 = summary(nls_result_ic50)
  ic50std_Error <- sumIC50$coefficients["IC50","Std. Error"]
  ic50std_resid <- round(sqrt(sum((sumIC50$residuals)^2)/(length(sumIC50$residuals)-1)), 1)
  max_signal <- max(mat_tbl$dose, na.rm = TRUE)
  min_signal <- min(mat_tbl$dose, na.rm = TRUE)
  mat_tbl$residuals = sumIC50$residuals
  mat_tbl$fitted = fitted(nls_result_ic50)
  
  coef_ic50 <- coef(nls_result_ic50)[c("IC50", "SLOPE","MAX","MIN")]
  coef_ic50["IC50"] <- 10^coef_ic50["IC50"]
  coef_ic50["IC50"] <- ifelse(coef_ic50["SLOPE"] < 0, max_signal, coef_ic50["IC50"])
  coef_ic50["IC50"] <- ifelse(coef_ic50["MAX"] < 0, max_signal, coef_ic50["IC50"])
  coef_ic50["IC50"] <- ifelse(coef_ic50["MAX"] < 10, max_signal, coef_ic50["IC50"])
  coef_ic50["MAX"] <- ifelse(coef_ic50["MAX"] < 0, 0, coef_ic50["MAX"])
  coef_ic50["IC50"] <- ifelse(all(c(max(mat_tbl$inhibition_percent, na.rm = TRUE),
                                    min(mat_tbl$inhibition_percent, na.rm = TRUE)) > 50),
                              min_signal, coef_ic50["IC50"])
  
  x <- seq(min(mat_tbl$logconc), max(mat_tbl$logconc), length = 100)
  yic <- predict(nls_result_ic50, data.frame(logconc = x))
  
  if(plot_dose_response){
    icpl <- ggplot2::ggplot(mat_tbl, aes(logconc, inhibition_percent)) +
      scale_x_continuous(breaks = mat_tbl$logconc, labels = mat_tbl$dose) +
      geom_point(size = 2.8) +
      geom_line(data = data.frame(x = x, y = yic), aes(x, yic),
                color = "steelblue", size = 0.8) +
      geom_vline(xintercept = log10(coef_ic50["IC50"]), colour = "grey", size = 0.8) +
      ggtitle(paste0(strtrim(drug_name, 15), " (SE:", as.numeric(ic50std_resid), ")\n")) +
      theme_bw() +
      labs(y = "% inhibition", x = "conc(nM)") +
      ylim(-25, 125) +
      geom_text(mapping = aes(x2, y2, label = text2),
                data = data.frame(x2 = log10(coef_ic50["IC50"]) * 0.95,
                                  y2 = 115, text2 = "IC50"),
                color = "grey", parse = TRUE) +
      theme(plot.background = element_rect(fill = "transparent", colour = NA),
            panel.background = element_rect(fill = "transparent", colour = NA),
            plot.title = element_text(hjust = 0.5, size = 12.5))
    
    # Check directory path
    clean_dir <- normalizePath(directory_curves, winslash = "/", mustWork = FALSE)
    if (!dir.exists(clean_dir)) {
      warning(sprintf("Directory '%s' does not exist. Creating it now.", clean_dir))
      dir.create(clean_dir, recursive = TRUE)
    }
    
    # Save the plot
    ggsave(
      filename = file.path(clean_dir, paste0(drug_name, "_", dataset, ".jpg")),
      plot = icpl,
      width = 500/72,
      height = 500/72,
      dpi = 300
    )
  }
  return(mat_tbl)
}




#' Create combined visualization
#' @param plate_data Processed plate data with residuals
#' @param plate_data_raw Raw processed plate data with controls preserved
#' @param file_path Output file path
#' @return None, saves plot to file
#' @export
create_combined_plot <- function(plate_data, plate_data_raw, file_path) {
  
  suppressWarnings({
    heatmap <- create_plate_heatmap(plate_data_raw, "inhibition_percent")
    heatmap2 <- create_plate_heatmap(plate_data, "residuals_norm")
    scatter <- create_row_scatter(plate_data_raw)
    
    layout_matrix <- matrix(c(1, 3,  # First column: heatmap and heatmap2
                              2, 2),   # Second column: scatter spans both rows
                            nrow = 2, byrow = FALSE)
    
    ggsave(file_path,
           plot = arrangeGrob(
             heatmap + theme(legend.position = "right"),
             scatter + theme(legend.position = "bottom"),
             heatmap2 + theme(legend.position = "right"),
             layout_matrix = layout_matrix,
             widths = c(1.2, 1),
             heights = c(1, 1)
           ),
           width = 20,
           height = 10,
           dpi = 300,
           units = "in")
  })
}



#' Process Plate Data for Drug Response Analysis
#' 
#' @description
#' Analyzes high-throughput screening plate data for drug response experiments, calculating 
#' quality control metrics including Normalized Residual Fit Error (NRFE), Z-factor, and 
#' Strictly Standardized Mean Difference (SSMD). The function processes plates identified 
#' by unique barcodes and can generate both dose-response curve visualizations and 
#' plate summary plots.
#'
#' @param plate_data A data frame containing plate assay data with required columns:
#'   \itemize{
#'     \item BARCODE: Unique identifier for each plate
#'     \item DRUG_NAME: Name of the drug or control (including "POS_CTRL" and "NEG_CTRL")
#'     \item CONC: Drug concentration in nM
#'     \item INTENSITY: Measured response intensity
#'     \item WELL: Well position identifier (e.g., "A1", "B2"). Optional for core analysis; Required for visualizations
#'   }
#' @param plot_dose_response Logical. If TRUE, generates dose-response curve plots. 
#'   Default is FALSE.
#' @param directory_curves Character string. Directory path for saving dose-response 
#'   curve plots. Default is "./fits".
#' @param plot_plate_summary Logical. If TRUE, generates summary plots for each plate including
#'   inhibition heatmap, residuals heatmap, and row-wise distribution. Default is FALSE.
#' @param remove_empty_wells Logical. If TRUE, removes wells with missing or empty 
#'   concentration values. Default is TRUE.
#' @param verbose Logical. If TRUE, prints processing status messages. Default is FALSE.
#' @param cores_n Integer. Number of cores to use for parallel processing of dose-response
#'   curves. Default is 1.
#'
#' @return Returns an object of class "plate_analysis" containing:
#'   \itemize{
#'     \item plate_statistics: Data frame with quality metrics for each plate including:
#'       \itemize{
#'         \item zfactor: Z-factor quality metric
#'         \item ssmd: Strictly Standardized Mean Difference
#'         \item robust_z_prime: Robust version of Z-factor using median and MAD
#'         \item signal_vs_bg: Signal to background ratio
#'         \item NRFE: Normalized Residual Fit Error
#'       }
#'     \item detailed_results: List of processed results for each plate containing:
#'       \itemize{
#'         \item processed_results: Processed plate data with calculated metrics
#'         \item dose_response_results: Dose-response curve fitting results
#'         \item quality_metrics: Plate-specific quality metrics
#'       }
#'     \item metadata: Processing information including:
#'       \itemize{
#'         \item processed_at: Timestamp of analysis
#'         \item total_plates: Number of plates in input
#'         \item successful_plates: Number of successfully processed plates
#'       }
#'   }
#'   Returns NULL if no plates were successfully processed.
#'
#' @importFrom dplyr filter %>%
#' @importFrom stats na.omit
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(plate_data)
#' 
#' # Basic usage - minimal processing
#' results <- process_plate_data(plate_data)
#' 
#' # Full analysis with all visualizations
#' results <- process_plate_data(
#'   plate_data,
#'   directory_curves = "./fits", 
#'   plot_dose_response = TRUE,
#'   plot_plate_summary = TRUE,
#'   verbose = TRUE,
#'   cores_n = 4
#' )
#' 
#' # Access quality metrics
#' plate_stats <- results$plate_statistics
#' print(plate_stats)
#' 
#' # Access processed data for a specific plate
#' first_plate <- results$detailed_results[[1]]
#' }
process_plate_data <- function(plate_data, 
                               plot_dose_response = FALSE,
                               directory_curves = "./fits", 
                               plot_plate_summary = FALSE,
                               remove_empty_wells = TRUE, 
                               verbose = FALSE,
                               cores_n = 1) {
  
  # Input validation
  if (!is.data.frame(plate_data)) {
    stop("Input must be a data frame", call. = FALSE)
  }
  
  # Define required columns without INTENSITY
  base_required_cols <- c("BARCODE", "DRUG_NAME", "WELL", "CONC")
  
  # Handle missing WELL column
  if (!"WELL" %in% colnames(plate_data)) {
    plate_data$WELL <- rep("A1", nrow(plate_data))
  }
  
  # Check if either INTENSITY or inhibition_percent is present
  has_intensity <- "INTENSITY" %in% colnames(plate_data)
  has_inhibition <- "inhibition_percent" %in% colnames(plate_data)
  
  if (!has_intensity && !has_inhibition) {
    stop("Either INTENSITY or inhibition_percent column must be present", call. = FALSE)
  }
  
  # Check other required columns
  missing_cols <- setdiff(base_required_cols, colnames(plate_data))
  if (length(missing_cols) > 0) {
    stop(
      sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }
  
  # Directory validation and creation
  if (plot_dose_response) {
    if (!dir.exists(directory_curves)) {
      dir.create(directory_curves, recursive = TRUE)
    }
    if (!dir.exists(directory_curves)) {
      stop(
        sprintf("Failed to create directory: %s", directory_curves),
        call. = FALSE
      )
    }
  }
  
  # Remove empty wells if specified
  if (remove_empty_wells) {
    plate_data <- plate_data %>%
      filter(!(is.null(CONC) | is.na(CONC) | CONC == "") | 
               DRUG_NAME %in% c("POS_CTRL", "NEG_CTRL"))
  }
  
  # Initialize results containers
  results_list <- list()
  all_plate_stats <- data.frame()
  
  # Get unique barcodes
  barcodes <- unique(plate_data$BARCODE)
  if (length(barcodes) == 0) {
    warning("No valid plates found in the data", call. = FALSE)
    return(NULL)
  }
  
  # Process each plate
  for (barcode in barcodes) {
    tryCatch({
      if (verbose) message(sprintf("\nAnalyzing %s plate:", barcode))
      
      # Get data for single plate
      current_plate <- plate_data[plate_data$BARCODE == barcode, ]
      
      # Validate plate data
      if (nrow(current_plate) == 0) {
        if (verbose) message("  Skipping - empty plate")
        warning(sprintf("Empty plate data for barcode: %s", barcode), call. = FALSE)
        next
      }
      
      if (verbose) message("  Processing controls...")
      processed_results <- process_plate(current_plate)
      if (is.null(processed_results)) {
        if (verbose) message(" Failed")
        next
      }
      if (verbose) message(" Done")
      
      if (verbose) message("  Fitting dose response curves...")
      dose_response_results <- analyze_dose_response(
        processed_results$processed_data, 
        plot_dose_response, 
        directory_curves,
        cores_n
      )
      if (verbose) message(" Done")
      
      if (verbose) message("  Computing quality metrics...")
      plate_stats <- processed_results$quality_metrics
      plate_stats$barcode <- barcode
      plate_stats$NRFE <- mean(dose_response_results$residuals_norm, na.rm = TRUE)
      if (verbose) message(" Done")
      
      # Generate plate summary visualization
      if (plot_plate_summary & !all(current_plate$WELL=="A1")) {
        if (verbose) message("  Creating summary plots...")
        output_path <- file.path(".", sprintf("plate_analysis_%s.pdf", barcode))
        tryCatch({
          create_combined_plot(dose_response_results, processed_results$processed_data, output_path)
          if (verbose) message(" Done")
        }, error = function(e) {
          if (verbose) message(" Failed")
          warning(
            sprintf("Failed to create plot for plate %s: %s", barcode, e$message),
            call. = FALSE
          )
        })
      }
      
      # Store results
      all_plate_stats <- rbind(all_plate_stats, plate_stats)
      results_list[[barcode]] <- list(
        processed_results = processed_results,
        dose_response_results = dose_response_results,
        quality_metrics = plate_stats
      )
      
      if (verbose) message(sprintf("  Plate %s completed successfully", barcode))
      
    }, error = function(e) {
      if (verbose) message(sprintf("  Error: %s", e$message))
      warning(sprintf("Error processing plate %s: %s", barcode, e$message), call. = FALSE)
    })
  }
  
  # Check if any plates were successfully processed
  if (length(results_list) == 0) {
    warning("No plates were successfully processed", call. = FALSE)
    return(NULL)
  }
  
  if (verbose) {
    message(sprintf(
      "\nProcessed %d/%d plates successfully",
      length(results_list),
      length(barcodes)
    ))
  }
  
  # Add processing metadata
  result <- list(
    plate_statistics = all_plate_stats,
    detailed_results = results_list,
    metadata = list(
      processed_at = Sys.time(),
      total_plates = length(barcodes),
      successful_plates = length(results_list),
      remove_empty_wells = remove_empty_wells
    )
  )
  
  class(result) <- c("plate_analysis", class(result))
  return(result)
}
