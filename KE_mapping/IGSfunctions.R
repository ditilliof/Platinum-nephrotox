if (!require("pacman", quietly = T)) {
  install.packages("pacman")
}
pacman::p_load(tidyverse, corrplot, readr)

compute_scorenum_vitro <- function(data, pear, genes) {
  if (length(data$l2fc) == 0 || length(genes) == 0) {
    return(0)
  } else {
    score <- data$l2fc * (1 / nrow(pear)) * sum(pear[, genes])
    return(score)
  }
}

compute_weight_vitro <- function(data, pear, genes) {
  if (length(data$l2fc) == 0 || length(genes) == 0) {
    return(0)
  } else {
    weight <- (1 / nrow(pear)) * sum(pear[, genes])
    return(weight)
  }
}


compute_score_vitro <- function(data, pear, genes, label) {
  scores <- data.frame(Score = numeric(),
                       REPLICATE = integer(),
                       TIME = numeric(),
                       StateVar = character(),
                       DOSE = numeric()
  )
  for (dose in c(0.1, 0.5, 1, 2.5, 5, 10, 20, 30, 50)) {
    data_dose <- data %>% dplyr::filter(DOSE == dose)
    for (replicate in unique(data_dose$REPLICATE)) {
      for (time in unique(data_dose$TIME)) {
        score_num <- sum(sapply(seq_along(genes), function(i) {
          compute_scorenum_vitro(data_dose[data_dose$TIME == time & data_dose$gene_symbol == genes[i] & data_dose$REPLICATE == replicate,], pear, genes[i])
        }))
        score_den <- sum(sapply(seq_along(genes), function(i) {
          compute_weight_vitro(data_dose[data_dose$TIME == time & data_dose$gene_symbol == genes[i] & data_dose$REPLICATE == replicate,],pear, genes[i])
        }))
        scores[nrow(scores) + 1,1] <- c(score_num/score_den)
        scores[nrow(scores),2] <- c(replicate)
        scores[nrow(scores),3] <- c(time)
        scores[nrow(scores),4] <- c(label)
        scores[nrow(scores),5] <- c(dose)
      }
    }
  }
  return(scores)
}

compute_score_woutlabel_vitro <- function(data, pear, genes) {
  scores <- data.frame(Score = numeric(),
                       REPLICATE = integer(),
                       TIME = numeric(),
                       DOSE = numeric()
  )
  if (length(data$l2fc) == 0 || length(genes) == 0) {
    return(0)
  } else {
    for (dose in c(0.1, 0.5, 1, 2.5, 5, 10, 20, 30, 50)) {
      data_dose <- data %>% dplyr::filter(DOSE == dose)
      for (replicate in unique(data_dose$REPLICATE)) {
        for (time in unique(data_dose$TIME)) {
          score_num <- sum(sapply(seq_along(genes), function(i) {
            compute_scorenum_vitro(data_dose[data_dose$TIME == time & data_dose$gene_symbol == genes[i] & data_dose$REPLICATE == replicate,], pear, genes[i])
          }))
          score_den <- sum(sapply(seq_along(genes), function(i) {
            compute_weight_vitro(data_dose[data_dose$TIME == time & data_dose$gene_symbol == genes[i] & data_dose$REPLICATE == replicate,],pear, genes[i])
          }))
          scores[nrow(scores) + 1,1] <- c(score_num/score_den)
          scores[nrow(scores),2] <- c(replicate)
          scores[nrow(scores),3] <- c(time)
          scores[nrow(scores),4] <- c(dose)
        }
      }
    }
    return(scores)
  }
}

compute_scorenum_vivo <- function(data, pear, genes) {
  if (length(data$log2FoldChange) == 0 || length(genes) == 0) {
    return(0)
  } else {
    min_row_index <- which.min(data$pvalue)
    score <- (data[min_row_index,])$log2FoldChange * (1 / nrow(pear)) * sum(pear[, genes])
    return(score)
  }
}

#(data$log2FoldChange / nrow(pear)) in the original formula

compute_weight_vivo <- function(data, pear, genes) {
  if (length(data$log2FoldChange) == 0 || length(genes) == 0) {
    return(0)
  } else {
    weight <- (1 / nrow(pear)) * sum(pear[, genes])
    return(weight)
  }
}

compute_score_vivo <- function(data, pear, genes, label) {
  scores <- data.frame(Score = numeric(),
                       REPLICATE = integer(),
                       TIMEPOINT = numeric(),
                       StateVar = character(),
                       CONCENTRATION = numeric()
  )
  for (dose in unique(data$CONCENTRATION)) {
    data_dose <- data %>% filter(CONCENTRATION == dose)
    for (replicate in unique(data_dose$REPLICATE)) {
      for (time in unique((data_dose%>%filter(GeneSymbol %in% genes, REPLICATE == replicate))$TIMEPOINT)) {
        score_num <- sum(sapply(seq_along(genes), function(i) {
          compute_scorenum_vivo(data_dose[data_dose$TIMEPOINT == time & data_dose$GeneSymbol == genes[i] & data_dose$REPLICATE == replicate,], pear, genes[i])
        }))
        score_den <- sum(sapply(seq_along(genes), function(i) {
          compute_weight_vivo(data_dose[data_dose$TIMEPOINT == time & data_dose$GeneSymbol == genes[i] & data_dose$REPLICATE == replicate,],pear, genes[i])
        }))
        scores[nrow(scores) + 1,1] <- c(score_num/score_den)
        scores[nrow(scores),2] <- c(replicate)
        scores[nrow(scores),3] <- c(time)
        scores[nrow(scores),4] <- c(label)
        scores[nrow(scores),5] <- c(dose)
      }
    }
  }
  return(scores)
}

compute_score_woutlabel_vivo <- function(data, pear, genes) {
  scores <- data.frame(Score = numeric(),
                       REPLICATE = integer(),
                       TIMEPOINT = numeric(),
                       CONCENTRATION = numeric()
  )
  if (length(data$log2FoldChange) == 0 || length(genes) == 0) {
    return(0)
  } else {
    for (dose in unique(data$CONCENTRATION)) {
      data_dose <- data %>% filter(CONCENTRATION == dose)
      for (replicate in unique(data_dose$REPLICATE)) {
        for (time in unique(data_dose$TIMEPOINT)) {
          score_num <- sum(sapply(seq_along(genes), function(i) {
            compute_scorenum_vivo(data_dose[data_dose$TIMEPOINT == time & data_dose$GeneSymbol == genes[i] & data_dose$REPLICATE == replicate,], pear, genes[i])
          }))
          score_den <- sum(sapply(seq_along(genes), function(i) {
            compute_weight_vivo(data_dose[data_dose$TIMEPOINT == time & data_dose$GeneSymbol == genes[i] & data_dose$REPLICATE == replicate,],pear, genes[i])
          }))
          scores[nrow(scores) + 1,1] <- c(score_num/score_den)
          scores[nrow(scores),2] <- c(replicate)
          scores[nrow(scores),3] <- c(time)
          scores[nrow(scores),4] <- c(dose)
        }
      }
    }
    return(scores)
  }
}

process_histopathology <- function(df, header_row, statevar, concentration = 5) {
  # header_row has: "Grade / Day", 1, 2, 4, ..., 29
  time_labels <- as.character(header_row[1, -1])  # only the numeric timepoints
  
  # 1 grade column + timepoints must equal ncol(df)
  stopifnot(ncol(df) == length(time_labels) + 1)
  
  # Rename data columns to: Grade, 1, 2, 4, ..., 29
  names(df) <- c("Grade", time_labels)
  
  df %>%
    pivot_longer(
      cols      = all_of(time_labels),
      names_to  = "TIMEPOINT",
      values_to = "score"
    ) %>%
    mutate(
      Grade_num = case_when(
        Grade == "absent"   ~ 0,
        Grade == "minimal"  ~ 0.125,
        Grade == "mild"     ~ 0.375,
        Grade == "moderate" ~ 0.625,
        Grade == "marked"   ~ 0.875,
        TRUE                ~ NA_real_
      ),
      TIMEPOINT = as.numeric(TIMEPOINT)
    ) %>%
    uncount(weights = score) %>%
    group_by(TIMEPOINT) %>%
    summarise(
      means = mean(Grade_num, na.rm = TRUE),
      sds   = sd(Grade_num,   na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      StateVar      = statevar,
      CONCENTRATION = concentration
    )
}

filter_and_summ = function(df, module_num) {
  df_name <- deparse(substitute(df))  # Capture the original name
  
  df = df %>% 
    filter(module == paste0('rKIDNEY_', module_num)) %>%
    group_by(time, DOSE) %>%
    summarise(meaneg = mean(eg_score), sdeg = sd(eg_score), .groups = "drop")
  
  df$DOSE <- factor(df$DOSE, levels = unique(df$DOSE), 
                    labels = case_when(
                      grepl("vitro", df_name) ~ paste0(unique(df$DOSE), " \u03BCM"),
                      TRUE                   ~ paste0(5, "mg/kg")
                    )
  )
  return(df)
}


plot_vitro = function(module_num){
  tiff(filename = paste0(PlotOXdir,"/RKID",module_num,"vitro.tiff"), width = 15, height = 7, units = "in", res = 700)
  
  print(ggplot(get(paste0("EGs_vitro_",module_num,"_summ")), aes(x = time, y = meaneg)) +
          geom_errorbar(aes(ymin = meaneg - sdeg, ymax = meaneg + sdeg), width = 1) +
          geom_line() +
          geom_point() +
          facet_wrap(~DOSE, labeller = labeller(DOSE = dose_labels), nrow = 2, ncol = 5) +  # Custom labels applied here
          theme_classic() + 
          labs(
            title = paste0("Eigengene Score of module rKID",module_num), 
            subtitle = "Uploaded data: in vitro RPTEC/TERT1", 
            x = "Time (h)", 
            y = "Eigengene Score"
          ) + 
          base_theme + 
          theme(
            aspect.ratio = 1,  # Keep square facets
            panel.spacing = unit(2, "lines"),  # Increase spacing between facets
            plot.margin = margin(50, 50, 50, 50),  # Increase plot margin
            plot.title = element_text(size = 24, face = "bold"),  # Increase title size
            plot.subtitle = element_text(size = 20),  # Increase subtitle size
            axis.text.x = element_text(angle = 45, hjust = 1, size = 14)
          ))
  print("Saving TIFF now...")
  # Close the TIFF device to save the file
  dev.off()
  print("TIFF saved successfully.")
}

plot_vitro_single_dose <- function(module_num, dose_value) {
  df_name <- paste0("EGs_vitro_", module_num, "_summ")
  df <- get(df_name)
  
  dose_label <- paste0(dose_value, " μM") 
  df_filtered <- df %>% filter(DOSE == dose_label)
  
  if (nrow(df_filtered) == 0) {
    warning(paste("No data found for module", module_num, "at dose", dose_label))
    return(invisible(NULL))
  }
  
  filename <- paste0(PlotOXdir,"/RKID", module_num, "_vitro_dose_", dose_value, "uM.tiff")
  tiff(filename = filename, width = 7, height = 7, units = "in", res = 700)
  
  print(
    ggplot(df_filtered, aes(x = time, y = meaneg)) +
      geom_errorbar(aes(ymin = meaneg - sdeg, ymax = sdeg + meaneg), width = 1) +
      geom_line() +
      geom_point() +
      theme_classic() +
      labs(
        title = paste0("Eigengene Score (Module rKID", module_num, ", Dose ", dose_label, ")"),
        subtitle = "in vitro RPTEC/TERT1",
        x = "Time (h)",
        y = "Eigengene Score"
      ) +
      base_theme +
      theme(
        plot.title = element_text(size = 20, face = "bold"),
        plot.subtitle = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14)
      )
  )
  
  dev.off()
  print(paste("TIFF saved to", filename))
}


plot_vivo = function(module_num){
  tiff(filename = paste0(PlotOXdir,"/RKID",module_num,"vivo.tiff"), 
       width = 15, height = 4, units = "in", res = 700)
  
  print(
    ggplot(get(paste0("EGs_vivo_",module_num,"_summ")), 
           aes(x = time, y = meaneg)) +
      geom_errorbar(aes(ymin = meaneg - sdeg, ymax = meaneg + sdeg), width = 1) +
      geom_line() +
      geom_point() +
      theme_classic() + 
      labs(
        title = paste0("Eigengene Score of module rKID",module_num), 
        subtitle = "Uploaded data: in vivo rat kidney PPT", 
        x = "Time (h)", 
        y = "Eigengene Score"
      ) + 
      base_theme + 
      theme(
        plot.title = element_text(size = 24, face = "bold"),
        plot.subtitle = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
        axis.text.y = element_text(size = 20)
      )
  )
  
  dev.off()
}

plot_oxgenes_tiff <- function(data, xvar, yvar, ymin, ymax, colorvar,
                      title, subtitle = NULL, xlabel, ylabel,
                      filename, width = 7, height = 7, facet_var = NULL,
                      base_theme_override = NULL) {
  
  tiff(filename = filename, width = width, height = height, units = "in", res = 700)
  
  p <- ggplot(data, aes_string(x = xvar, y = yvar, color = colorvar)) +
    geom_line() +
    geom_errorbar(aes_string(ymin = ymin, ymax = ymax), width = 0.2) +
    labs(title = title, subtitle = subtitle, x = xlabel, y = ylabel, color = colorvar) +
    theme_minimal()
  
  if (!is.null(facet_var)) {
    p <- p + facet_wrap(as.formula(paste("~", facet_var)), scales = "free")
  }
  
  if (!is.null(base_theme_override)) {
    p <- p + base_theme_override
  }
  
  print(p)
  dev.off()
}

