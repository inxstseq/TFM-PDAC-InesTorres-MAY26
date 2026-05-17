# --- 09_analisis_supervivencia.R ---
# Autor: Ines Torres
# Fecha: 17/02/2026
# Kaplan-Meier y regresion de Cox para genes clave

# --- librerias ---
if (!require("survival")) install.packages("survival")
if (!require("survminer")) install.packages("survminer")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("gridExtra")) install.packages("gridExtra")

library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(gridExtra)

# --- configuracion ---
setwd("~/TFM_PDAC")

out_dir <- "06_integracion/survival"
fig_dir <- "06_integracion/survival/figures"
res_dir <- "resultados_TFM_v2"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

# --- cargar datos ---
key_genes <- read.csv("06_integracion/results/genes_clave_DEG_intersect_Hub.csv",
                       stringsAsFactors = FALSE)

matriz_vst <- read.table(
  "02_preprocesamiento/outputs_FINAL/matriz_vst_normalizada.txt",
  header = TRUE, row.names = 1, sep = "\t"
)

traits <- read.table(
  "01_datos/clinical/TCGA_PAAD_traits_numeric.txt",
  header = TRUE, row.names = 1, sep = "\t"
)

cat("[OK] Datos cargados\n")
cat("  Genes clave:", nrow(key_genes), "\n")
cat("  Matriz:", nrow(matriz_vst), "genes x", ncol(matriz_vst), "muestras\n")
cat("  Muestras clinicas:", nrow(traits), "\n\n")

# --- preparar datos de supervivencia ---
rownames(traits) <- gsub("-", ".", rownames(traits))

valid_survival <- !is.na(traits$survival) & !is.na(traits$death) & traits$survival > 0
traits_surv <- traits[valid_survival, ]

common_samples <- intersect(rownames(traits_surv), colnames(matriz_vst))

traits_final <- traits_surv[common_samples, ]
expr_final <- matriz_vst[, common_samples]

cat("[OK] Muestras con supervivencia valida:", length(common_samples), "\n")
cat("  Mediana supervivencia:", round(median(traits_final$survival, na.rm = TRUE), 1), "meses\n")
cat("  Eventos (muertes):", sum(traits_final$death), "de", nrow(traits_final), "\n\n")

# --- seleccionar genes para analisis ---
# score combinado: |kME| x |GS| x |log2FC|
key_genes$score <- abs(key_genes$kME) * abs(key_genes$GS_death) * abs(key_genes$log2FoldChange)
key_genes <- key_genes[order(-key_genes$score), ]

# analizar todos los genes clave con nombre conocido
genes_to_analyze <- key_genes[!is.na(key_genes$gene_symbol) &
                                key_genes$gene_symbol != "NA", ]

# verificar disponibilidad en la matriz
genes_available <- genes_to_analyze$ensembl_id[genes_to_analyze$ensembl_id %in% rownames(expr_final)]
genes_to_analyze <- genes_to_analyze[genes_to_analyze$ensembl_id %in% genes_available, ]

cat("[OK] Genes a analizar:", nrow(genes_to_analyze), "\n")
for (i in 1:nrow(genes_to_analyze)) {
  cat("  ", genes_to_analyze$gene_symbol[i], "(", genes_to_analyze$module_color[i], ")\n")
}
cat("\n")

# --- Kaplan-Meier por gen individual ---
km_results <- list()

for (i in 1:nrow(genes_to_analyze)) {

  gene_ensembl <- genes_to_analyze$ensembl_id[i]
  gene_symbol <- genes_to_analyze$gene_symbol[i]
  gene_module <- genes_to_analyze$module_color[i]

  cat("  Analizando", gene_symbol, "...\n")

  gene_expr <- as.numeric(expr_final[gene_ensembl, common_samples])

  # dicotomizo por mediana
  median_expr <- median(gene_expr, na.rm = TRUE)
  group <- ifelse(gene_expr >= median_expr, "Alta expresión", "Baja expresión")

  surv_data <- data.frame(
    time = traits_final$survival,
    event = traits_final$death,
    expression_group = factor(group, levels = c("Baja expresión", "Alta expresión")),
    expression = gene_expr,
    age = traits_final$age,
    stage = traits_final$stage,
    stringsAsFactors = FALSE
  )

  surv_data <- surv_data[complete.cases(surv_data[, c("time", "event", "expression_group")]), ]

  surv_obj <- Surv(surv_data$time, surv_data$event)
  km_fit <- survfit(surv_obj ~ expression_group, data = surv_data)

  # log-rank test
  logrank <- survdiff(surv_obj ~ expression_group, data = surv_data)
  logrank_p <- 1 - pchisq(logrank$chisq, df = 1)

  median_surv <- surv_median(km_fit)

  cat("    Log-rank p:", format(logrank_p, digits = 3),
      "| Mediana alta:", round(median_surv$median[2], 1),
      "| baja:", round(median_surv$median[1], 1), "\n")

  km_results[[gene_symbol]] <- list(
    gene = gene_symbol,
    ensembl = gene_ensembl,
    module = gene_module,
    logrank_p = logrank_p,
    median_high = median_surv$median[2],
    median_low = median_surv$median[1],
    n_high = sum(surv_data$expression_group == "Alta expresión"),
    n_low = sum(surv_data$expression_group == "Baja expresión")
  )

  # solo guardar grafico KM para genes significativos (evitar 200+ graficos)
  if (logrank_p < 0.05) {
    color_palette <- if (gene_module == "yellow") c("#4ECDC4", "#E6AB02") else c("#4ECDC4", "#E7298A")

    p <- ggsurvplot(
      km_fit,
      data = surv_data,
      pval = TRUE,
      pval.coord = c(0, 0.1),
      conf.int = TRUE,
      risk.table = TRUE,
      risk.table.col = "strata",
      risk.table.height = 0.25,
      palette = color_palette,
      xlab = "Tiempo (meses)",
      ylab = "Probabilidad de supervivencia",
      title = paste0("Kaplan-Meier: ", gene_symbol, " (", gene_module, ")\n",
                     "Dicotomizado por mediana de expresión"),
      legend.title = gene_symbol,
      legend.labs = c("Baja expresión", "Alta expresión"),
      surv.median.line = "hv",
      ggtheme = theme_classic() + theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "right"
      )
    )

    filename <- sprintf("%s/KM_%02d_%s.png", fig_dir, i, gene_symbol)
    ggsave(filename, plot = print(p), width = 10, height = 8, dpi = 300)
    cat("    [OK] Significativo - guardado:", basename(filename), "\n\n")
  } else {
    cat("    No significativo - sin grafico\n\n")
  }
}

# --- tabla resumen Kaplan-Meier ---
km_summary <- do.call(rbind, lapply(km_results, function(x) {
  data.frame(
    Gene = x$gene,
    Module = x$module,
    LogRank_p = x$logrank_p,
    Median_High = x$median_high,
    Median_Low = x$median_low,
    n_High = x$n_high,
    n_Low = x$n_low,
    Significant = ifelse(x$logrank_p < 0.05, "SI", "no"),
    stringsAsFactors = FALSE
  )
}))

km_summary <- km_summary[order(km_summary$LogRank_p), ]

# correccion por comparaciones multiples (Benjamini-Hochberg)
km_summary$FDR <- p.adjust(km_summary$LogRank_p, method = "BH")
km_summary$Significant_FDR <- ifelse(km_summary$FDR < 0.05, "SI", "no")

print(km_summary[, c("Gene", "Module", "LogRank_p", "FDR", "Significant_FDR")])

n_sig <- sum(km_summary$LogRank_p < 0.05)
n_fdr <- sum(km_summary$FDR < 0.05)
cat("\nGenes significativos (p < 0.05):", n_sig, "de", nrow(km_summary), "\n")
cat("Genes significativos tras FDR < 0.05:", n_fdr, "\n")
cat("Genes con FDR < 0.01:", sum(km_summary$FDR < 0.01), "\n")

write.csv(km_summary,
          file.path(out_dir, "kaplan_meier_summary.csv"),
          row.names = FALSE)
write.csv(km_summary,
          file.path(res_dir, "tabla_kaplan_meier_summary.csv"),
          row.names = FALSE)
cat("[OK] Guardado: kaplan_meier_summary.csv\n\n")

# --- Cox univariante ---
cox_results <- data.frame()

for (i in 1:nrow(genes_to_analyze)) {

  gene_ensembl <- genes_to_analyze$ensembl_id[i]
  gene_symbol <- genes_to_analyze$gene_symbol[i]

  gene_expr <- as.numeric(expr_final[gene_ensembl, common_samples])

  surv_data <- data.frame(
    time = traits_final$survival,
    event = traits_final$death,
    expression = scale(gene_expr)[, 1],
    stringsAsFactors = FALSE
  )
  surv_data <- surv_data[complete.cases(surv_data), ]

  cox_fit <- coxph(Surv(time, event) ~ expression, data = surv_data)
  cox_sum <- summary(cox_fit)

  hr <- cox_sum$conf.int[1, 1]
  hr_lower <- cox_sum$conf.int[1, 3]
  hr_upper <- cox_sum$conf.int[1, 4]
  p_value <- cox_sum$coefficients[1, 5]
  concordance <- cox_sum$concordance[1]

  cox_results <- rbind(cox_results, data.frame(
    Gene = gene_symbol,
    Ensembl = gene_ensembl,
    Module = genes_to_analyze$module_color[i],
    HR = round(hr, 3),
    HR_CI_lower = round(hr_lower, 3),
    HR_CI_upper = round(hr_upper, 3),
    Cox_pvalue = p_value,
    Concordance = round(concordance, 3),
    Risk = ifelse(hr > 1, "Riesgo", "Protector"),
    Significant = ifelse(p_value < 0.05, "SI", "no"),
    stringsAsFactors = FALSE
  ))

  cat("  ", gene_symbol, ": HR =", round(hr, 3),
      paste0("(", round(hr_lower, 2), "-", round(hr_upper, 2), ")"),
      "p =", format(p_value, digits = 3),
      ifelse(p_value < 0.05, "***", ""),
      ifelse(hr > 1, "Riesgo", "Protector"), "\n")
}

cox_results <- cox_results[order(cox_results$Cox_pvalue), ]

# correccion FDR para Cox
cox_results$FDR <- p.adjust(cox_results$Cox_pvalue, method = "BH")
cox_results$Significant_FDR <- ifelse(cox_results$FDR < 0.05, "SI", "no")

write.csv(cox_results,
          file.path(out_dir, "cox_univariate_results.csv"),
          row.names = FALSE)
write.csv(cox_results,
          file.path(res_dir, "tabla_cox_univariate_results.csv"),
          row.names = FALSE)
cat("\n[OK] Guardado: cox_univariate_results.csv\n\n")

# --- forest plot (top 10 yellow + top 5 red por FDR) ---
cox_named <- cox_results[cox_results$Gene != "" & cox_results$FDR < 0.05, ]
cox_named <- cox_named[!duplicated(cox_named$Gene), ]
cox_yellow <- head(cox_named[cox_named$Module == "yellow", ], 10)
cox_red <- head(cox_named[cox_named$Module == "red", ], 5)
cox_sig <- rbind(cox_yellow, cox_red)
cox_sig <- cox_sig[order(cox_sig$FDR), ]

if (nrow(cox_sig) > 0) {
  cox_sig$Gene <- factor(cox_sig$Gene, levels = rev(cox_sig$Gene))
  cox_sig$color <- ifelse(cox_sig$Module == "yellow", "#E6AB02", "#E7298A")

  p_forest <- ggplot(cox_sig, aes(x = HR, y = Gene)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
    geom_errorbarh(aes(xmin = HR_CI_lower, xmax = HR_CI_upper),
                   height = 0.25, linewidth = 0.8, color = cox_sig$color) +
    geom_point(size = 4, shape = 18, color = cox_sig$color) +
    geom_text(aes(label = sprintf("HR=%.2f (p=%s)", HR,
                                   ifelse(Cox_pvalue < 0.001, "<0.001",
                                          format(round(Cox_pvalue, 3), nsmall = 3)))),
              hjust = -0.15, vjust = -0.7, size = 3.5) +
    labs(
      title = "Forest Plot: Hazard Ratios (Cox Univariante)\nTop 10 yellow (riesgo) + Top 5 red (protector)",
      x = "Hazard Ratio (IC 95%)",
      y = ""
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.y = element_text(face = "bold", size = 11),
      axis.text.x = element_text(size = 10)
    ) +
    annotate("text", x = max(cox_sig$HR_CI_upper) * 0.9, y = 0.5,
             label = "Yellow = ME4 (muerte)\nRed = ME6 (muerte)",
             size = 3.5, hjust = 1, fontface = "italic")

  plot_height <- max(6, nrow(cox_sig) * 0.5)
  ggsave(file.path(fig_dir, "forest_plot_cox.png"), p_forest,
         width = 12, height = plot_height, dpi = 300, limitsize = FALSE)
  ggsave(file.path(res_dir, "fig_19_forest_plot_cox.png"), p_forest,
         width = 12, height = plot_height, dpi = 300, limitsize = FALSE)
  cat("[OK] Guardado: forest_plot_cox.png\n\n")
} else {
  cat("[AVISO] Ningun gen significativo en Cox para forest plot\n\n")
}

# --- Cox multivariante (top 10 genes con FDR < 0.05) ---
sig_genes <- cox_results[cox_results$FDR < 0.05 & cox_results$Gene != "", ]
sig_genes <- head(sig_genes, 10)  # limitar a top 10 para evitar sobreajuste

if (nrow(sig_genes) >= 2) {
  cat("Genes para modelo multivariante (top 10):", nrow(sig_genes), "\n")
  cat("Genes:", paste(sig_genes$Gene, collapse = ", "), "\n\n")

  multi_data <- data.frame(
    time = traits_final$survival,
    event = traits_final$death,
    age = traits_final$age,
    stage = traits_final$stage
  )

  # usar nombres seguros para R (reemplazar caracteres problematicos)
  safe_names <- make.names(sig_genes$Gene)
  for (j in 1:nrow(sig_genes)) {
    gene_ensembl <- sig_genes$Ensembl[j]
    gene_expr <- as.numeric(expr_final[gene_ensembl, common_samples])
    multi_data[[safe_names[j]]] <- scale(gene_expr)[, 1]
  }

  multi_data <- multi_data[complete.cases(multi_data), ]

  gene_vars <- paste(safe_names, collapse = " + ")
  formula_multi <- as.formula(paste("Surv(time, event) ~ age + stage +", gene_vars))

  cat("Formula:", deparse(formula_multi), "\n\n")

  cox_multi <- coxph(formula_multi, data = multi_data)
  print(summary(cox_multi))

  multi_summary <- as.data.frame(summary(cox_multi)$coefficients)
  multi_summary$Variable <- rownames(multi_summary)
  multi_summary$HR <- exp(multi_summary$coef)
  multi_ci <- as.data.frame(summary(cox_multi)$conf.int)
  multi_summary$HR_lower <- multi_ci[, 3]
  multi_summary$HR_upper <- multi_ci[, 4]

  write.csv(multi_summary,
            file.path(out_dir, "cox_multivariate_results.csv"),
            row.names = FALSE)
  cat("\n[OK] Guardado: cox_multivariate_results.csv\n")

} else if (nrow(sig_genes) == 1) {
  cat("Solo 1 gen significativo. Modelo multivariante con covariables clinicas:\n\n")

  gene_ensembl <- sig_genes$Ensembl[1]
  gene_symbol <- sig_genes$Gene[1]
  gene_expr <- as.numeric(expr_final[gene_ensembl, common_samples])

  multi_data <- data.frame(
    time = traits_final$survival,
    event = traits_final$death,
    age = traits_final$age,
    stage = traits_final$stage,
    gene = scale(gene_expr)[, 1]
  )
  colnames(multi_data)[5] <- gene_symbol
  multi_data <- multi_data[complete.cases(multi_data), ]

  formula_multi <- as.formula(paste("Surv(time, event) ~ age + stage +", gene_symbol))
  cox_multi <- coxph(formula_multi, data = multi_data)
  print(summary(cox_multi))

  write.csv(
    as.data.frame(summary(cox_multi)$coefficients),
    file.path(out_dir, "cox_multivariate_results.csv")
  )
  cat("\n[OK] Guardado: cox_multivariate_results.csv\n")

} else {
  cat("[AVISO] Ningun gen significativo en Cox univariante (p < 0.05).\n")

  trend_genes <- cox_results[cox_results$Cox_pvalue < 0.1, ]
  if (nrow(trend_genes) > 0) {
    cat("  Genes con tendencia (p < 0.1):", paste(trend_genes$Gene, collapse = ", "), "\n")
  }
}

# --- resumen ---
cat("\nResumen supervivencia:\n")
cat("  Muestras:", length(common_samples), "| Genes:", nrow(genes_to_analyze), "\n")
cat("  KM significativos:", sum(km_summary$LogRank_p < 0.05), "\n")
cat("  Cox significativos:", sum(cox_results$Cox_pvalue < 0.05), "\n")
cat("  Archivos en:", out_dir, "\n")
