# --- 08_wgcna_correlacion_traits.R ---
# Autor: Ines Torres
# Fecha: 10/02/2026
# Correlaciono modulos WGCNA con traits clinicos

library(WGCNA)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

allowWGCNAThreads()

setwd("~/TFM_PDAC")

if (!dir.exists("04_red_coexpresion/WGCNA/traits")) {
  dir.create("04_red_coexpresion/WGCNA/traits", recursive = TRUE)
}
if (!dir.exists("resultados_TFM_v2")) {
  dir.create("resultados_TFM_v2", recursive = TRUE)
}

# --- cargar datos ---
MEs_colors <- read.table(
  "04_red_coexpresion/WGCNA/outputs/module_eigengenes.txt",
  header = TRUE, row.names = 1, sep = "\t"
)

module_assignment <- read.table(
  "04_red_coexpresion/WGCNA/outputs/module_assignment.txt",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)

geneModuleMembership <- read.table(
  "04_red_coexpresion/WGCNA/outputs/gene_module_membership.txt",
  header = TRUE, row.names = 1, sep = "\t"
)

traits_raw <- read.table(
  "01_datos/clinical/TCGA_PAAD_traits_numeric.txt",
  header = TRUE, row.names = 1, sep = "\t"
)

# alinear muestras
rownames(traits_raw) <- gsub("-", ".", rownames(traits_raw))
muestras_comunes <- intersect(rownames(MEs_colors), rownames(traits_raw))

MEs_final <- MEs_colors[muestras_comunes, ]
datTraits <- traits_raw[muestras_comunes, ]

cat("[OK] Datos cargados. Muestras:", nrow(MEs_final),
    "| Modulos:", ncol(MEs_final), "| Traits:", ncol(datTraits), "\n\n")

# --- correlaciones modulo-trait ---
moduleTraitCor <- cor(MEs_final, datTraits, use = "pairwise.complete.obs")
moduleTraitPvalue_raw <- corPvalueStudent(moduleTraitCor, nrow(MEs_final))

# correccion BH: testamos N modulos x 5 traits, hay que corregir por multiples comparaciones
pvals_vector <- as.vector(moduleTraitPvalue_raw)
pvals_adj <- p.adjust(pvals_vector, method = "BH")
moduleTraitPvalue <- matrix(pvals_adj,
                            nrow = nrow(moduleTraitPvalue_raw),
                            ncol = ncol(moduleTraitPvalue_raw),
                            dimnames = dimnames(moduleTraitPvalue_raw))

cat("[OK] P-valores corregidos por Benjamini-Hochberg (",
    length(pvals_vector), "tests)\n\n")

# --- heatmap modulo-trait ---
textMatrix <- paste0(
  round(moduleTraitCor, 2), "\n(",
  round(moduleTraitPvalue, 3), ")"
)
dim(textMatrix) <- dim(moduleTraitCor)

trait_labels <- c("Edad", "Sexo\n(M=1)", "Estadio\n(I-IV)", "Superv.\n(meses)", "Muerte\n(0/1)")

png(
  "resultados_TFM_v2/fig_10_module_trait_heatmap.png",
  width = 10, height = 12, units = "in", res = 300
)

par(mar = c(6, 8.5, 3, 3))

labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = trait_labels,
  yLabels = rownames(moduleTraitCor),
  ySymbols = rownames(moduleTraitCor),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.5,
  zlim = c(-1, 1),
  main = "Correlacion Modulo-Trait\nTCGA-PAAD (n=165)"
)

dev.off()

cat("[OK] Heatmap guardado\n\n")

# --- modulos significativos (tres niveles de criterio) ---
get_sig_modules <- function(cor_matrix, pval_matrix, cor_thresh, pval_thresh) {
  results <- list()

  for (trait in colnames(cor_matrix)) {
    cors <- cor_matrix[, trait]
    pvals <- pval_matrix[, trait]
    sig_idx <- which(abs(cors) >= cor_thresh & pvals < pval_thresh)

    if (length(sig_idx) > 0) {
      df <- data.frame(
        Module = names(cors)[sig_idx],
        Correlation = round(cors[sig_idx], 3),
        Pvalue = round(pvals[sig_idx], 4),
        Direction = ifelse(cors[sig_idx] > 0, "Positiva", "Negativa"),
        Trait = trait,
        stringsAsFactors = FALSE
      )
      df <- df[order(-abs(df$Correlation)), ]
      results[[trait]] <- df
    }
  }
  return(results)
}

# funcion auxiliar para contar asociaciones (maneja listas vacias)
count_sig <- function(sig_list) {
  if (length(sig_list) == 0) return(0)
  sum(sapply(sig_list, nrow))
}

sig_nivel1 <- get_sig_modules(moduleTraitCor, moduleTraitPvalue, 0.30, 0.05)
n_nivel1 <- count_sig(sig_nivel1)

sig_nivel2 <- get_sig_modules(moduleTraitCor, moduleTraitPvalue, 0.20, 0.05)
n_nivel2 <- count_sig(sig_nivel2)

sig_nivel3 <- get_sig_modules(moduleTraitCor, moduleTraitPvalue, 0.15, 0.10)
n_nivel3 <- count_sig(sig_nivel3)

cat("Asociaciones significativas:\n")
cat("  Nivel 1 (|r|>0.30, p<0.05):", n_nivel1, "\n")
cat("  Nivel 2 (|r|>0.20, p<0.05):", n_nivel2, "\n")
cat("  Nivel 3 (|r|>0.15, p<0.10):", n_nivel3, "\n\n")

# uso nivel 2 para analisis principal
sig_modules_main <- sig_nivel2

for (trait in names(sig_modules_main)) {
  cat("TRAIT:", trait, "\n")
  print(sig_modules_main[[trait]])
  cat("\n")
}

if (length(sig_modules_main) > 0) {
  sig_modules_df <- do.call(rbind, sig_modules_main)
  rownames(sig_modules_df) <- NULL

  write.table(
    sig_modules_df,
    "04_red_coexpresion/WGCNA/traits/significant_modules_all_levels.txt",
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  write.table(
    sig_modules_df,
    "resultados_TFM_v2/tabla_significant_modules_all_levels.txt",
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  cat("[OK] Guardado: significant_modules_all_levels.txt\n\n")
}

# --- tabla completa de correlaciones ---
cor_table <- as.data.frame(moduleTraitCor)
cor_table$Module <- rownames(cor_table)

pval_table <- as.data.frame(moduleTraitPvalue)

sig_table <- cor_table
for (trait in colnames(datTraits)) {
  for (mod in rownames(moduleTraitCor)) {
    r <- round(moduleTraitCor[mod, trait], 3)
    p <- moduleTraitPvalue[mod, trait]
    stars <- ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
    sig_table[mod, trait] <- paste0(r, stars)
  }
}

write.table(
  sig_table,
  "04_red_coexpresion/WGCNA/traits/module_trait_correlations_full.txt",
  sep = "\t", quote = FALSE, row.names = TRUE
)
write.table(
  sig_table,
  "resultados_TFM_v2/tabla_module_trait_correlations_full.txt",
  sep = "\t", quote = FALSE, row.names = TRUE
)

cat("[OK] Guardado: module_trait_correlations_full.txt (* p<0.05, ** p<0.01, *** p<0.001)\n\n")

# --- genes hub con criterios ajustados ---
# cargo matriz de expresion
matriz_vst <- read.table(
  "02_preprocesamiento/outputs_FINAL/matriz_vst_normalizada.txt",
  header = TRUE, row.names = 1, sep = "\t"
)

datExpr <- t(matriz_vst)
datExpr <- datExpr[muestras_comunes, ]
genes_en_red <- module_assignment$Gene
genes_disponibles <- intersect(genes_en_red, colnames(datExpr))
datExpr_filtrado <- datExpr[, genes_disponibles]

cat("[OK] Expresion cargada:", nrow(datExpr_filtrado), "muestras,",
    ncol(datExpr_filtrado), "genes\n\n")

# calculo GS para todos los traits
geneSignificance <- list()

for (trait in colnames(datTraits)) {
  trait_values <- datTraits[, trait]
  valid_idx <- !is.na(trait_values)

  gs_cors <- cor(
    datExpr_filtrado[valid_idx, ],
    trait_values[valid_idx],
    use = "pairwise.complete.obs"
  )
  geneSignificance[[trait]] <- gs_cors
}

# --- construir mapa ME_numero -> color ---
# module_assignment usa colores ("turquoise","blue",...) pero los eigengenes
# usan etiquetas numericas ("ME1","ME2",...). Necesito la correspondencia.
# WGCNA asigna: 0=grey, 1=modulo mas grande, 2=segundo mas grande, etc.
module_sizes <- sort(table(module_assignment$Module), decreasing = TRUE)
module_sizes <- module_sizes[names(module_sizes) != "grey"]  # grey siempre es 0

# mapa: ME1 -> color mas grande, ME2 -> segundo, etc.
me_to_color <- c("0" = "grey")
for (idx in seq_along(module_sizes)) {
  me_to_color[as.character(idx)] <- names(module_sizes)[idx]
}

cat("[OK] Mapa ME -> color:\n")
for (me_num in names(me_to_color)) {
  cat("  ME", me_num, " -> ", me_to_color[me_num], "\n", sep = "")
}
cat("\n")

# --- identifico hubs: kME > 0.5, GS > 0.15 ---
if (length(sig_modules_main) > 0) {
  modules_to_analyze <- sig_modules_df
} else if (length(sig_nivel3) > 0) {
  cat("[AVISO] Sin modulos en nivel 2, usando nivel 3 (criterios relajados)\n")
  sig_modules_main <- sig_nivel3
  modules_to_analyze <- do.call(rbind, sig_nivel3)
} else {
  cat("[AVISO] Sin modulos significativos en ningun nivel.\n")
  cat("  Usando todos los modulos no-grey para buscar hub genes.\n")
  all_modules <- rownames(moduleTraitCor)
  modules_to_analyze <- data.frame(
    Module = all_modules,
    Correlation = moduleTraitCor[, "death"],
    Pvalue = moduleTraitPvalue[, "death"],
    Direction = ifelse(moduleTraitCor[, "death"] > 0, "Positiva", "Negativa"),
    Trait = "death",
    stringsAsFactors = FALSE
  )
  modules_to_analyze <- modules_to_analyze[modules_to_analyze$Module != "ME0", ]
}

hub_genes_final <- data.frame()

if (!is.null(modules_to_analyze) && nrow(modules_to_analyze) > 0) {

  for (i in 1:nrow(modules_to_analyze)) {

    trait_name <- modules_to_analyze$Trait[i]
    mod_name <- modules_to_analyze$Module[i]
    me_num <- gsub("ME", "", mod_name)

    # convertir numero a color usando el mapa
    if (me_num %in% names(me_to_color)) {
      mod_color <- me_to_color[me_num]
    } else {
      cat("  [WARN] No se encontro color para", mod_name, "- saltando\n")
      next
    }

    genes_modulo <- module_assignment$Gene[module_assignment$Module == mod_color]
    genes_modulo <- intersect(genes_modulo, colnames(datExpr_filtrado))

    if (length(genes_modulo) == 0) next

    kME_col <- mod_name
    if (!kME_col %in% colnames(geneModuleMembership)) next
    kME_vals <- geneModuleMembership[genes_modulo, kME_col]

    gs_vals <- geneSignificance[[trait_name]][genes_modulo, 1]

    is_hub <- abs(kME_vals) > 0.5 & abs(gs_vals) > 0.15

    if (sum(is_hub, na.rm = TRUE) > 0) {
      hub_df <- data.frame(
        Gene = genes_modulo[is_hub],
        Module = mod_color,
        ME = mod_name,
        Trait = trait_name,
        kME = round(kME_vals[is_hub], 3),
        GS = round(gs_vals[is_hub], 3),
        Trait_correlation = modules_to_analyze$Correlation[i],
        stringsAsFactors = FALSE
      )
      hub_df <- hub_df[order(-abs(hub_df$kME)), ]
      hub_genes_final <- rbind(hub_genes_final, hub_df)
    }
  }
}

if (nrow(hub_genes_final) > 0) {
  cat("[OK] Genes HUB identificados:", nrow(hub_genes_final), "\n")
  cat("Top 30:\n")
  print(head(hub_genes_final[order(-abs(hub_genes_final$kME)), ], 30))
  cat("\n")

  write.table(
    hub_genes_final,
    "04_red_coexpresion/WGCNA/traits/hub_genes_with_traits.txt",
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  write.table(
    hub_genes_final,
    "resultados_TFM_v2/tabla_hub_genes_with_traits.txt",
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  cat("[OK] Guardado: hub_genes_with_traits.txt\n\n")

} else {
  # criterios minimos: kME > 0.4, GS > 0.1
  cat("[AVISO] Sin hubs con kME>0.5, GS>0.15. Usando criterios minimos.\n\n")

  hub_genes_final <- data.frame()

  for (i in 1:nrow(modules_to_analyze)) {
    trait_name <- modules_to_analyze$Trait[i]
    mod_name <- modules_to_analyze$Module[i]
    me_num <- gsub("ME", "", mod_name)

    if (me_num %in% names(me_to_color)) {
      mod_color <- me_to_color[me_num]
    } else {
      next
    }

    genes_modulo <- module_assignment$Gene[module_assignment$Module == mod_color]
    genes_modulo <- intersect(genes_modulo, colnames(datExpr_filtrado))
    if (length(genes_modulo) == 0) next

    kME_col <- mod_name
    if (!kME_col %in% colnames(geneModuleMembership)) next
    kME_vals <- geneModuleMembership[genes_modulo, kME_col]
    gs_vals <- geneSignificance[[trait_name]][genes_modulo, 1]

    is_hub <- abs(kME_vals) > 0.4 & abs(gs_vals) > 0.1

    if (sum(is_hub, na.rm = TRUE) > 0) {
      hub_df <- data.frame(
        Gene = genes_modulo[is_hub],
        Module = mod_color,
        ME = mod_name,
        Trait = trait_name,
        kME = round(kME_vals[is_hub], 3),
        GS = round(gs_vals[is_hub], 3),
        stringsAsFactors = FALSE
      )
      hub_df <- hub_df[order(-abs(hub_df$kME)), ]
      hub_genes_final <- rbind(hub_genes_final, hub_df)
    }
  }

  if (nrow(hub_genes_final) > 0) {
    cat("[OK] Genes HUB (criterios minimos):", nrow(hub_genes_final), "\n")
    print(head(hub_genes_final[order(-abs(hub_genes_final$kME)), ], 20))

    write.table(
      hub_genes_final,
      "04_red_coexpresion/WGCNA/traits/hub_genes_with_traits.txt",
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    write.table(
      hub_genes_final,
      "resultados_TFM_v2/tabla_hub_genes_with_traits.txt",
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    cat("\n[OK] Guardado: hub_genes_with_traits.txt\n\n")
  }
}

# --- interseccion con DEGs ---
if (file.exists("03_expresion_diferencial/outputs/DEGs_UP_regulated.txt") &&
    nrow(hub_genes_final) > 0) {

  degs_up <- read.table(
    "03_expresion_diferencial/outputs/DEGs_UP_regulated.txt",
    header = TRUE, sep = "\t", stringsAsFactors = FALSE
  )

  cat("DEGs UP:", nrow(degs_up), "| Hub genes:", nrow(hub_genes_final), "\n")

  final_candidates <- hub_genes_final[hub_genes_final$Gene %in% degs_up$Gene, ]

  if (nrow(final_candidates) > 0) {

    final_candidates <- merge(
      final_candidates,
      degs_up[, c("Gene", "log2FoldChange", "padj")],
      by = "Gene"
    )

    final_candidates <- final_candidates[order(-abs(final_candidates$kME)), ]

    cat("[OK] Genes candidatos finales:", nrow(final_candidates), "\n")
    print(head(final_candidates, 20))
    cat("\n")

    write.table(
      final_candidates,
      "04_red_coexpresion/WGCNA/traits/GENES_CANDIDATOS_FINALES.txt",
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    cat("[OK] Guardado: GENES_CANDIDATOS_FINALES.txt\n\n")

  } else {
    cat("[AVISO] No hay interseccion hub-DEG\n\n")
  }
}

# --- validacion genes conocidos ---
genes_pancreas <- c("KRAS", "TP53", "CDKN2A", "SMAD4",
                    "MUC1", "MUC4", "MSLN", "S100A4",
                    "CEACAM5", "CEACAM6", "KRT19", "VIM",
                    "S100P", "ANXA1", "ANXA2")

genes_en_modulos <- module_assignment[module_assignment$Gene %in% genes_pancreas, ]

if (nrow(genes_en_modulos) > 0) {
  cat("Genes conocidos en modulos:\n")
  print(genes_en_modulos)
  cat("\n")
}

if (nrow(hub_genes_final) > 0) {
  hub_conocidos <- hub_genes_final[hub_genes_final$Gene %in% genes_pancreas, ]

  if (nrow(hub_conocidos) > 0) {
    cat("Genes conocidos que son HUB:\n")
    print(hub_conocidos)
    cat("\n")
  }
}

# --- resumen ---
cat("Resumen:\n")
cat("  Modulos:", ncol(MEs_final), "| Muestras:", nrow(MEs_final), "\n")
cat("  Sig nivel1:", n_nivel1, "| nivel2:", n_nivel2, "| nivel3:", n_nivel3, "\n")
if (nrow(hub_genes_final) > 0) cat("  Hub genes:", nrow(hub_genes_final), "\n")
if (exists("final_candidates") && nrow(final_candidates) > 0) {
  cat("  Candidatos finales:", nrow(final_candidates), "\n")
}
