# --- 07_wgcna_red_coexpresion.R ---
# Autor: Ines Torres
# Fecha: 07/02/2026
# Construyo la red de coexpresion WGCNA e identifico modulos

# --- librerias ---
library(WGCNA)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(dplyr)

allowWGCNAThreads()

# --- configurar directorio ---
setwd("~/TFM_PDAC")

if (!dir.exists("04_red_coexpresion/WGCNA/outputs")) {
  dir.create("04_red_coexpresion/WGCNA/outputs", recursive = TRUE)
}
if (!dir.exists("resultados_TFM_v2")) {
  dir.create("resultados_TFM_v2", recursive = TRUE)
}

# --- cargar datos ---
matriz_vst <- read.table(
  "02_preprocesamiento/outputs_FINAL/matriz_vst_normalizada.txt",
  header = TRUE,
  row.names = 1,
  sep = "\t"
)

metadata <- read.table(
  "02_preprocesamiento/outputs_FINAL/metadata_muestras.txt",
  header = TRUE,
  row.names = 1,
  sep = "\t"
)

cat("[OK] Datos cargados\n")
cat("  Genes:", nrow(matriz_vst), "| Muestras matriz:", ncol(matriz_vst),
    "| Muestras metadata:", nrow(metadata), "\n\n")

# --- sincronizar nombres de muestras (guiones a puntos) ---
rownames(metadata) <- gsub("-", ".", rownames(metadata))
metadata$sample_id <- gsub("-", ".", metadata$sample_id)

coincidencias <- sum(rownames(metadata) %in% colnames(matriz_vst))
cat("[OK] Nombres convertidos. Coincidencias:", coincidencias, "de", nrow(metadata), "\n\n")

if (coincidencias == 0) {
  stop("ERROR: No coinciden los nombres. Revisar manualmente.")
}

# --- seleccionar muestras tumorales ---
muestras_tumor <- rownames(metadata[metadata$grupo == "TCGA_Tumor", ])
muestras_tumor_disponibles <- muestras_tumor[muestras_tumor %in% colnames(matriz_vst)]

if (length(muestras_tumor_disponibles) == 0) {
  stop("ERROR: No hay muestras tumorales en la matriz")
}

# WGCNA necesita muestras en filas, genes en columnas
datExpr <- t(matriz_vst[, muestras_tumor_disponibles])

cat("[OK] Matriz expresion: ", nrow(datExpr), "muestras x", ncol(datExpr), "genes\n\n")

# --- control de calidad - outliers ---
sampleTree <- hclust(dist(datExpr), method = "average")

# umbral automatico con IQR
alturas <- sampleTree$height
Q1 <- quantile(alturas, 0.25)
Q3 <- quantile(alturas, 0.75)
IQR_val <- Q3 - Q1
umbral_outlier <- Q3 + 1.5 * IQR_val

png("resultados_TFM_v2/fig_06_sample_clustering.png",
    width = 12, height = 6, units = "in", res = 300)

par(cex = 0.6)
par(mar = c(4, 5, 2, 1))
plot(sampleTree,
     main = "Clustering de Muestras - Detección de Outliers",
     sub = "",
     xlab = "",
     cex.lab = 1.5,
     cex.axis = 1.5,
     cex.main = 2)

abline(h = umbral_outlier, col = "red", lwd = 2, lty = 2)
text(50, umbral_outlier * 1.1,
     paste("Umbral:", round(umbral_outlier)),
     col = "red", cex = 1.2)

dev.off()

# identificar y eliminar outliers
clust <- cutreeStatic(sampleTree, cutHeight = umbral_outlier, minSize = 10)
keepSamples <- (clust == 1)
datExpr_clean <- datExpr[keepSamples, ]
n_outliers <- sum(!keepSamples)

cat("[OK] QC: ", nrow(datExpr_clean), "muestras mantenidas,", n_outliers, "outliers eliminados\n")

if (n_outliers > 0) {
  cat("  Muestras eliminadas:\n")
  print(rownames(datExpr)[!keepSamples])
}
cat("\n")

datExpr <- datExpr_clean

# --- seleccion de genes: DEGs estrictos + filtrado por MAD ---
# primero filtro los DEGs con criterios mas estrictos (padj < 0.01, |log2FC| > 2)
# y luego me quedo solo con los que tienen alta variabilidad entre pacientes (MAD > P75)
# asi reduzco los ~26000 DEGs a un numero manejable para WGCNA

degs_all <- read.table(
  "03_expresion_diferencial/outputs/resultados_completos_DESeq2.txt",
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

# PASO 1: filtrado estricto por significancia y magnitud de cambio
degs_strict <- degs_all[!is.na(degs_all$padj) &
                          degs_all$padj < 0.01 &
                          abs(degs_all$log2FoldChange) > 2, ]

cat("[OK] DEGs estrictos (padj<0.01, |log2FC|>2):", nrow(degs_strict), "\n")
cat("  UP:", sum(degs_strict$log2FoldChange > 0),
    "/ DOWN:", sum(degs_strict$log2FoldChange < 0), "\n")

# filtrar la matriz de expresion a estos DEGs
genes_degs <- degs_strict$gene_id
genes_disponibles <- intersect(genes_degs, colnames(datExpr))
cat("  Presentes en matriz tumoral:", length(genes_disponibles), "\n")

datExpr_degs <- datExpr[, genes_disponibles]

# PASO 2: filtrado por MAD (Median Absolute Deviation)
# MAD mide la variabilidad de cada gen entre pacientes tumorales
# retenemos genes por encima del percentil 75 de MAD
gene_mad <- apply(datExpr_degs, 2, mad)
mad_threshold <- quantile(gene_mad, probs = 0.75)
genes_high_mad <- names(gene_mad[gene_mad >= mad_threshold])

datExpr_filtered <- datExpr_degs[, genes_high_mad]

cat("\n[OK] Filtrado por MAD (percentil 75):\n")
cat("  Umbral MAD:", round(mad_threshold, 3), "\n")
cat("  Genes retenidos:", ncol(datExpr_filtered), "de", length(genes_disponibles), "\n")
cat("  (criterio: MAD >= P75, no es un numero fijo arbitrario)\n\n")

# verificar calidad
gsg <- goodSamplesGenes(datExpr_filtered, verbose = 3)

if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0) cat("  Genes problematicos:", sum(!gsg$goodGenes), "\n")
  if (sum(!gsg$goodSamples) > 0) cat("  Muestras problematicas:", sum(!gsg$goodSamples), "\n")
  datExpr_filtered <- datExpr_filtered[gsg$goodSamples, gsg$goodGenes]
}

cat("[OK] Datos filtrados:", nrow(datExpr_filtered), "muestras,", ncol(datExpr_filtered), "genes\n\n")

# --- seleccion de soft threshold (beta) ---
powers <- c(seq(1, 10, by = 1), seq(12, 20, by = 2))

sft <- pickSoftThreshold(
  datExpr_filtered,
  powerVector = powers,
  verbose = 5,
  networkType = "signed"
)

# grafico de seleccion
png("resultados_TFM_v2/fig_07_soft_threshold_selection.png",
    width = 12, height = 6, units = "in", res = 300)

par(mfrow = c(1, 2))

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n",
     main = "Scale independence")

text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers,
     cex = 0.9,
     col = "red")

abline(h = 0.80, col = "blue", lty = 2)
text(5, 0.85, "R^2 = 0.80", col = "blue", cex = 0.8)

plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = "Mean connectivity")

text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = 0.9,
     col = "red")

dev.off()

# selecciono beta optimo
soft_power <- sft$powerEstimate

if (is.na(soft_power)) {
  soft_power <- sft$fitIndices[which.max(sft$fitIndices[, 2]), 1]
  cat("[AVISO] R^2 no alcanzo 0.80, usando mejor disponible\n")
}

cat("[OK] Soft threshold (beta):", soft_power, "\n")
cat("  R^2:", round(sft$fitIndices[sft$fitIndices[,1] == soft_power, 2], 3), "\n")
cat("  Mean connectivity:", round(sft$fitIndices[sft$fitIndices[,1] == soft_power, 5], 1), "\n\n")

# --- construccion de red y deteccion de modulos ---
cat("[OK] Construyendo red de co-expresion...\n")

net <- blockwiseModules(
  datExpr_filtered,
  power = soft_power,
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "04_red_coexpresion/WGCNA/outputs/TOM",
  verbose = 3
)

cat("[OK] Red construida\n\n")

# --- visualizar modulos ---
moduleLabels <- net$colors
moduleColors <- labels2colors(moduleLabels)
n_modulos <- length(unique(moduleColors)) - 1  # sin gris

cat("[OK] Modulos identificados:", n_modulos, "\n")
cat("  Genes asignados:", sum(moduleColors != "grey"), "\n")
cat("  Genes no asignados (grey):", sum(moduleColors == "grey"), "\n\n")

module_table <- table(moduleColors)
print(module_table)
cat("\n")

# dendrograma
png("resultados_TFM_v2/fig_08_dendrogram_modules.png",
    width = 14, height = 8, units = "in", res = 300)

plotDendroAndColors(
  net$dendrograms[[1]],
  moduleColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Dendrograma de Genes y Asignacion de Modulos"
)

dev.off()

# --- calcular module eigengenes ---
MEs <- net$MEs
MEs_colors <- orderMEs(MEs)

cat("[OK] Eigengenes calculados para", ncol(MEs_colors), "modulos\n\n")

# --- identificar genes hub (preliminar, kME > 0.8) ---
# nota: este es un filtro preliminar para explorar, el filtro definitivo
# con kME > 0.5 y |GS_death| > 0.15 se aplica en el script 08
geneModuleMembership <- as.data.frame(cor(datExpr_filtered, MEs_colors, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(
  as.matrix(geneModuleMembership),
  nrow(datExpr_filtered)
))

hub_genes_list <- list()

for (module in unique(moduleColors)) {
  if (module == "grey") next

  inModule <- (moduleColors == module)
  modGenes <- colnames(datExpr_filtered)[inModule]

  if (length(modGenes) == 0) next

  module_column <- paste0("ME", module)

  if (!module_column %in% colnames(geneModuleMembership)) next

  kME_values <- geneModuleMembership[modGenes, module_column, drop = TRUE]

  if (length(kME_values) == 0) next

  isHub <- (abs(kME_values) > 0.8)

  if (sum(isHub, na.rm = TRUE) > 0) {
    hub_genes <- data.frame(
      Gene = modGenes[isHub],
      Module = module,
      kME = kME_values[isHub],
      stringsAsFactors = FALSE
    )
    hub_genes <- hub_genes[order(-abs(hub_genes$kME)), ]
    hub_genes_list[[module]] <- hub_genes
  }
}

if (length(hub_genes_list) > 0) {
  all_hub_genes <- do.call(rbind, hub_genes_list)
  rownames(all_hub_genes) <- NULL
} else {
  all_hub_genes <- data.frame()
}

cat("[OK] Genes HUB identificados:", nrow(all_hub_genes), "\n")
if (nrow(all_hub_genes) > 0) {
  print(table(all_hub_genes$Module))
  cat("\nTop 10 genes HUB:\n")
  print(head(all_hub_genes[order(-abs(all_hub_genes$kME)), ], 10))
}
cat("\n")

# --- heatmap de eigengenes ---
png("resultados_TFM_v2/fig_09_eigengene_heatmap.png",
    width = 10, height = 8, units = "in", res = 300)

plotEigengeneNetworks(
  MEs_colors,
  "Eigengene adjacency heatmap",
  marHeatmap = c(3,4,2,2),
  plotDendrograms = FALSE,
  xLabelsAngle = 90
)

dev.off()

# --- guardar resultados ---
module_assignment <- data.frame(
  Gene = colnames(datExpr_filtered),
  Module = moduleColors,
  stringsAsFactors = FALSE
)

write.table(
  module_assignment,
  "04_red_coexpresion/WGCNA/outputs/module_assignment.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  MEs_colors,
  "04_red_coexpresion/WGCNA/outputs/module_eigengenes.txt",
  sep = "\t", quote = FALSE, row.names = TRUE
)

write.table(
  geneModuleMembership,
  "04_red_coexpresion/WGCNA/outputs/gene_module_membership.txt",
  sep = "\t", quote = FALSE, row.names = TRUE
)

if (nrow(all_hub_genes) > 0) {
  write.table(
    all_hub_genes,
    "04_red_coexpresion/WGCNA/outputs/hub_genes.txt",
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

cat("[OK] Resultados guardados en 04_red_coexpresion/WGCNA/outputs/\n\n")

# --- resumen ---
cat("Resumen WGCNA:\n")
cat("  Muestras:", nrow(datExpr_filtered), "| Genes:", ncol(datExpr_filtered), "\n")
cat("  Beta:", soft_power, "| Modulos:", n_modulos, "| Hubs:", nrow(all_hub_genes), "\n")
