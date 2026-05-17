# 03_preprocesamiento_integracion.R
# Autor: Ines Torres
# Fecha: 06/02/2026
# Preproceso e integro las matrices TCGA + GTEx, corrijo batch effects

library(sva)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(reshape2)

# --- directorio de trabajo y carpetas ---

setwd("~/TFM_PDAC")

if (!dir.exists("02_preprocesamiento/outputs_FINAL")) {
  dir.create("02_preprocesamiento/outputs_FINAL", recursive = TRUE)
}

# carpeta unica para todos los resultados finales (figuras, tablas, datos)
if (!dir.exists("resultados_TFM_v2")) {
  dir.create("resultados_TFM_v2", recursive = TRUE)
}

# --- cargar datos ---

cat("Cargando datos...\n")

tcga <- read.table(
  "01_datos/processed/TCGA_PAAD_raw_counts.txt",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

cat("TCGA:", nrow(tcga), "genes x", ncol(tcga), "muestras\n")

gtex <- read.table(
  "01_datos/processed/GTEx_Pancreas_raw_counts.txt",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

cat("GTEx:", nrow(gtex), "genes x", ncol(gtex), "muestras\n")

# --- quitar versiones de ENSEMBL IDs y agregar duplicados ---
# sin esto KRAS se pierde porque TCGA tiene .13 y GTEx tiene .11

library(data.table)

# proceso TCGA
cat("\nEliminando versiones de ENSEMBL IDs...\n")

tcga$gene_id_base <- gsub("\\..*", "", rownames(tcga))

duplicados_tcga <- table(tcga$gene_id_base)
genes_duplicados_tcga <- names(duplicados_tcga[duplicados_tcga > 1])
cat("TCGA - duplicados por version:", length(genes_duplicados_tcga), "\n")

tcga_dt <- as.data.table(tcga, keep.rownames = TRUE)
setnames(tcga_dt, "rn", "gene_id_original")

cols_to_sum <- setdiff(colnames(tcga_dt), c("gene_id_original", "gene_id_base"))

tcga_agregado <- tcga_dt[, lapply(.SD, sum),
                         by = gene_id_base,
                         .SDcols = cols_to_sum]

tcga_final <- as.data.frame(tcga_agregado)
rownames(tcga_final) <- tcga_final$gene_id_base
tcga_final$gene_id_base <- NULL

cat("TCGA genes unicos:", nrow(tcga_final), "\n")

# proceso GTEx
gtex$gene_id_base <- gsub("\\..*", "", rownames(gtex))

duplicados_gtex <- table(gtex$gene_id_base)
genes_duplicados_gtex <- names(duplicados_gtex[duplicados_gtex > 1])
cat("GTEx - duplicados por version:", length(genes_duplicados_gtex), "\n")

gtex_dt <- as.data.table(gtex, keep.rownames = TRUE)
setnames(gtex_dt, "rn", "gene_id_original")

cols_to_sum_gtex <- setdiff(colnames(gtex_dt), c("gene_id_original", "gene_id_base"))

gtex_agregado <- gtex_dt[, lapply(.SD, sum),
                         by = gene_id_base,
                         .SDcols = cols_to_sum_gtex]

gtex_final <- as.data.frame(gtex_agregado)
rownames(gtex_final) <- gtex_final$gene_id_base
gtex_final$gene_id_base <- NULL

cat("GTEx genes unicos:", nrow(gtex_final), "\n")

# reemplazo variables
tcga <- tcga_final
gtex <- gtex_final

# --- genes comunes ---

genes_comunes <- intersect(rownames(tcga), rownames(gtex))
cat("\nGenes en comun:", length(genes_comunes), "\n")

# --- filtrado de genes ---

tcga_comun <- tcga[genes_comunes, ]
gtex_comun <- gtex[genes_comunes, ]

# combino antes de filtrar
matriz_combinada_sin_filtrar <- cbind(tcga_comun, gtex_comun)

cat("Matriz combinada sin filtrar:", nrow(matriz_combinada_sin_filtrar), "genes x",
    ncol(matriz_combinada_sin_filtrar), "muestras\n")

# filtro 1: row sum >= 10
filtro1 <- rowSums(matriz_combinada_sin_filtrar) >= 10

# filtro 2: expresado en al menos 10% de muestras
umbral_muestras <- 0.10 * ncol(matriz_combinada_sin_filtrar)
filtro2 <- rowSums(matriz_combinada_sin_filtrar > 0) >= umbral_muestras

filtro_final <- filtro1 & filtro2

cat("Genes tras filtrado:", sum(filtro_final), "(",
    round(sum(filtro_final)/length(filtro_final)*100, 1), "%)\n")

matriz_combinada <- matriz_combinada_sin_filtrar[filtro_final, ]

# --- verificar genes clave de cancer pancreatico ---

genes_clave <- c("KRAS", "TP53", "CDKN2A", "SMAD4", "MUC1", "MUC4", "MSLN", "S100A4", "CEACAM6")
ensembl_ids <- c(
  "ENSG00000133703",  # KRAS
  "ENSG00000141510",  # TP53
  "ENSG00000147889",  # CDKN2A
  "ENSG00000141646",  # SMAD4
  "ENSG00000185499",  # MUC1
  "ENSG00000145113",  # MUC4
  "ENSG00000102854",  # MSLN
  "ENSG00000196154",  # S100A4
  "ENSG00000086548"   # CEACAM6
)

nombres_genes <- rownames(matriz_combinada)

resultados_verificacion <- data.frame(
  Gen = genes_clave,
  ENSEMBL = ensembl_ids,
  Presente = "NO",
  stringsAsFactors = FALSE
)

for (i in 1:length(ensembl_ids)) {
  if (ensembl_ids[i] %in% nombres_genes) {
    resultados_verificacion$Presente[i] <- "SI"
  }
}

cat("\nVerificacion genes clave:\n")
print(resultados_verificacion)

genes_encontrados <- sum(resultados_verificacion$Presente == "SI")
cat("Encontrados:", genes_encontrados, "de", length(genes_clave), "\n")

# --- vectores de grupo y batch ---

grupo <- c(
  rep("TCGA_Tumor", ncol(tcga_comun)),
  rep("GTEx_Normal", ncol(gtex_comun))
)

batch <- c(
  rep(1, ncol(tcga_comun)),
  rep(2, ncol(gtex_comun))
)

# --- correccion batch con ComBat ---

cat("\nAplicando ComBat...\n")

matriz_log <- log2(matriz_combinada + 1)

# uso mod = NULL porque batch y condicion estan confundidos
# (TCGA = tumor = batch1, GTEx = normal = batch2)
# no se puede meter mod = ~grupo porque ComBat no distingue batch de biologia
# es lo que se hace siempre con TCGA+GTEx (Leek et al. 2012)

matriz_combat <- ComBat(
  dat = matriz_log,
  batch = batch,
  mod = NULL,
  par.prior = TRUE,
  prior.plots = FALSE
)

cat("ComBat completado.\n")

# --- normalizacion VST ---

cat("Aplicando VST...\n")

metadata <- data.frame(
  sample_id = colnames(matriz_combinada),
  grupo = grupo,
  batch = batch,
  row.names = colnames(matriz_combinada)
)

matriz_enteros <- round(matriz_combinada)

dds <- DESeqDataSetFromMatrix(
  countData = matriz_enteros,
  colData = metadata,
  design = ~ grupo
)

vst_data <- vst(dds, blind = FALSE)
matriz_vst <- assay(vst_data)

cat("VST completado.\n")

# --- guardar datos ---

cat("\nGuardando archivos...\n")

write.table(
  matriz_combinada,
  file = "02_preprocesamiento/outputs_FINAL/matriz_combinada_raw.txt",
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

write.table(
  matriz_combat,
  file = "02_preprocesamiento/outputs_FINAL/matriz_combat.txt",
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

write.table(
  matriz_vst,
  file = "02_preprocesamiento/outputs_FINAL/matriz_vst_normalizada.txt",
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

write.table(
  metadata,
  file = "02_preprocesamiento/outputs_FINAL/metadata_muestras.txt",
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

saveRDS(dds, "02_preprocesamiento/outputs_FINAL/DESeqDataSet.rds")

cat("Archivos guardados en 02_preprocesamiento/outputs_FINAL/\n")

# --- graficos ---

cat("Generando graficos...\n")

# boxplot
set.seed(123)
genes_muestra <- sample(1:nrow(matriz_combinada), min(1000, nrow(matriz_combinada)))

datos_raw <- log2(matriz_combinada[genes_muestra, ] + 1)
datos_raw_long <- reshape2::melt(as.matrix(datos_raw))
colnames(datos_raw_long) <- c("Gene", "Sample", "Expression")
datos_raw_long$Grupo <- ifelse(grepl("TCGA", datos_raw_long$Sample), "TCGA", "GTEx")
datos_raw_long$Tipo <- "Raw (log2)"

datos_vst_long <- reshape2::melt(as.matrix(matriz_vst[genes_muestra, ]))
colnames(datos_vst_long) <- c("Gene", "Sample", "Expression")
datos_vst_long$Grupo <- ifelse(grepl("TCGA", datos_vst_long$Sample), "TCGA", "GTEx")
datos_vst_long$Tipo <- "VST"

datos_plot <- rbind(datos_raw_long, datos_vst_long)

p1 <- ggplot(datos_plot, aes(x = Sample, y = Expression, fill = Grupo)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~ Tipo, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c("TCGA" = "#E64B35", "GTEx" = "#4DBBD5")) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top"
  ) +
  labs(
    title = "Distribucion de expresion genica",
    subtitle = "Antes y despues de normalizacion VST",
    x = "Muestras",
    y = "Expresion",
    fill = "Dataset"
  )

ggsave(
  "resultados_TFM_v2/fig_01_boxplot_normalizacion.png",
  plot = p1,
  width = 10,
  height = 8,
  dpi = 300
)

# PCA
pca_antes <- prcomp(t(matriz_log), scale. = FALSE)
var_antes <- round(100 * summary(pca_antes)$importance[2, 1:2], 1)

df_pca_antes <- data.frame(
  PC1 = pca_antes$x[, 1],
  PC2 = pca_antes$x[, 2],
  Grupo = grupo,
  Tipo = "Antes de ComBat"
)

pca_despues <- prcomp(t(matriz_combat), scale. = FALSE)
var_despues <- round(100 * summary(pca_despues)$importance[2, 1:2], 1)

df_pca_despues <- data.frame(
  PC1 = pca_despues$x[, 1],
  PC2 = pca_despues$x[, 2],
  Grupo = grupo,
  Tipo = "Despues de ComBat"
)

df_pca_total <- rbind(df_pca_antes, df_pca_despues)

p2 <- ggplot(df_pca_total, aes(x = PC1, y = PC2, color = Grupo)) +
  geom_point(size = 2, alpha = 0.7) +
  facet_wrap(~ Tipo, scales = "free") +
  scale_color_manual(values = c("TCGA_Tumor" = "#E64B35", "GTEx_Normal" = "#4DBBD5")) +
  theme_bw() +
  theme(legend.position = "top") +
  labs(
    title = "PCA: Efecto de correccion batch",
    x = paste0("PC1 (", var_antes[1], "% / ", var_despues[1], "%)"),
    y = paste0("PC2 (", var_antes[2], "% / ", var_despues[2], "%)"),
    color = "Grupo"
  )

ggsave(
  "resultados_TFM_v2/fig_02_PCA_batch_correction.png",
  plot = p2,
  width = 12,
  height = 6,
  dpi = 300
)

# heatmap de correlacion
set.seed(456)
genes_heatmap <- sample(1:nrow(matriz_vst), min(500, nrow(matriz_vst)))
cor_matrix <- cor(matriz_vst[genes_heatmap, ], method = "pearson")

annotation_col <- data.frame(
  Grupo = grupo,
  row.names = colnames(matriz_vst)
)

annotation_colors <- list(
  Grupo = c("TCGA_Tumor" = "#E64B35", "GTEx_Normal" = "#4DBBD5")
)

png(
  "resultados_TFM_v2/fig_03_heatmap_correlacion.png",
  width = 10,
  height = 10,
  units = "in",
  res = 300
)

pheatmap(
  cor_matrix,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  show_rownames = FALSE,
  show_colnames = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Correlacion entre muestras"
)

dev.off()

cat("Graficos guardados en resultados_TFM_v2/\n")

# --- resumen final ---

cat("\nResumen:\n")
cat("- Matriz final:", nrow(matriz_vst), "genes x", ncol(matriz_vst), "muestras\n")
cat("  Tumor:", sum(grupo == "TCGA_Tumor"), "| Normal:", sum(grupo == "GTEx_Normal"), "\n")
cat("- Genes clave encontrados:", genes_encontrados, "de", length(genes_clave), "\n")
