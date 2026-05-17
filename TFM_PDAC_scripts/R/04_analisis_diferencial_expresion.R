# 04_analisis_diferencial_expresion.R
# Autor: Ines Torres
# Fecha: 06/02/2026
# Analisis diferencial de expresion con DESeq2 (tumor vs normal)

# --- librerias ---

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(ggrepel)

# --- configurar directorio ---

setwd("~/TFM_PDAC")

if (!dir.exists("03_expresion_diferencial/outputs")) {
  dir.create("03_expresion_diferencial/outputs", recursive = TRUE)
}

if (!dir.exists("resultados_TFM_v2")) {
  dir.create("resultados_TFM_v2", recursive = TRUE)
}

# --- cargar datos de fase 2 ---

dds <- readRDS("02_preprocesamiento/outputs_FINAL/DESeqDataSet.rds")

cat("[OK] Datos cargados -", nrow(dds), "genes,", ncol(dds), "muestras\n")
cat("Grupos:", paste(unique(dds$grupo), collapse = ", "), "\n")

# --- analisis diferencial ---

cat("Ejecutando DESeq2... (puede tardar 10-20 min)\n")
dds <- DESeq(dds)
cat("[OK] Analisis diferencial completado\n")

# --- extraer resultados ---

# tumor vs normal: log2FC positivo = mas alto en tumor
res <- results(dds,
               contrast = c("grupo", "TCGA_Tumor", "GTEx_Normal"),
               alpha = 0.05)

summary(res)

# --- filtrar DEGs ---

# criterios: padj < 0.05, |log2FC| > 1
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

# quito NAs
res_df <- res_df[!is.na(res_df$padj), ]

degs <- res_df[res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, ]

# clasifico UP y DOWN
degs_up <- degs[degs$log2FoldChange > 1, ]
degs_down <- degs[degs$log2FoldChange < -1, ]

# ordeno por significancia
degs_up <- degs_up[order(degs_up$padj), ]
degs_down <- degs_down[order(degs_down$padj), ]

cat("[OK] DEGs:", nrow(degs), "(UP:", nrow(degs_up), "/ DOWN:", nrow(degs_down), ")\n")

# --- verificar genes clave de PAAD ---

genes_clave <- c(
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

nombres_genes <- c("KRAS", "TP53", "CDKN2A", "SMAD4", "MUC1",
                   "MUC4", "MSLN", "S100A4", "CEACAM6")

verificacion <- data.frame(
  Gen = nombres_genes,
  ENSEMBL = genes_clave,
  log2FC = NA,
  padj = NA,
  Estado = NA,
  stringsAsFactors = FALSE
)

for (i in 1:length(genes_clave)) {
  gene_id <- genes_clave[i]

  if (gene_id %in% rownames(res_df)) {
    gene_data <- res_df[gene_id, ]
    verificacion$log2FC[i] <- round(gene_data$log2FoldChange, 2)
    verificacion$padj[i] <- format(gene_data$padj, scientific = TRUE, digits = 3)

    if (gene_data$padj < 0.05 & abs(gene_data$log2FoldChange) > 1) {
      if (gene_data$log2FoldChange > 1) {
        verificacion$Estado[i] <- "UP (tumor)"
      } else {
        verificacion$Estado[i] <- "DOWN (tumor)"
      }
    } else if (gene_data$padj < 0.05) {
      verificacion$Estado[i] <- "Sig. pero FC < 2x"
    } else {
      verificacion$Estado[i] <- "No significativo"
    }
  } else {
    verificacion$Estado[i] <- "No encontrado"
  }
}

print(verificacion)

# --- volcano plot ---

res_df$diffexpressed <- "NO"
res_df$diffexpressed[res_df$log2FoldChange > 1 & res_df$padj < 0.05] <- "UP"
res_df$diffexpressed[res_df$log2FoldChange < -1 & res_df$padj < 0.05] <- "DOWN"

res_df$log10padj <- -log10(res_df$padj)

# etiquetas para genes clave
res_df$gene_label <- ""
res_df$gene_label[res_df$gene_id %in% genes_clave] <- nombres_genes[match(
  res_df$gene_id[res_df$gene_id %in% genes_clave],
  genes_clave
)]

p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = log10padj, color = diffexpressed)) +
  geom_point(alpha = 0.4, size = 1.5) +
  scale_color_manual(
    values = c("UP" = "#E64B35", "DOWN" = "#4DBBD5", "NO" = "grey"),
    labels = c("DOWN" = paste0("DOWN (", nrow(degs_down), ")"),
               "NO" = "No significativo",
               "UP" = paste0("UP (", nrow(degs_up), ")"))
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_text_repel(
    aes(label = gene_label),
    size = 3,
    max.overlaps = 20,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50"
  ) +
  labs(
    title = "Volcano Plot: Expresion Diferencial en Cancer de Pancreas",
    subtitle = paste0("DEGs: ", nrow(degs), " (UP: ", nrow(degs_up), ", DOWN: ", nrow(degs_down), ")"),
    x = "log2 Fold Change (Tumor vs Normal)",
    y = "-log10(adjusted p-value)",
    color = "Estado"
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  xlim(c(-10, 10))

ggsave(
  "resultados_TFM_v2/fig_04_volcano_plot.png",
  plot = p_volcano,
  width = 12,
  height = 10,
  dpi = 300
)

cat("[OK] Volcano plot guardado\n")

# --- heatmap top 50 DEGs ---

# cargo matriz VST normalizada
matriz_vst <- read.table(
  "02_preprocesamiento/outputs_FINAL/matriz_vst_normalizada.txt",
  header = TRUE,
  row.names = 1,
  sep = "\t"
)

# top 25 UP + 25 DOWN
top_up <- head(degs_up$gene_id, 25)
top_down <- head(degs_down$gene_id, 25)
top_genes <- c(top_up, top_down)

heatmap_data <- matriz_vst[top_genes, ]

# anotaciones
metadata <- read.table(
  "02_preprocesamiento/outputs_FINAL/metadata_muestras.txt",
  header = TRUE,
  row.names = 1,
  sep = "\t"
)

annotation_col <- data.frame(
  Grupo = metadata$grupo,
  row.names = rownames(metadata)
)

annotation_colors <- list(
  Grupo = c("TCGA_Tumor" = "#E64B35", "GTEx_Normal" = "#4DBBD5")
)

annotation_row <- data.frame(
  Regulacion = c(rep("UP", 25), rep("DOWN", 25)),
  row.names = top_genes
)

annotation_colors$Regulacion <- c("UP" = "#E64B35", "DOWN" = "#4DBBD5")

png(
  "resultados_TFM_v2/fig_05_heatmap_top50_DEGs.png",
  width = 12,
  height = 14,
  units = "in",
  res = 300
)

pheatmap(
  heatmap_data,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  annotation_col = annotation_col,
  annotation_row = annotation_row,
  annotation_colors = annotation_colors,
  show_rownames = TRUE,
  show_colnames = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Top 50 Genes Diferencialmente Expresados\n(25 UP + 25 DOWN)",
  fontsize_row = 8,
  fontsize = 10
)

dev.off()

cat("[OK] Heatmap guardado\n")

# --- guardar resultados ---

write.table(
  res_df,
  file = "03_expresion_diferencial/outputs/resultados_completos_DESeq2.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  degs,
  file = "03_expresion_diferencial/outputs/DEGs_significativos.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  degs_up,
  file = "03_expresion_diferencial/outputs/DEGs_UP_regulated.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  degs_down,
  file = "03_expresion_diferencial/outputs/DEGs_DOWN_regulated.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  verificacion,
  file = "03_expresion_diferencial/outputs/verificacion_genes_clave.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("[OK] Todos los resultados guardados en 03_expresion_diferencial/\n")

# --- resumen final ---

cat("\n--- RESUMEN ---\n")
cat("Genes analizados:", nrow(dds), "\n")
cat("Genes validos:", nrow(res_df), "\n")
cat("DEGs totales:", nrow(degs), "\n")
cat("  UP:", nrow(degs_up), "/ DOWN:", nrow(degs_down), "\n")
cat("Figuras: volcano_plot.png, heatmap_top50_DEGs.png\n")
