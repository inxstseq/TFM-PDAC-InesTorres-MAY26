# 01_descargar_datos_TCGA.R
# Autor: Ines Torres
# Fecha: 06/02/2026
# Descargo datos de expresion genica y clinicos de TCGA-PAAD

library(TCGAbiolinks)
library(SummarizedExperiment)

# --- directorio de trabajo ---

setwd("~/TFM_PDAC")

# --- buscar y descargar datos de TCGA-PAAD ---

proyecto <- "TCGA-PAAD"

# query para RNA-seq, conteos crudos de STAR
query_paad <- GDCquery(
  project = proyecto,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

num_muestras <- nrow(query_paad$results[[1]])
cat("Muestras encontradas:", num_muestras, "\n")

# --- descargar archivos ---

cat("Descargando archivos (tarda 15-30 min)...\n")

GDCdownload(
  query = query_paad,
  directory = "01_datos/raw/"
)

cat("Descarga completada.\n")

# --- preparar los datos ---

cat("Preparando datos...\n")

# convierte archivos individuales en un SummarizedExperiment
datos_paad <- GDCprepare(
  query = query_paad,
  save = TRUE,
  save.filename = "01_datos/processed/TCGA_PAAD_SummarizedExperiment.rda",
  directory = "01_datos/raw/"
)

# --- extraer la matriz de expresion ---

# conteos crudos: filas = genes (ENSEMBL), columnas = muestras
matriz_expresion <- assay(datos_paad, "unstranded")

cat("Matriz de expresion:", nrow(matriz_expresion), "genes x", ncol(matriz_expresion), "muestras\n")
print(matriz_expresion[1:5, 1:3])

# --- guardar la matriz ---

write.table(
  matriz_expresion,
  file = "01_datos/processed/TCGA_PAAD_raw_counts.txt",
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

cat("Matriz guardada: 01_datos/processed/TCGA_PAAD_raw_counts.txt\n")

# --- procesar la informacion clinica ---

info_clinica <- colData(datos_paad)

cat("Info clinica:", nrow(info_clinica), "muestras,", ncol(info_clinica), "variables\n")

# algunas columnas son listas y write.table no las acepta,
# las convierto a texto
info_clinica_limpia <- data.frame(row.names = rownames(info_clinica))

for (col_name in colnames(info_clinica)) {

  columna <- info_clinica[[col_name]]

  if (is.list(columna)) {
    info_clinica_limpia[[col_name]] <- sapply(columna, function(x) {
      if (is.null(x) || length(x) == 0) {
        return(NA)
      } else {
        return(paste(x, collapse = "; "))
      }
    })
  } else {
    info_clinica_limpia[[col_name]] <- columna
  }
}

# guardo info clinica completa
write.table(
  info_clinica_limpia,
  file = "01_datos/metadata/TCGA_PAAD_clinical_complete.txt",
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

cat("Info clinica guardada: 01_datos/metadata/TCGA_PAAD_clinical_complete.txt\n")

# --- resumen con variables importantes ---

columnas_clave <- c(
  "barcode",
  "patient",
  "sample_type",
  "gender",
  "vital_status",
  "age_at_diagnosis",
  "tumor_stage"
)

columnas_existentes <- columnas_clave[columnas_clave %in% colnames(info_clinica_limpia)]

if (length(columnas_existentes) > 0) {
  resumen_clinico <- info_clinica_limpia[, columnas_existentes, drop = FALSE]

  write.table(
    resumen_clinico,
    file = "01_datos/metadata/TCGA_PAAD_clinical_summary.txt",
    sep = "\t",
    quote = FALSE,
    row.names = TRUE,
    col.names = NA
  )

  cat("Resumen clinico guardado: 01_datos/metadata/TCGA_PAAD_clinical_summary.txt\n")
}

# --- estadisticas basicas ---

if ("sample_type" %in% colnames(info_clinica_limpia)) {
  cat("\nTipos de muestra:\n")
  print(table(info_clinica_limpia$sample_type))
}

cat("\nResumen de conteos:\n")
cat("- Min:", min(matriz_expresion), "\n")
cat("- Max:", max(matriz_expresion), "\n")
cat("- Media:", round(mean(matriz_expresion), 2), "\n")
cat("- Mediana:", median(matriz_expresion), "\n")

genes_cero <- sum(rowSums(matriz_expresion) == 0)
cat("- Genes sin expresion:", genes_cero, "(", round(genes_cero/nrow(matriz_expresion)*100, 2), "%)\n")

# --- resumen final ---

cat("\nArchivos generados:\n")
cat("  01_datos/raw/TCGA-PAAD/ (archivos originales GDC)\n")
cat("  01_datos/processed/TCGA_PAAD_SummarizedExperiment.rda\n")
cat("  01_datos/processed/TCGA_PAAD_raw_counts.txt (", nrow(matriz_expresion), "x", ncol(matriz_expresion), ")\n")
cat("  01_datos/metadata/TCGA_PAAD_clinical_complete.txt\n")
cat("  01_datos/metadata/TCGA_PAAD_clinical_summary.txt\n")
