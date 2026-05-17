# 02_procesar_datos_GTEx.R
# Autor: Ines Torres
# Fecha: 06/02/2026
# Proceso los datos de GTEx descargados de UCSC Xena, extraigo solo pancreas normal

library(data.table)

# --- directorio de trabajo ---

setwd("~/TFM_PDAC")

# --- verifico que estan los archivos ---

archivos_necesarios <- c(
  "01_datos/raw/GTEx/gtex_gene_expected_count",
  "01_datos/raw/GTEx/GTEX_phenotype"
)

for (archivo in archivos_necesarios) {
  if (!file.exists(archivo)) {
    stop("Falta archivo: ", basename(archivo), ". Descargar de UCSC Xena primero.")
  }
}

# --- procesar phenotype ---

cat("Cargando phenotype...\n")

phenotype <- fread("01_datos/raw/GTEx/GTEX_phenotype",
                   header = TRUE,
                   sep = "\t")

cat("Phenotype:", nrow(phenotype), "muestras,", ncol(phenotype), "variables\n")

# busco columna de tejido
columnas_tejido <- grep("site|tissue|body_site",
                        colnames(phenotype),
                        ignore.case = TRUE,
                        value = TRUE)

col_tejido <- columnas_tejido[1]
cat("Columna de tejido:", col_tejido, "\n")

# top 10 tejidos
cat("\nTop 10 tejidos:\n")
tabla_tejidos <- sort(table(phenotype[[col_tejido]]), decreasing = TRUE)
print(head(tabla_tejidos, 10))

# filtro muestras de pancreas
idx_pancreas <- grep("pancreas",
                     phenotype[[col_tejido]],
                     ignore.case = TRUE)

cat("\nMuestras de pancreas:", length(idx_pancreas), "\n")

muestras_pancreas <- phenotype[idx_pancreas, ]

# IDs de pancreas
col_id <- colnames(phenotype)[1]
ids_pancreas <- muestras_pancreas[[col_id]]

# guardo IDs
write.table(
  ids_pancreas,
  file = "01_datos/metadata/GTEx_Pancreas_sample_IDs.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# guardo metadata de pancreas
write.table(
  muestras_pancreas,
  file = "01_datos/metadata/GTEx_Pancreas_metadata_complete.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Metadata guardada: 01_datos/metadata/GTEx_Pancreas_metadata_complete.txt\n")

# --- extraer matriz de expresion de pancreas ---

cat("\nCargando matriz de expresion completa (tarda 10-20 min)...\n")

expresion_gtex <- fread("01_datos/raw/GTEx/gtex_gene_expected_count",
                        header = TRUE,
                        sep = "\t",
                        showProgress = TRUE)

cat("Dimensiones totales:", nrow(expresion_gtex), "genes x", ncol(expresion_gtex) - 1, "muestras\n")

col_genes <- colnames(expresion_gtex)[1]

# --- filtrar solo muestras de pancreas ---

ids_encontrados <- ids_pancreas[ids_pancreas %in% colnames(expresion_gtex)]
cat("IDs de pancreas en la matriz:", length(ids_encontrados), "\n")

# puede haber diferencia entre IDs en phenotype y en la matriz por QC, es normal

if (length(ids_encontrados) > 100) {

  # filtro columnas: genes + pancreas
  matriz_pancreas <- expresion_gtex[, c(col_genes, ids_encontrados), with = FALSE]

  cat("Matriz filtrada:", nrow(matriz_pancreas), "genes x", ncol(matriz_pancreas) - 1, "muestras\n")

  # --- convertir a formato estandar y guardar ---

  matriz_df <- as.data.frame(matriz_pancreas)
  rownames(matriz_df) <- matriz_df[[col_genes]]
  matriz_df[[col_genes]] <- NULL

  cat("Guardando matriz...\n")

  write.table(
    matriz_df,
    file = "01_datos/processed/GTEx_Pancreas_raw_counts.txt",
    sep = "\t",
    quote = FALSE,
    row.names = TRUE,
    col.names = NA
  )

  saveRDS(matriz_df, "01_datos/processed/GTEx_Pancreas_matrix.rds")

  cat("Matriz guardada: 01_datos/processed/GTEx_Pancreas_raw_counts.txt\n")
  cat("Objeto RDS guardado: 01_datos/processed/GTEx_Pancreas_matrix.rds\n")

  # --- estadisticas ---

  cat("\nResumen matriz pancreas:\n")
  cat("- Genes:", nrow(matriz_df), "\n")
  cat("- Muestras:", ncol(matriz_df), "\n")
  cat("- Min:", min(matriz_df), "\n")
  cat("- Max:", max(matriz_df), "\n")
  cat("- Media:", round(mean(as.matrix(matriz_df)), 2), "\n")
  cat("- Mediana:", median(as.matrix(matriz_df)), "\n")

  cat("\nVista previa (5 genes x 3 muestras):\n")
  print(matriz_df[1:5, 1:3])

  # --- resumen final ---

  cat("\nArchivos generados:\n")
  cat("  01_datos/metadata/GTEx_Pancreas_sample_IDs.txt (", length(ids_encontrados), "IDs)\n")
  cat("  01_datos/metadata/GTEx_Pancreas_metadata_complete.txt\n")
  cat("  01_datos/processed/GTEx_Pancreas_raw_counts.txt (", nrow(matriz_df), "x", ncol(matriz_df), ")\n")
  cat("  01_datos/processed/GTEx_Pancreas_matrix.rds\n")

  cat("\nNota: diferencia entre phenotype (203) y matriz (", ncol(matriz_df), ") es por QC de GTEx.\n")

} else {
  cat("ERROR: No se encontraron suficientes muestras de pancreas\n")
  cat("Encontradas:", length(ids_encontrados), "(minimo esperado: 100)\n")
  stop("Error en el procesamiento")
}
