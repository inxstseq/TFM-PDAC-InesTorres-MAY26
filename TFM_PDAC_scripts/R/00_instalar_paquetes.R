# 00_instalar_paquetes.R
# Autor: Ines Torres
# Fecha: 06/02/2026
# Instalo todos los paquetes necesarios para el TFM (ejecutar solo una vez)

# --- instalar BiocManager ---

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# --- paquetes de Bioconductor ---

paquetes_bioc <- c(
  "TCGAbiolinks",
  "SummarizedExperiment",
  "DESeq2",
  "sva",
  "WGCNA",
  "biomaRt"
)

for (paq in paquetes_bioc) {
  if (!require(paq, character.only = TRUE, quietly = TRUE)) {
    BiocManager::install(paq, update = FALSE)
  }
}

# --- paquetes de CRAN ---

paquetes_cran <- c(
  "data.table",
  "dplyr",
  "ggplot2",
  "pheatmap",
  "RColorBrewer",
  "reshape2",
  "ggrepel",
  "survival",
  "survminer",
  "gridExtra"
)

for (paq in paquetes_cran) {
  if (!require(paq, character.only = TRUE, quietly = TRUE)) {
    install.packages(paq)
  }
}

# --- verificacion final ---

cat("\nVerificacion de paquetes:\n")

todos_paquetes <- c(paquetes_bioc, paquetes_cran)
instalados_ok <- TRUE

for (paq in todos_paquetes) {
  if (require(paq, character.only = TRUE, quietly = TRUE)) {
    cat("[OK]", paq, "\n")
  } else {
    cat("[ERROR]", paq, "\n")
    instalados_ok <- FALSE
  }
}

cat("\n")
if (instalados_ok) {
  cat("Todos los paquetes instalados correctamente.\n")
} else {
  cat("Algunos paquetes fallaron. Revisar mensajes de error.\n")
}
