# 05_convertir_ensembl_simbolos.R
# Autor: Ines Torres
# Fecha: 07/02/2026
# Conversion de ENSEMBL IDs a simbolos genicos (top 500 UP/DOWN) con biomaRt

# --- librerias ---

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!require("biomaRt", quietly = TRUE)) {
  BiocManager::install("biomaRt")
}

library(biomaRt)

# --- configurar ---

setwd("~/TFM_PDAC")

if (!dir.exists("05_enriquecimiento_DEGs")) {
  dir.create("05_enriquecimiento_DEGs", recursive = TRUE)
}

# --- conectar a Ensembl ---

options(timeout = 300)

cat("Conectando a Ensembl...\n")
tryCatch({
  mart <- useMart("ensembl",
                  dataset = "hsapiens_gene_ensembl",
                  host = "useast.ensembl.org")
  cat("[OK] Conectado (mirror USA)\n")
}, error = function(e1) {
  cat("[ERROR] USA no disponible, probando Asia...\n")
  tryCatch({
    mart <<- useMart("ensembl",
                     dataset = "hsapiens_gene_ensembl",
                     host = "asia.ensembl.org")
    cat("[OK] Conectado (mirror Asia)\n")
  }, error = function(e2) {
    cat("[ERROR] Asia no disponible, probando principal...\n")
    mart <<- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    cat("[OK] Conectado (servidor principal)\n")
  })
})

# --- cargar DEGs ---

degs_up <- read.table(
  "03_expresion_diferencial/outputs/DEGs_UP_regulated.txt",
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

degs_down <- read.table(
  "03_expresion_diferencial/outputs/DEGs_DOWN_regulated.txt",
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

cat("[OK] Cargados:", nrow(degs_up), "UP,", nrow(degs_down), "DOWN\n")

# --- seleccionar top 500 ---

# ordeno por padj y cojo los 500 mas significativos de cada grupo
degs_up_sorted <- degs_up[order(degs_up$padj), ]
degs_down_sorted <- degs_down[order(degs_down$padj), ]

if (nrow(degs_up_sorted) >= 500) {
  degs_up_top500 <- degs_up_sorted[1:500, ]
} else {
  degs_up_top500 <- degs_up_sorted
  cat("[AVISO] Genes UP: solo hay", nrow(degs_up_sorted), "(menos de 500)\n")
}

if (nrow(degs_down_sorted) >= 500) {
  degs_down_top500 <- degs_down_sorted[1:500, ]
} else {
  degs_down_top500 <- degs_down_sorted
  cat("[AVISO] Genes DOWN: solo hay", nrow(degs_down_sorted), "(menos de 500)\n")
}

# --- funcion de conversion ---

convert_ensembl_to_symbol <- function(ensembl_ids, mart_obj, description) {
  ensembl_clean <- gsub("\\..*", "", ensembl_ids)

  cat("Consultando biomaRt para", description, "...\n")

  all_results <- getBM(
    attributes = c('ensembl_gene_id', 'hgnc_symbol', 'gene_biotype', 'description'),
    filters = 'ensembl_gene_id',
    values = ensembl_clean,
    mart = mart_obj
  )

  return(all_results)
}

# --- convertir UP ---

conversion_up <- convert_ensembl_to_symbol(degs_up_top500$gene_id, mart, "genes UP")

# filtro solo protein-coding con simbolo valido
conversion_up_clean <- conversion_up[
  conversion_up$gene_biotype == "protein_coding" &
    conversion_up$hgnc_symbol != "",
]

cat("[OK] UP: enviados", nrow(degs_up_top500), "-> protein-coding:", nrow(conversion_up_clean), "\n")

# --- convertir DOWN ---

conversion_down <- convert_ensembl_to_symbol(degs_down_top500$gene_id, mart, "genes DOWN")

conversion_down_clean <- conversion_down[
  conversion_down$gene_biotype == "protein_coding" &
    conversion_down$hgnc_symbol != "",
]

cat("[OK] DOWN: enviados", nrow(degs_down_top500), "-> protein-coding:", nrow(conversion_down_clean), "\n")

# --- guardar tablas de conversion ---
# estas tablas las usan los scripts de Python para mapear ENSEMBL -> simbolo

write.table(
  conversion_up_clean,
  file = "05_enriquecimiento_DEGs/conversion_table_UP.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  conversion_down_clean,
  file = "05_enriquecimiento_DEGs/conversion_table_DOWN.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("[OK] Tablas de conversion guardadas en 05_enriquecimiento_DEGs/\n")

# --- resumen ---

symbols_up <- unique(conversion_up_clean$hgnc_symbol)
symbols_down <- unique(conversion_down_clean$hgnc_symbol)

cat("\n--- RESUMEN ---\n")
cat("UP protein-coding:", length(symbols_up), "simbolos\n")
cat("DOWN protein-coding:", length(symbols_down), "simbolos\n")
cat("Tasa conversion UP:", round(length(symbols_up)/nrow(degs_up_top500)*100, 1), "%\n")
cat("Tasa conversion DOWN:", round(length(symbols_down)/nrow(degs_down_top500)*100, 1), "%\n")
