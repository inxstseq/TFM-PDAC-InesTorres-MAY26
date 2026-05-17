# 06_procesar_datos_clinicos.R
# Autor: Ines Torres
# Fecha: 10/02/2026
# Extraigo y limpio datos clinicos TCGA-PAAD (supervivencia, estadio, etc.)

library(dplyr, warn.conflicts = FALSE)

setwd("~/TFM_PDAC")

if (!dir.exists("01_datos/clinical")) {
  dir.create("01_datos/clinical", recursive = TRUE)
}

# --- cargar archivo ---

clinical_raw <- read.table(
  "01_datos/clinical/paad_tcga_clinical_data.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

cat("[OK] Archivo cargado:", nrow(clinical_raw), "pacientes\n")

# --- procesar variables clinicas ---

clinical_clean <- data.frame(
  patient_id = clinical_raw$patientId,
  age_years = as.numeric(clinical_raw$AGE),
  sex_male = ifelse(clinical_raw$SEX == "Male", 1, 0),
  stage_numeric = case_when(
    grepl("^Stage I$|^Stage I[^IV]", clinical_raw$AJCC_PATHOLOGIC_TUMOR_STAGE, ignore.case = TRUE) ~ 1,
    grepl("Stage II", clinical_raw$AJCC_PATHOLOGIC_TUMOR_STAGE, ignore.case = TRUE) ~ 2,
    grepl("Stage III", clinical_raw$AJCC_PATHOLOGIC_TUMOR_STAGE, ignore.case = TRUE) ~ 3,
    grepl("Stage IV", clinical_raw$AJCC_PATHOLOGIC_TUMOR_STAGE, ignore.case = TRUE) ~ 4,
    TRUE ~ NA_real_
  ),
  survival_months = as.numeric(clinical_raw$OS_MONTHS),
  death_event = ifelse(grepl("1:DECEASED", clinical_raw$OS_STATUS), 1, 0),
  stringsAsFactors = FALSE
)

# --- matching con metadata ---

metadata <- read.table(
  "02_preprocesamiento/outputs_FINAL/metadata_muestras.txt",
  header = TRUE,
  row.names = 1,
  sep = "\t"
)

# extraigo patient_id de los nombres de muestra
sample_names <- rownames(metadata)
sample_names_hyphens <- gsub("\\.", "-", sample_names)

patient_ids <- sapply(strsplit(sample_names_hyphens, "-"), function(x) {
  if (length(x) >= 3) {
    paste(x[1:3], collapse = "-")
  } else {
    NA
  }
})

metadata_df <- data.frame(
  original_rowname = sample_names,
  patient_id = patient_ids,
  grupo = metadata$grupo,
  stringsAsFactors = FALSE
)

clinical_merged <- merge(
  metadata_df,
  clinical_clean,
  by = "patient_id",
  all.x = TRUE
)

cat("[OK] Merge:", nrow(clinical_merged), "filas\n")

# --- eliminar duplicados en rownames ---

dup_rownames <- duplicated(clinical_merged$original_rowname)

if (any(dup_rownames)) {
  n_dup <- sum(dup_rownames)
  cat("[AVISO]", n_dup, "duplicados encontrados, anado sufijos\n")

  dup_names <- unique(clinical_merged$original_rowname[dup_rownames])

  for (dup_name in dup_names) {
    idx <- which(clinical_merged$original_rowname == dup_name)
    for (i in seq_along(idx)) {
      if (i > 1) {
        clinical_merged$original_rowname[idx[i]] <- paste0(dup_name, ".", i-1)
      }
    }
  }
}

if (any(duplicated(clinical_merged$original_rowname))) {
  stop("[ERROR] Aun hay duplicados en rownames")
}

rownames(clinical_merged) <- clinical_merged$original_rowname

n_matched <- sum(!is.na(clinical_merged$age_years))
cat("[OK] Muestras:", nrow(clinical_merged), "- con datos clinicos:", n_matched, "\n")

# --- crear traits numericos ---

tumor_samples <- clinical_merged[clinical_merged$grupo == "TCGA_Tumor", ]

traits_numeric <- data.frame(
  row.names = rownames(tumor_samples),
  age = tumor_samples$age_years,
  sex_male = tumor_samples$sex_male,
  stage = tumor_samples$stage_numeric,
  survival = tumor_samples$survival_months,
  death = tumor_samples$death_event
)

if (any(duplicated(rownames(traits_numeric)))) {
  stop("[ERROR] Duplicados en traits")
}

cat("[OK] Traits creados -", nrow(tumor_samples), "muestras tumorales\n")

valid_counts <- colSums(!is.na(traits_numeric))
print(valid_counts)
print(summary(traits_numeric))

# --- guardar ---

write.table(
  clinical_merged,
  "01_datos/clinical/TCGA_PAAD_clinical_complete.txt",
  sep = "\t", quote = FALSE, row.names = TRUE
)

write.table(
  traits_numeric,
  "01_datos/clinical/TCGA_PAAD_traits_numeric.txt",
  sep = "\t", quote = FALSE, row.names = TRUE
)

cat("[OK] Archivos guardados en 01_datos/clinical/\n")

# --- resumen ---

cat("\n--- RESUMEN ---\n")
cat("Muestras totales:", nrow(clinical_merged), "\n")
cat("Muestras tumorales:", nrow(tumor_samples), "\n")
cat("Valores validos por trait:\n")
for (trait in names(valid_counts)) {
  cat("  ", trait, ":", valid_counts[trait], "\n")
}
