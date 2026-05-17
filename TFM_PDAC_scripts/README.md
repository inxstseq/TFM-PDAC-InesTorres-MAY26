# TFM - Biomarcadores pronósticos en PDAC

Scripts del Trabajo Fin de Máster del Máster en Bioinformática (VIU, curso 2025-2026).

**Autora:** Inés Torres Marcos  
**Directora:** Dra. Dulcenombre de María del Saz Navarro

## Sobre el proyecto

Análisis bioinformático para identificar genes candidatos a biomarcadores pronósticos en adenocarcinoma ductal pancreático (PDAC). Se usan datos de RNA-seq de TCGA-PAAD (183 muestras tumorales) y GTEx (167 muestras de páncreas normal).

El análisis incluye expresión diferencial (DESeq2), red de coexpresión génica (WGCNA) y un modelo de supervivencia LASSO-Cox que selecciona una firma de 12 genes.

## Requisitos

R (>= 4.2) con los paquetes: TCGAbiolinks, SummarizedExperiment, DESeq2, sva, WGCNA, biomaRt, data.table, dplyr, ggplot2, pheatmap, RColorBrewer, reshape2, ggrepel, survival, survminer, gridExtra.

Python (>= 3.9) con: pandas, numpy, matplotlib, seaborn, gseapy, networkx.

El script `R/00_instalar_paquetes.R` instala todo lo necesario de R.

## Orden de ejecución

1. `R/00_instalar_paquetes.R` - instala las librerías (solo hace falta una vez)
2. `R/01_descargar_datos_TCGA.R` - descarga RNA-seq de TCGA-PAAD
3. `R/02_procesar_datos_GTEx.R` - procesa los datos de GTEx (UCSC Xena)
4. `R/03_preprocesamiento_integracion.R` - integra TCGA + GTEx, ComBat y VST
5. `R/04_analisis_diferencial_expresion.R` - expresión diferencial con DESeq2
6. `R/05_convertir_ensembl_simbolos.R` - convierte IDs ENSEMBL a símbolos génicos
7. `R/06_procesar_datos_clinicos.R` - datos clínicos de TCGA-PAAD
8. `R/07_wgcna_red_coexpresion.R` - red de coexpresión y módulos
9. `R/08_wgcna_correlacion_traits.R` - correlación módulo-trait y genes hub
10. `Python/11_integracion_DEG_hub_genes.py` - intersección DEGs + hub genes
11. `R/09_analisis_supervivencia.R` - Kaplan-Meier y Cox univariante
12. `Python/15_machine_learning.py` - regresión LASSO-Cox (firma pronóstica)
13. `Python/12_enriquecimiento_genes_clave.py` - enriquecimiento GO/KEGG
14. `Python/13_visualizaciones_finales.py` - figuras para el manuscrito

La numeración de los scripts no es consecutiva porque algunos se descartaron durante el desarrollo, pero el orden de arriba es el correcto.

## Antes de ejecutar

Hay que cambiar la ruta del directorio de trabajo en cada script. Por defecto usan `~/TFM_PDAC`:

```r
setwd("~/TFM_PDAC")  # en R
```
```python
BASE_DIR = os.path.expanduser("~/TFM_PDAC")  # en Python
```

Los datos de GTEx se descargan manualmente desde UCSC Xena antes de ejecutar el script 02. Los datos clínicos se obtienen de cBioPortal.

## Estructura de carpetas esperada

El directorio de trabajo (`TFM_PDAC/`) debe tener esta estructura para que los scripts funcionen:

```
TFM_PDAC/
├── 01_datos/
│   ├── raw/
│   ├── processed/
│   ├── metadata/
│   └── clinical/
├── 02_preprocesamiento/
├── 03_expresion_diferencial/
├── 04_red_coexpresion/
├── 05_enriquecimiento_DEGs/
├── 06_integracion/
└── resultados_TFM_v2/
```

Las carpetas se crean automáticamente al ejecutar los scripts si no existen.
