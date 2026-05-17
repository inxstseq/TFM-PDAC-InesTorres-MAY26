#!/usr/bin/env python3
"""
11_integracion_DEG_hub_genes.py
Autor: Ines Torres
Fecha: 17/02/2026
Integro DEGs con hub genes de WGCNA para identificar genes clave candidatos.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Circle
import seaborn as sns
from math import sqrt, pi
import os
import warnings
warnings.filterwarnings('ignore')

# --- configuracion ---

BASE_DIR = os.path.expanduser("~/TFM_PDAC")

# umbrales
KME_THRESHOLD = 0.5
GS_THRESHOLD = 0.15
DEG_PADJ = 0.05
DEG_LOG2FC = 1.0

# carpetas de salida
INTEGRATION_DIR = os.path.join(BASE_DIR, "06_integracion")
VALIDATION_DIR = os.path.join(BASE_DIR, "07_validacion")
RESULTS_DIR = os.path.join(BASE_DIR, "resultados_TFM_v2")

for d in [os.path.join(INTEGRATION_DIR, "results"),
          os.path.join(INTEGRATION_DIR, "figures"),
          os.path.join(VALIDATION_DIR, "results"),
          RESULTS_DIR]:
    os.makedirs(d, exist_ok=True)

# --- cargo datos ---

# DEGs completos
degs_all = pd.read_csv(
    os.path.join(BASE_DIR, "03_expresion_diferencial/outputs/resultados_completos_DESeq2.txt"),
    sep="\t", index_col=0
)
print(f"[OK] DESeq2: {len(degs_all)} genes")

# DEGs significativos
degs_sig = degs_all[
    (degs_all['padj'] < DEG_PADJ) &
    (degs_all['log2FoldChange'].abs() > DEG_LOG2FC)
].copy()
degs_up = degs_sig[degs_sig['log2FoldChange'] > 0]
degs_down = degs_sig[degs_sig['log2FoldChange'] < 0]
print(f"[OK] DEGs significativos: {len(degs_sig)} (UP:{len(degs_up)}, DOWN:{len(degs_down)})")

# WGCNA: asignacion de modulos
module_assign = pd.read_csv(
    os.path.join(BASE_DIR, "04_red_coexpresion/WGCNA/outputs/module_assignment.txt"),
    sep="\t"
)
print(f"[OK] WGCNA: {len(module_assign)} genes, {module_assign['Module'].nunique()} modulos")

# WGCNA: module membership (kME)
kme_df = pd.read_csv(
    os.path.join(BASE_DIR, "04_red_coexpresion/WGCNA/outputs/gene_module_membership.txt"),
    sep="\t", index_col=0
)
print(f"[OK] kME: {len(kme_df)} genes, {kme_df.shape[1]} modulos")

# module eigengenes
me_df = pd.read_csv(
    os.path.join(BASE_DIR, "04_red_coexpresion/WGCNA/outputs/module_eigengenes.txt"),
    sep="\t", index_col=0
)
print(f"[OK] Eigengenes: {len(me_df)} muestras, {me_df.shape[1]} modulos")

# traits clinicos
traits = pd.read_csv(
    os.path.join(BASE_DIR, "01_datos/clinical/TCGA_PAAD_traits_numeric.txt"),
    sep="\t", index_col=0
)
traits.index = traits.index.str.replace("-", ".", regex=False)
print(f"[OK] Traits: {len(traits)} muestras, {list(traits.columns)}")

# tablas de conversion ENSEMBL -> Symbols
conv_up = pd.read_csv(
    os.path.join(BASE_DIR, "05_enriquecimiento_DEGs/conversion_table_UP.txt"), sep="\t"
)
conv_down = pd.read_csv(
    os.path.join(BASE_DIR, "05_enriquecimiento_DEGs/conversion_table_DOWN.txt"), sep="\t"
)
conv_all = pd.concat([conv_up, conv_down]).drop_duplicates(subset='ensembl_gene_id')
ensembl_to_symbol = dict(zip(conv_all['ensembl_gene_id'], conv_all['hgnc_symbol']))
print(f"[OK] Conversiones ENSEMBL->symbol: {len(ensembl_to_symbol)}")

# --- mapeo modulo numerico <-> color ---

module_sizes = module_assign[module_assign['Module'] != 'grey']['Module'].value_counts().sort_values(ascending=False)

color_order = ['grey'] + list(module_sizes.index)
number_to_color = {i: color for i, color in enumerate(color_order)}
color_to_number = {color: i for i, color in enumerate(color_order)}

# leo modulos significativos desde los resultados de script 08
# (en lugar de hardcodearlos, los detecto dinamicamente)
sig_modules_file = os.path.join(BASE_DIR, "04_red_coexpresion/WGCNA/traits/hub_genes_with_traits.txt")

if os.path.exists(sig_modules_file):
    hub_traits = pd.read_csv(sig_modules_file, sep="\t")
    # los modulos con hub genes son los relevantes
    sig_colors = hub_traits['Module'].unique().tolist()
    sig_colors = [c for c in sig_colors if c != 'grey']
    sig_modules_num = {}
    for color in sig_colors:
        num = color_to_number.get(color)
        if num is not None:
            sig_modules_num[f"ME{num}"] = color
    print(f"[OK] Modulos con hub genes (leidos de script 08): {sig_modules_num}")
else:
    # fallback: usar todos los modulos no-grey
    print("[WARN] No se encontro hub_genes_with_traits.txt, usando todos los modulos")
    sig_modules_num = {f"ME{num}": color for num, color in number_to_color.items() if color != 'grey'}

for num, color in sorted(number_to_color.items()):
    n_genes = len(module_assign[module_assign['Module'] == color])
    marker = " <- CON HUB GENES" if color in sig_modules_num.values() else ""
    print(f"  ME{num} = {color} ({n_genes} genes){marker}")

# --- calculo gene significance (GS) ---

# alineo muestras entre eigengenes y traits
common_samples = me_df.index.intersection(traits.index)
print(f"\n[OK] Muestras comunes (eigengenes + traits): {len(common_samples)}")

# cargo matriz de expresion VST (solo genes WGCNA)
wgcna_genes = set(module_assign['Gene'].values)

with open(os.path.join(BASE_DIR, "02_preprocesamiento/outputs_FINAL/matriz_vst_normalizada.txt")) as f:
    header = f.readline().strip().split('\t')

print(f"[OK] Total muestras en VST: {len(header)}")
print("  Filtrando genes WGCNA de la matriz...")

chunks = []
chunk_size = 5000
reader = pd.read_csv(
    os.path.join(BASE_DIR, "02_preprocesamiento/outputs_FINAL/matriz_vst_normalizada.txt"),
    sep="\t", index_col=0, chunksize=chunk_size
)
for chunk in reader:
    filtered = chunk[chunk.index.isin(wgcna_genes)]
    if len(filtered) > 0:
        chunks.append(filtered)

expr_matrix = pd.concat(chunks)
print(f"[OK] Genes cargados: {len(expr_matrix)} de {len(wgcna_genes)} esperados")

# transpongo: muestras en filas, genes en columnas
expr_T = expr_matrix.T
expr_T.index = expr_T.index.str.replace("-", ".", regex=False)

# alineo con traits
common_expr = expr_T.index.intersection(traits.index).intersection(common_samples)
print(f"[OK] Muestras para GS (expresion + traits): {len(common_expr)}")

expr_aligned = expr_T.loc[common_expr]
traits_aligned = traits.loc[common_expr]

# calculo correlacion gen-trait (death)
print("  Calculando correlacion gen-trait (death)...")
death_values = traits_aligned['death'].values

def _norm_cdf(x):
    """Aproximacion de CDF normal estandar."""
    return 0.5 * (1 + np.tanh(0.7978845608 * (x + 0.044715 * x**3)))

gs_death = {}
gs_pval_death = {}

for gene in expr_aligned.columns:
    gene_expr = expr_aligned[gene].values
    valid = ~np.isnan(gene_expr) & ~np.isnan(death_values)
    if valid.sum() < 10:
        continue

    x = gene_expr[valid]
    y = death_values[valid]

    x_mean = np.mean(x)
    y_mean = np.mean(y)
    x_diff = x - x_mean
    y_diff = y - y_mean

    numerator = np.sum(x_diff * y_diff)
    denominator = np.sqrt(np.sum(x_diff**2) * np.sum(y_diff**2))

    if denominator == 0:
        continue

    r = numerator / denominator
    n_valid = valid.sum()

    if abs(r) < 1:
        t_stat = r * sqrt((n_valid - 2) / (1 - r**2))
        p_val = 2 * (1 - _norm_cdf(abs(t_stat)))
    else:
        p_val = 0.0

    gs_death[gene] = r
    gs_pval_death[gene] = p_val

gs_series = pd.Series(gs_death)
print(f"[OK] GS calculado para {len(gs_series)} genes")
print(f"  GS medio: {gs_series.mean():.4f}, rango: [{gs_series.min():.4f}, {gs_series.max():.4f}]")
print(f"  Genes con |GS| > {GS_THRESHOLD}: {(gs_series.abs() > GS_THRESHOLD).sum()}")

# --- identifico hub genes ---
# criterios: kME > KME_THRESHOLD, |GS| > GS_THRESHOLD

hub_genes_all = []

for me_name, color in sig_modules_num.items():
    print(f"\n  Modulo {me_name} ({color}):")

    genes_in_module = module_assign[module_assign['Module'] == color]['Gene'].values
    genes_with_data = [g for g in genes_in_module if g in kme_df.index and g in gs_death]
    print(f"    Genes en modulo: {len(genes_in_module)}, con datos kME+GS: {len(genes_with_data)}")

    if len(genes_with_data) == 0:
        print(f"    [WARN] Sin genes con datos")
        continue

    kme_col = me_name
    if kme_col not in kme_df.columns:
        print(f"    [WARN] Columna {kme_col} no encontrada en kME")
        continue

    kme_vals = kme_df.loc[genes_with_data, kme_col]
    gs_vals = pd.Series({g: gs_death[g] for g in genes_with_data})

    is_hub = (kme_vals.abs() > KME_THRESHOLD) & (gs_vals.abs() > GS_THRESHOLD)
    hub_genes = [g for g in genes_with_data if is_hub.get(g, False)]

    print(f"    |kME|>{KME_THRESHOLD}: {(kme_vals.abs() > KME_THRESHOLD).sum()}")
    print(f"    |GS|>{GS_THRESHOLD}: {(gs_vals.abs() > GS_THRESHOLD).sum()}")
    print(f"    Hub genes (ambos criterios): {len(hub_genes)}")

    for gene in hub_genes:
        hub_genes_all.append({
            'Gene': gene,
            'Module_ME': me_name,
            'Module_color': color,
            'kME': kme_vals[gene],
            'GS_death': gs_vals[gene],
            'GS_pval': gs_pval_death.get(gene, np.nan)
        })

hub_df = pd.DataFrame(hub_genes_all)
if len(hub_df) > 0:
    hub_df = hub_df.sort_values('kME', ascending=False, key=abs)
    print(f"\n[OK] TOTAL HUB GENES: {len(hub_df)}")
else:
    print("\n[WARN] Sin hub genes con los umbrales actuales")

# --- interseccion DEGs y hub genes ---

if 'gene_id' in degs_sig.columns:
    deg_genes = set(degs_sig['gene_id'].dropna().values)
else:
    idx_str = degs_sig.index.astype(str)
    if idx_str.str.startswith('ENSG').any():
        deg_genes = set(idx_str)
    else:
        deg_genes = set(degs_sig.index)

hub_gene_set = set(hub_df['Gene'].values) if len(hub_df) > 0 else set()
key_genes = deg_genes & hub_gene_set

print(f"\n[OK] DEGs: {len(deg_genes)} | Hub genes: {len(hub_gene_set)} | Interseccion: {len(key_genes)}")

# si la interseccion esta vacia, intento con umbrales relajados
if len(key_genes) == 0 and len(hub_df) > 0:
    print("[WARN] Interseccion vacia. Probando umbral relajado...")
    degs_relaxed = degs_all[degs_all['padj'] < DEG_PADJ]
    if 'gene_id' in degs_relaxed.columns:
        deg_genes_relaxed = set(degs_relaxed['gene_id'].dropna().values)
    else:
        deg_genes_relaxed = set(degs_relaxed.index.astype(str))

    key_genes_relaxed = deg_genes_relaxed & hub_gene_set
    print(f"  DEGs relajados (solo padj<0.05): {len(deg_genes_relaxed)}")
    print(f"  Interseccion relajada: {len(key_genes_relaxed)}")

# --- tabla de genes clave ---

if len(key_genes) > 0:
    key_genes_data = []

    if 'gene_id' in degs_sig.columns:
        degs_sig_indexed = degs_sig.set_index('gene_id', drop=False)
    else:
        degs_sig_indexed = degs_sig.copy()

    for gene in key_genes:
        if gene in degs_sig_indexed.index:
            deg_info = degs_sig_indexed.loc[gene]
            if isinstance(deg_info, pd.DataFrame):
                deg_info = deg_info.iloc[0]
        else:
            continue

        hub_info = hub_df[hub_df['Gene'] == gene].iloc[0]
        symbol = ensembl_to_symbol.get(gene, "NA")

        key_genes_data.append({
            'ensembl_id': gene,
            'gene_symbol': symbol,
            'log2FoldChange': deg_info.get('log2FoldChange', np.nan),
            'padj': deg_info.get('padj', np.nan),
            'direction': 'UP' if deg_info.get('log2FoldChange', 0) > 0 else 'DOWN',
            'module_ME': hub_info['Module_ME'],
            'module_color': hub_info['Module_color'],
            'kME': hub_info['kME'],
            'GS_death': hub_info['GS_death'],
            'GS_pval': hub_info['GS_pval']
        })

    key_genes_df = pd.DataFrame(key_genes_data)
    key_genes_df = key_genes_df.sort_values('kME', ascending=False, key=abs)

    print(f"\n[OK] Genes clave identificados: {len(key_genes_df)}")
    print(key_genes_df[['gene_symbol', 'ensembl_id', 'log2FoldChange', 'padj',
                         'module_color', 'kME', 'GS_death', 'direction']].to_string(index=False))

    key_genes_df.to_csv(
        os.path.join(INTEGRATION_DIR, "results", "genes_clave_DEG_intersect_Hub.csv"),
        index=False
    )
    key_genes_df.to_csv(
        os.path.join(RESULTS_DIR, "tabla_genes_clave_DEG_intersect_Hub.csv"),
        index=False
    )
    print(f"[OK] Guardado: genes_clave_DEG_intersect_Hub.csv")

else:
    print("\n  Interseccion estricta vacia, guardo hub genes como candidatos...")
    key_genes_df = hub_df.copy()
    key_genes_df['gene_symbol'] = key_genes_df['Gene'].map(ensembl_to_symbol)

    # verifico cuales son DEGs (sin filtro de FC)
    degs_padj = degs_all[degs_all['padj'] < DEG_PADJ]
    if 'gene_id' in degs_padj.columns:
        deg_set_padj = set(degs_padj['gene_id'].dropna().values)
    else:
        deg_set_padj = set(degs_padj.index.astype(str))

    key_genes_df['is_DEG_padj'] = key_genes_df['Gene'].isin(deg_set_padj)
    key_genes_df['is_DEG_strict'] = key_genes_df['Gene'].isin(deg_genes)

    if 'gene_id' in degs_all.columns:
        log2fc_map = dict(zip(degs_all['gene_id'], degs_all['log2FoldChange']))
        padj_map = dict(zip(degs_all['gene_id'], degs_all['padj']))
    else:
        log2fc_map = dict(zip(degs_all.index.astype(str), degs_all['log2FoldChange']))
        padj_map = dict(zip(degs_all.index.astype(str), degs_all['padj']))

    key_genes_df['log2FoldChange'] = key_genes_df['Gene'].map(log2fc_map)
    key_genes_df['padj'] = key_genes_df['Gene'].map(padj_map)
    key_genes_df['direction'] = key_genes_df['log2FoldChange'].apply(
        lambda x: 'UP' if x > 0 else 'DOWN' if x < 0 else 'NS'
    )

    key_genes_df.to_csv(
        os.path.join(INTEGRATION_DIR, "results", "hub_genes_candidatos.csv"),
        index=False
    )
    key_genes_df.to_csv(
        os.path.join(RESULTS_DIR, "tabla_hub_genes_candidatos.csv"),
        index=False
    )
    print(f"[OK] Guardado: hub_genes_candidatos.csv")
    print(f"  Hub genes: {len(key_genes_df)}, DEG(padj<0.05): {key_genes_df['is_DEG_padj'].sum()}, DEG estricto: {key_genes_df['is_DEG_strict'].sum()}")

# --- validacion con genes conocidos ---

known_genes = {
    'KRAS': 'ENSG00000133703',
    'TP53': 'ENSG00000141510',
    'CDKN2A': 'ENSG00000147889',
    'SMAD4': 'ENSG00000141646',
    'MUC1': 'ENSG00000185499',
    'MUC4': 'ENSG00000145113',
    'S100A4': 'ENSG00000196154',
    'CEACAM6': 'ENSG00000086548',
    'MSLN': 'ENSG00000102854',
    'KRT19': 'ENSG00000171345',
    'VIM': 'ENSG00000026025',
    'S100P': 'ENSG00000163993',
    'CEACAM5': 'ENSG00000105388',
}

validation_results = []
for symbol, ensembl in known_genes.items():
    in_degs = ensembl in deg_genes
    in_hub = ensembl in hub_gene_set
    in_key = ensembl in key_genes
    in_wgcna = ensembl in set(module_assign['Gene'].values)

    module_color = "NA"
    if in_wgcna:
        mod = module_assign[module_assign['Gene'] == ensembl]['Module'].values
        module_color = mod[0] if len(mod) > 0 else "NA"

    kme_val = np.nan
    gs_val = np.nan
    log2fc_val = np.nan
    padj_val = np.nan

    if ensembl in kme_df.index:
        # busco kME en el modulo al que pertenece este gen
        if module_color in color_to_number:
            me_name = f"ME{color_to_number[module_color]}"
            if me_name in kme_df.columns:
                kme_val = kme_df.loc[ensembl, me_name]

    if ensembl in gs_death:
        gs_val = gs_death[ensembl]

    if 'gene_id' in degs_all.columns:
        deg_match = degs_all[degs_all['gene_id'] == ensembl]
        if len(deg_match) > 0:
            log2fc_val = deg_match.iloc[0]['log2FoldChange']
            padj_val = deg_match.iloc[0]['padj']
    elif ensembl in degs_all.index:
        log2fc_val = degs_all.loc[ensembl, 'log2FoldChange']
        padj_val = degs_all.loc[ensembl, 'padj']

    validation_results.append({
        'gene_symbol': symbol,
        'ensembl_id': ensembl,
        'log2FoldChange': log2fc_val,
        'padj': padj_val,
        'in_DEGs': in_degs,
        'in_WGCNA': in_wgcna,
        'module': module_color,
        'in_Hub': in_hub,
        'in_Key_Genes': in_key,
        'kME': kme_val,
        'GS_death': gs_val
    })

val_df = pd.DataFrame(validation_results)

print("\nValidacion con genes conocidos:")
print(val_df[['gene_symbol', 'log2FoldChange', 'padj', 'in_DEGs', 'in_WGCNA',
              'module', 'in_Hub', 'in_Key_Genes']].to_string(index=False))

validated = val_df[val_df['in_Key_Genes'] == True]
in_hub_not_deg = val_df[(val_df['in_Hub'] == True) & (val_df['in_Key_Genes'] == False)]

print(f"\n  Validados (en genes clave): {len(validated)}")
if len(validated) > 0:
    print(f"    {', '.join(validated['gene_symbol'].values)}")

print(f"  En hub pero no DEG estricto: {len(in_hub_not_deg)}")
if len(in_hub_not_deg) > 0:
    print(f"    {', '.join(in_hub_not_deg['gene_symbol'].values)}")

# nuevos candidatos (no en lista conocida)
known_ensembl = set(known_genes.values())
if len(key_genes_df) > 0:
    gene_col = 'Gene' if 'Gene' in key_genes_df.columns else 'ensembl_id'
    novel_candidates = key_genes_df[~key_genes_df[gene_col].isin(known_ensembl)]
    print(f"\n  Nuevos candidatos a biomarcadores: {len(novel_candidates)}")
    if len(novel_candidates) > 0 and 'gene_symbol' in novel_candidates.columns:
        print(f"    {', '.join(novel_candidates['gene_symbol'].dropna().values[:20])}")

# guardo validacion
val_df.to_csv(
    os.path.join(VALIDATION_DIR, "results", "validacion_genes_conocidos.csv"),
    index=False
)
val_df.to_csv(
    os.path.join(RESULTS_DIR, "tabla_validacion_genes_conocidos.csv"),
    index=False
)
print(f"[OK] Guardado: validacion_genes_conocidos.csv")

# --- guardo datos complementarios ---

# GS para todos los genes WGCNA
gs_df = pd.DataFrame({
    'Gene': list(gs_death.keys()),
    'GS_death': list(gs_death.values()),
    'GS_pval': [gs_pval_death.get(g, np.nan) for g in gs_death.keys()]
})
gs_df.to_csv(
    os.path.join(INTEGRATION_DIR, "results", "gene_significance_death.csv"),
    index=False
)

# hub genes con anotaciones
hub_df_annotated = hub_df.copy()
if len(hub_df_annotated) > 0:
    hub_df_annotated['gene_symbol'] = hub_df_annotated['Gene'].map(ensembl_to_symbol)
    hub_df_annotated.to_csv(
        os.path.join(INTEGRATION_DIR, "results", "hub_genes_complete.csv"),
        index=False
    )

print("[OK] Guardados: gene_significance_death.csv, hub_genes_complete.csv")

# --- resumen ---

print(f"\n# --- resumen ---")
print(f"  DESeq2: {len(degs_all)} genes | DEGs: {len(degs_sig)}")
print(f"  WGCNA: {len(module_assign)} genes | Modulos sig: {len(sig_modules_num)}")
print(f"  Hub genes: {len(hub_df)} | Genes clave (DEG+Hub): {len(key_genes)}")
print(f"  Validados: {len(validated)}")
if len(key_genes_df) > 0:
    gene_col = 'Gene' if 'Gene' in key_genes_df.columns else 'ensembl_id'
    novel = key_genes_df[~key_genes_df[gene_col].isin(known_ensembl)]
    print(f"  Nuevos candidatos: {len(novel)}")

print("\n[OK] Pipeline completado")
