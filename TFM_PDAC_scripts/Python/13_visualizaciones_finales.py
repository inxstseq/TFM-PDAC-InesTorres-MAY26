#!/usr/bin/env python3
"""
13_visualizaciones_finales.py
Autora: Ines Torres
Fecha: 17/02/2026
Genera figuras de integracion DEG-WGCNA para la seccion de resultados.
- fig_11: diagrama de filtrado (funnel) del pipeline (Figura 11 del TFM)
- fig_12: scatter kME vs GS por modulo (Figura 12 del TFM)
- fig_13: heatmap top 50 genes clave (Figura 13 del TFM)
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import seaborn as sns
import os
import warnings
warnings.filterwarnings('ignore')

# --- configuracion ---
BASE_DIR = os.path.expanduser("~/TFM_PDAC")
INTEGRATION_DIR = os.path.join(BASE_DIR, "06_integracion")
RESULTS_DIR = os.path.join(BASE_DIR, "resultados_TFM_v2")
os.makedirs(RESULTS_DIR, exist_ok=True)

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 11,
    'axes.titlesize': 14,
    'axes.labelsize': 12,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.3,
})

print("[OK] Iniciando generacion de visualizaciones")

# --- cargar datos ---
key_genes = pd.read_csv(os.path.join(INTEGRATION_DIR, "results", "genes_clave_DEG_intersect_Hub.csv"))

# deteccion dinamica de la columna de modulo
module_col = 'module_color' if 'module_color' in key_genes.columns else 'Module_color'
if module_col not in key_genes.columns:
    module_col = 'Module'

# modulos presentes en los genes clave (sin grey)
active_modules = [m for m in key_genes[module_col].unique() if m != 'grey' and pd.notna(m)]

# leo asignacion de modulos para calcular ME numbers
module_assign = pd.read_csv(
    os.path.join(BASE_DIR, "04_red_coexpresion/WGCNA/outputs/module_assignment.txt"),
    sep="\t"
)
module_sizes = module_assign[module_assign['Module'] != 'grey']['Module'].value_counts()
color_to_num = {color: i+1 for i, color in enumerate(module_sizes.index)}
color_to_num['grey'] = 0

module_me_map = {}
for mod in active_modules:
    me_num = color_to_num.get(mod)
    if me_num is not None:
        module_me_map[mod] = f"ME{me_num}"

# correlaciones modulo-trait (muerte)
trait_corr_file = os.path.join(BASE_DIR, "04_red_coexpresion/WGCNA/traits/module_trait_correlations_full.txt")
trait_corr_map = {}
if os.path.exists(trait_corr_file):
    corr_df = pd.read_csv(trait_corr_file, sep="\t")
    for _, row in corr_df.iterrows():
        me = row.get('Module', '')
        if 'death' in corr_df.columns:
            death_val = str(row['death']).replace('*', '')
            try:
                trait_corr_map[me] = float(death_val)
            except:
                pass

# paleta de colores para los modulos
color_palette = {
    'turquoise': '#1B9E77', 'blue': '#4472C4', 'brown': '#D95F02',
    'yellow': '#E6AB02', 'green': '#66A61E', 'red': '#E7298A',
    'black': '#333333', 'pink': '#F781BF', 'magenta': '#984EA3',
    'purple': '#9370DB', 'greenyellow': '#A6D854', 'tan': '#D2B48C',
    'salmon': '#FA8072', 'cyan': '#00BFC4', 'midnightblue': '#191970',
    'lightcyan': '#B0E0E6', 'grey60': '#999999', 'lightgreen': '#90EE90'
}

# datos de GS y kME
gs_data = pd.read_csv(os.path.join(INTEGRATION_DIR, "results", "gene_significance_death.csv"))
kme_df = pd.read_csv(
    os.path.join(BASE_DIR, "04_red_coexpresion/WGCNA/outputs/gene_module_membership.txt"),
    sep="\t", index_col=0
)

# mapeo ensembl -> simbolo para anotar los graficos
conv_up = pd.read_csv(os.path.join(BASE_DIR, "05_enriquecimiento_DEGs/conversion_table_UP.txt"), sep="\t")
conv_down = pd.read_csv(os.path.join(BASE_DIR, "05_enriquecimiento_DEGs/conversion_table_DOWN.txt"), sep="\t")
conv_all = pd.concat([conv_up, conv_down]).drop_duplicates(subset='ensembl_gene_id')
ensembl_to_symbol = dict(zip(conv_all['ensembl_gene_id'], conv_all['hgnc_symbol']))

print(f"[OK] Modulos activos: {module_me_map}")
print(f"[OK] Genes clave: {len(key_genes)}")

# --- figura 1: diagrama de filtrado del pipeline (Figura 11 del TFM) ---
print("[1/3] Diagrama de filtrado (funnel)...")

# leo los DEGs para obtener los numeros reales del pipeline
degs_all = pd.read_csv(
    os.path.join(BASE_DIR, "03_expresion_diferencial/outputs/resultados_completos_DESeq2.txt"),
    sep="\t", index_col=0
)
n_input = len(degs_all)
n_degs = len(degs_all[(degs_all['padj'] < 0.05) & (degs_all['log2FoldChange'].abs() > 1)])
n_degs_strict = len(degs_all[(degs_all['padj'] < 0.01) & (degs_all['log2FoldChange'].abs() > 2)])

# genes en WGCNA, modulos significativos, hub genes
hub_genes = pd.read_csv(os.path.join(INTEGRATION_DIR, "results", "hub_genes_complete.csv"))
n_wgcna = len(module_assign)
n_sig_modules = len(key_genes)  # genes en modulos significativos antes de filtrar hub
n_hub = len(hub_genes)
n_key = len(key_genes)

fig, ax = plt.subplots(figsize=(14, 11))

# datos del pipeline (de arriba a abajo)
steps = [
    (f"{n_input:,} genes analizados", "Input DESeq2 (TCGA + GTEx)", "#B0C4DE", n_input),
    (f"{n_degs:,} DEGs", "padj < 0.05, |log2FC| > 1", "#7BA7CC", n_degs),
    (f"{n_degs_strict:,} DEGs estrictos", "padj < 0.01, |log2FC| > 2", "#4A90D9", n_degs_strict),
    (f"{n_wgcna:,} genes WGCNA", "Filtro MAD >= P75", "#3B7DD8", n_wgcna),
    ("854 genes en modulos\nsignificativos", "ME4 (yellow) + ME6 (red)\n|r| > 0.20, p < 0.05 (BH)", "#E8A838", 854),
    (f"{n_hub} hub genes", "kME > 0.5, |GS| > 0.15", "#E07B24", n_hub),
]

n_steps = len(steps)
max_width = 0.70
y_start = 0.93
y_end = 0.38
y_positions = np.linspace(y_start, y_end, n_steps)
bar_height = (y_start - y_end) / n_steps * 0.65

log_vals = [np.log10(s[3] + 1) for s in steps]
max_log = max(log_vals)
widths = [max_width * (0.20 + 0.80 * (lv / max_log)) for lv in log_vals]

x_center = 0.45

for i, (label, desc, color, count) in enumerate(steps):
    y = y_positions[i]
    w = widths[i]

    rect = FancyBboxPatch(
        (x_center - w/2, y - bar_height/2), w, bar_height,
        boxstyle="round,pad=0.008",
        facecolor=color, edgecolor='white', linewidth=2, alpha=0.92
    )
    ax.add_patch(rect)

    ax.text(x_center, y + 0.002, label,
            ha='center', va='center', fontsize=12, fontweight='bold',
            color='white' if i >= 2 else '#1a1a1a',
            transform=ax.transAxes, zorder=5)

    ax.text(x_center + w/2 + 0.02, y,
            desc, ha='left', va='center', fontsize=9,
            color='#444444', style='italic',
            transform=ax.transAxes, zorder=5)

    # trapecio conector entre barras
    if i < n_steps - 1:
        w_next = widths[i + 1]
        y_next = y_positions[i + 1]
        y_bottom = y - bar_height/2
        y_top_next = y_next + bar_height/2

        trap_x = [x_center - w/2, x_center + w/2,
                   x_center + w_next/2, x_center - w_next/2]
        trap_y = [y_bottom, y_bottom, y_top_next, y_top_next]

        trap = plt.Polygon(list(zip(trap_x, trap_y)),
                           facecolor='#E8E8E8', edgecolor='none',
                           alpha=0.3, transform=ax.transAxes, zorder=1)
        ax.add_patch(trap)

# --- dos ramas paralelas desde hub genes ---
hub_y = y_positions[-1]
branch_y = 0.16
branch_height = 0.08
branch_w = 0.30

# rama izquierda: KM + Cox
left_x = x_center - 0.22
rect_left = FancyBboxPatch(
    (left_x - branch_w/2, branch_y - branch_height/2), branch_w, branch_height,
    boxstyle="round,pad=0.008",
    facecolor='#D94F30', edgecolor='white', linewidth=2, alpha=0.92
)
ax.add_patch(rect_left)
ax.text(left_x, branch_y + 0.005, '6 genes pronosticos',
        ha='center', va='center', fontsize=11, fontweight='bold',
        color='white', transform=ax.transAxes, zorder=5)
ax.text(left_x, branch_y - branch_height/2 - 0.025,
        'Kaplan-Meier p < 0.05',
        ha='center', va='top', fontsize=8.5, color='#555555', style='italic',
        transform=ax.transAxes)
ax.text(left_x, branch_y - branch_height/2 - 0.045,
        'Cox univariante (individual)',
        ha='center', va='top', fontsize=8.5, color='#555555', style='italic',
        transform=ax.transAxes)

# rama derecha: LASSO-Cox
right_x = x_center + 0.22
rect_right = FancyBboxPatch(
    (right_x - branch_w/2, branch_y - branch_height/2), branch_w, branch_height,
    boxstyle="round,pad=0.008",
    facecolor='#C0392B', edgecolor='white', linewidth=2, alpha=0.92
)
ax.add_patch(rect_right)
ax.text(right_x, branch_y + 0.005, '12 genes firma LASSO-Cox',
        ha='center', va='center', fontsize=11, fontweight='bold',
        color='white', transform=ax.transAxes, zorder=5)
ax.text(right_x, branch_y - branch_height/2 - 0.025,
        'Regularizacion L1 (multivariante)',
        ha='center', va='top', fontsize=8.5, color='#555555', style='italic',
        transform=ax.transAxes)
ax.text(right_x, branch_y - branch_height/2 - 0.045,
        'C-index = 0.58, p = 0.030',
        ha='center', va='top', fontsize=8.5, color='#555555', style='italic',
        transform=ax.transAxes)

# flechas conectoras
ax.annotate('', xy=(left_x, branch_y + branch_height/2 + 0.01),
            xytext=(x_center - 0.05, hub_y - bar_height/2),
            xycoords='axes fraction', textcoords='axes fraction',
            arrowprops=dict(arrowstyle='->', color='#D94F30', lw=2.5,
                          connectionstyle='arc3,rad=0.1'))

ax.annotate('', xy=(right_x, branch_y + branch_height/2 + 0.01),
            xytext=(x_center + 0.05, hub_y - bar_height/2),
            xycoords='axes fraction', textcoords='axes fraction',
            arrowprops=dict(arrowstyle='->', color='#C0392B', lw=2.5,
                          connectionstyle='arc3,rad=-0.1'))

ax.text(x_center, 0.285, 'Validacion clinica',
        ha='center', va='center', fontsize=10, fontweight='bold',
        color='#666666', transform=ax.transAxes,
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#F0F0F0',
                  edgecolor='#CCCCCC', alpha=0.8))

# titulo
ax.text(0.45, 0.99, 'Pipeline de Seleccion de Genes Candidatos a Biomarcadores',
        ha='center', va='top', fontsize=15, fontweight='bold',
        transform=ax.transAxes)
ax.text(0.45, 0.965, 'Adenocarcinoma Ductal de Pancreas (PDAC)',
        ha='center', va='top', fontsize=11, color='#555555',
        transform=ax.transAxes)

# flecha lateral
ax.annotate('', xy=(0.04, 0.35), xytext=(0.04, y_start),
            xycoords='axes fraction', textcoords='axes fraction',
            arrowprops=dict(arrowstyle='->', color='#888888', lw=2))
ax.text(0.02, (y_start + 0.35)/2, 'Filtrado\nprogresivo',
        ha='center', va='center', fontsize=9, color='#888888',
        rotation=90, transform=ax.transAxes)

ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis('off')

plt.tight_layout()
plt.savefig(os.path.join(RESULTS_DIR, "fig_11_funnel_pipeline_genes.png"),
            facecolor='white')
plt.close()
print("[OK] fig_11_funnel_pipeline_genes.png")

# --- figura 2: scatter kME vs GS por modulo (Figura 12 del TFM) ---
print("[2/3] Scatter kME vs GS...")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for idx, (me_name, color, ax) in enumerate(zip(
    list(module_me_map.values()),
    list(module_me_map.keys()),
    axes)):

    genes_mod = module_assign[module_assign['Module'] == color]['Gene'].values
    genes_with_data = [g for g in genes_mod if g in kme_df.index and g in gs_data['Gene'].values]

    if len(genes_with_data) == 0:
        continue

    gs_lookup = dict(zip(gs_data['Gene'], gs_data['GS_death']))

    kme_vals = kme_df.loc[genes_with_data, me_name].values
    gs_vals = [gs_lookup.get(g, 0) for g in genes_with_data]

    key_gene_set = set(key_genes['ensembl_id'].values)
    is_key = [g in key_gene_set for g in genes_with_data]

    not_key_idx = [i for i, k in enumerate(is_key) if not k]
    key_idx = [i for i, k in enumerate(is_key) if k]

    ax.scatter([kme_vals[i] for i in not_key_idx],
               [gs_vals[i] for i in not_key_idx],
               c='lightgray', alpha=0.5, s=30, label='Otros genes del modulo')

    ax.scatter([kme_vals[i] for i in key_idx],
               [gs_vals[i] for i in key_idx],
               c=color if color != 'blue' else '#4472C4', alpha=0.8, s=60,
               edgecolors='black', linewidth=0.5,
               label=f'Genes clave (n={len(key_idx)})')

    # anotar top 8 genes por kME
    top_n = 8
    if len(key_idx) > 0:
        key_kme = [(kme_vals[i], gs_vals[i], genes_with_data[i]) for i in key_idx]
        key_kme.sort(key=lambda x: abs(x[0]), reverse=True)

        for kme, gs, gene in key_kme[:top_n]:
            symbol = ensembl_to_symbol.get(gene, '')
            if symbol and str(symbol) != 'nan':
                ax.annotate(symbol, (kme, gs), fontsize=7,
                           fontweight='bold', color='black',
                           xytext=(5, 5), textcoords='offset points',
                           bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))

    # umbrales de seleccion de hub genes
    ax.axhline(y=0.15, color='red', linestyle='--', alpha=0.5, linewidth=1)
    ax.axhline(y=-0.15, color='red', linestyle='--', alpha=0.5, linewidth=1)
    ax.axvline(x=0.5, color='red', linestyle='--', alpha=0.5, linewidth=1)

    mod_corr = f"r={trait_corr_map.get(me_name, 0):.3f}"
    ax.set_title(f'Modulo {me_name} ({color})\nCorrelacion con muerte: {mod_corr}',
                fontsize=12, fontweight='bold')
    ax.set_xlabel('Module Membership (kME)', fontsize=11)
    ax.set_ylabel('Gene Significance (GS) - Death', fontsize=11)
    ax.legend(fontsize=9, loc='best')
    ax.grid(True, alpha=0.3)

plt.suptitle('Scatter Plot: kME vs GS para Modulos Significativos',
             fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(RESULTS_DIR, "fig_12_scatter_kME_vs_GS.png"))
plt.close()
print("[OK] fig_12_scatter_kME_vs_GS.png")

# --- figura 2: heatmap de genes clave top 50 (Figura 13 del TFM) ---
print("[3/3] Heatmap genes clave...")

top50 = key_genes.head(50).copy()

heatmap_data = top50[['kME', 'GS_death', 'log2FoldChange']].copy()
heatmap_data.index = top50.apply(
    lambda row: row['gene_symbol'] if pd.notna(row['gene_symbol']) and row['gene_symbol'] != 'NA'
    else row['ensembl_id'][:15], axis=1
)
heatmap_data.columns = ['kME\n(Module Membership)', 'GS\n(Gene Significance)', 'log2FC\n(Fold Change)']

fig, ax = plt.subplots(figsize=(8, 16))

sns.heatmap(
    heatmap_data,
    cmap='RdBu_r',
    center=0,
    annot=True,
    fmt='.2f',
    linewidths=0.5,
    linecolor='white',
    cbar_kws={'label': 'Valor', 'shrink': 0.5},
    ax=ax,
    annot_kws={'fontsize': 8}
)

# barra lateral con color del modulo
colors_side = []
for _, row in top50.iterrows():
    colors_side.append(color_palette.get(row[module_col], '#999999'))

for i, c in enumerate(colors_side):
    ax.add_patch(plt.Rectangle((-0.15, i), 0.1, 1,
                                transform=ax.get_yaxis_transform(),
                                color=c, clip_on=False))

ax.set_title('Top 50 Genes Clave (DEG ∩ Hub)\nOrdenados por |kME|',
             fontsize=14, fontweight='bold', pad=20)
ax.set_ylabel('')

legend_patches = [mpatches.Patch(color=color_palette.get(c, '#999999'),
                                    label=f'{module_me_map.get(c, c)} ({c})')
                    for c in module_me_map.keys()]
ax.legend(handles=legend_patches, loc='upper right',
          bbox_to_anchor=(1.0, 1.08), fontsize=9)

plt.tight_layout()
plt.savefig(os.path.join(RESULTS_DIR, "fig_13_heatmap_genes_clave_top50.png"))
plt.close()
print("[OK] fig_13_heatmap_genes_clave_top50.png")

# --- resumen ---
print(f"\n[OK] 3 figuras generadas en {RESULTS_DIR}/")
print("[OK] Visualizaciones finalizadas")
