#!/usr/bin/env python3
"""
13_visualizaciones_finales.py
Autora: Ines Torres
Fecha: 17/02/2026
Genera figuras finales para la tesis (seccion 3.9).
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Circle, FancyArrowPatch
import seaborn as sns
from math import sqrt, pi, cos, sin
import os
import warnings
warnings.filterwarnings('ignore')

# --- configuracion ---
BASE_DIR = os.path.expanduser("~/TFM_PDAC")
INTEGRATION_DIR = os.path.join(BASE_DIR, "06_integracion")
RESULTS_DIR = os.path.join(BASE_DIR, "resultados_TFM_v2")
os.makedirs(RESULTS_DIR, exist_ok=True)

# --- estilo de las figuras ---
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
hub_genes = pd.read_csv(os.path.join(INTEGRATION_DIR, "results", "hub_genes_complete.csv"))

# --- deteccion dinamica de modulos ---
module_col = 'module_color' if 'module_color' in key_genes.columns else 'Module_color'
if module_col not in key_genes.columns:
    module_col = 'Module'

# modulos presentes en los genes clave (excluyendo grey)
active_modules = [m for m in key_genes[module_col].unique() if m != 'grey' and pd.notna(m)]

# leo correlaciones modulo-trait para obtener ME numbers y correlaciones
trait_corr_file = os.path.join(BASE_DIR, "04_red_coexpresion/WGCNA/traits/module_trait_correlations_full.txt")
module_assign = pd.read_csv(os.path.join(BASE_DIR, "04_red_coexpresion/WGCNA/outputs/module_assignment.txt"), sep="\t")
module_sizes = module_assign[module_assign['Module'] != 'grey']['Module'].value_counts()
color_to_num = {color: i+1 for i, color in enumerate(module_sizes.index)}
color_to_num['grey'] = 0

module_me_map = {}  # color -> ME name
trait_corr_map = {}  # ME name -> correlation with death
for mod in active_modules:
    me_num = color_to_num.get(mod)
    if me_num is not None:
        me_name = f"ME{me_num}"
        module_me_map[mod] = me_name

# leer correlaciones si existe el archivo
if os.path.exists(trait_corr_file):
    corr_df = pd.read_csv(trait_corr_file, sep="\t")
    for _, row in corr_df.iterrows():
        me = row.get('Module', row.name if hasattr(row, 'name') else '')
        if 'death' in corr_df.columns:
            # extraer numero de la correlacion (puede tener asteriscos)
            death_val = str(row['death']).replace('*', '')
            try:
                trait_corr_map[me] = float(death_val)
            except:
                pass

# paleta de colores
color_palette = {
    'turquoise': '#1B9E77', 'blue': '#4472C4', 'brown': '#D95F02',
    'yellow': '#E6AB02', 'green': '#66A61E', 'red': '#E7298A',
    'black': '#333333', 'pink': '#F781BF', 'magenta': '#984EA3',
    'purple': '#9370DB', 'greenyellow': '#A6D854', 'tan': '#D2B48C',
    'salmon': '#FA8072', 'cyan': '#00BFC4', 'midnightblue': '#191970',
    'lightcyan': '#B0E0E6', 'grey60': '#999999', 'lightgreen': '#90EE90'
}

print(f"[OK] Modulos activos: {module_me_map}")
gs_data = pd.read_csv(os.path.join(INTEGRATION_DIR, "results", "gene_significance_death.csv"))
validation = pd.read_csv(os.path.join(BASE_DIR, "07_validacion", "results", "validacion_genes_conocidos.csv"))

degs_all = pd.read_csv(
    os.path.join(BASE_DIR, "03_expresion_diferencial/outputs/resultados_completos_DESeq2.txt"),
    sep="\t", index_col=0
)

module_assign = pd.read_csv(
    os.path.join(BASE_DIR, "04_red_coexpresion/WGCNA/outputs/module_assignment.txt"),
    sep="\t"
)

kme_df = pd.read_csv(
    os.path.join(BASE_DIR, "04_red_coexpresion/WGCNA/outputs/gene_module_membership.txt"),
    sep="\t", index_col=0
)

conv_up = pd.read_csv(os.path.join(BASE_DIR, "05_enriquecimiento_DEGs/conversion_table_UP.txt"), sep="\t")
conv_down = pd.read_csv(os.path.join(BASE_DIR, "05_enriquecimiento_DEGs/conversion_table_DOWN.txt"), sep="\t")
conv_all = pd.concat([conv_up, conv_down]).drop_duplicates(subset='ensembl_gene_id')
ensembl_to_symbol = dict(zip(conv_all['ensembl_gene_id'], conv_all['hgnc_symbol']))

print(f"[OK] Genes clave: {len(key_genes)}, Hub genes: {len(hub_genes)}")

# --- figura 1: diagrama de Venn (DEGs inter Hub Genes) ---
print("[1/6] Diagrama de Venn...")

fig, ax = plt.subplots(1, 1, figsize=(8, 6))

if 'gene_id' in degs_all.columns:
    degs_sig = degs_all[(degs_all['padj'] < 0.05) & (degs_all['log2FoldChange'].abs() > 1)]
    n_degs = len(degs_sig)
    deg_genes_set = set(degs_sig['gene_id'].dropna().values)
else:
    degs_sig = degs_all[(degs_all['padj'] < 0.05) & (degs_all['log2FoldChange'].abs() > 1)]
    n_degs = len(degs_sig)
    deg_genes_set = set(degs_sig.index)

hub_genes_set = set(hub_genes['Gene'].values)
intersection = deg_genes_set & hub_genes_set

n_only_degs = len(deg_genes_set - hub_genes_set)
n_only_hub = len(hub_genes_set - deg_genes_set)
n_intersection = len(intersection)

circle1 = Circle((-0.3, 0), 1.0, alpha=0.35, facecolor='#4472C4', edgecolor='#2F5496', linewidth=2)
circle2 = Circle((0.3, 0), 0.6, alpha=0.35, facecolor='#ED7D31', edgecolor='#C55A11', linewidth=2)
ax.add_patch(circle1)
ax.add_patch(circle2)

ax.text(-0.65, 0, f'{n_only_degs:,}', fontsize=18, fontweight='bold', ha='center', va='center', color='#2F5496')
ax.text(0.05, 0, f'{n_intersection}', fontsize=22, fontweight='bold', ha='center', va='center', color='#7B2D0D')
ax.text(0.6, 0, f'{n_only_hub}', fontsize=18, fontweight='bold', ha='center', va='center', color='#C55A11')

ax.text(-0.8, 1.15, f'DEGs\n(padj<0.05, |log2FC|>1)\nn={n_degs:,}',
        fontsize=11, ha='center', va='center', fontweight='bold', color='#2F5496')
ax.text(0.7, 0.8, f'Hub Genes\n(kME>0.5, |GS|>0.15)\nn={len(hub_genes_set)}',
        fontsize=11, ha='center', va='center', fontweight='bold', color='#C55A11')

ax.set_xlim(-1.6, 1.3)
ax.set_ylim(-1.3, 1.4)
ax.set_aspect('equal')
ax.axis('off')
ax.set_title('Integracion de Resultados: DEGs ∩ Hub Genes WGCNA\nGenes Clave Candidatos a Biomarcadores',
             fontsize=14, fontweight='bold', pad=15)

textstr = f'Genes clave identificados: {n_intersection}\n(presentes en ambos conjuntos)'
props = dict(boxstyle='round,pad=0.5', facecolor='#FFF2CC', edgecolor='#7B2D0D', alpha=0.8)
ax.text(0.05, -1.1, textstr, fontsize=11, ha='center', va='center', bbox=props, fontweight='bold')

plt.savefig(os.path.join(RESULTS_DIR, "fig_11_venn_DEGs_HubGenes.png"))
plt.savefig(os.path.join(INTEGRATION_DIR, "figures", "01_venn_DEGs_HubGenes.png"))
plt.close()
print("[OK] 01_venn_DEGs_HubGenes.png")

# --- figura 2: scatter kME vs GS (por modulo) ---
print("[2/6] Scatter kME vs GS...")

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
plt.savefig(os.path.join(INTEGRATION_DIR, "figures", "02_scatter_kME_vs_GS.png"))
plt.close()
print("[OK] 02_scatter_kME_vs_GS.png")

# --- figura 3: heatmap de genes clave (top 50) ---
print("[3/6] Heatmap genes clave...")

top50 = key_genes.head(50).copy()

heatmap_data = top50[['kME', 'GS_death', 'log2FoldChange']].copy()
heatmap_data.index = top50.apply(
    lambda row: row['gene_symbol'] if pd.notna(row['gene_symbol']) and row['gene_symbol'] != 'NA'
    else row['ensembl_id'][:15], axis=1
)
heatmap_data.columns = ['kME\n(Module Membership)', 'GS\n(Gene Significance)', 'log2FC\n(Fold Change)']

fig, ax = plt.subplots(figsize=(8, 16))

from matplotlib.colors import TwoSlopeNorm

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
plt.savefig(os.path.join(INTEGRATION_DIR, "figures", "03_heatmap_genes_clave_top50.png"))
plt.close()
print("[OK] 03_heatmap_genes_clave_top50.png")

# --- figura 4: red de coexpresion (genes clave) ---
print("[4/6] Red de coexpresion...")

fig, ax = plt.subplots(figsize=(14, 12))

# selecciono top genes por modulo dinamicamente
network_parts = []
for mod_color in module_me_map.keys():
    mod_genes = key_genes[key_genes[module_col] == mod_color].head(20)
    network_parts.append(mod_genes)
network_genes = pd.concat(network_parts) if network_parts else key_genes.head(40)

# posiciones dinamicas por modulo (circular)
positions = {}
labels = {}
mod_list = list(module_me_map.keys())
n_mods = len(mod_list)
spacing = 7.0 / max(n_mods, 1)

for mod_idx, mod_color in enumerate(mod_list):
    mod_genes = network_genes[network_genes[module_col] == mod_color]
    n_genes = len(mod_genes)
    if n_genes == 0:
        continue
    # distribuir clusters horizontalmente
    center_x = -3.5 + mod_idx * spacing
    center_y = 0
    for i, (_, row) in enumerate(mod_genes.iterrows()):
        angle = 2 * pi * i / n_genes
        r = 1.5 + 0.3 * abs(row['kME'])
        x = center_x + r * cos(angle)
        y = center_y + r * sin(angle)
        gene_id = row.get('ensembl_id', row.get('Gene', ''))
        positions[gene_id] = (x, y)
        symbol = row.get('gene_symbol', '')
        if pd.isna(symbol) or symbol == 'NA':
            symbol = ''
        labels[gene_id] = symbol

# Aristas (conexiones intra-modulo por similitud de kME)
gene_id_col = 'ensembl_id' if 'ensembl_id' in network_genes.columns else 'Gene'
for idx1, (_, row1) in enumerate(network_genes.iterrows()):
    for idx2, (_, row2) in enumerate(network_genes.iterrows()):
        if idx2 <= idx1:
            continue
        if row1[module_col] != row2[module_col]:
            continue
        if abs(row1['kME']) > 0.7 and abs(row2['kME']) > 0.7:
            g1 = row1[gene_id_col]
            g2 = row2[gene_id_col]
            if g1 in positions and g2 in positions:
                pos1 = positions[g1]
                pos2 = positions[g2]
                alpha = min(abs(row1['kME'] * row2['kME']), 1.0) * 0.3
                c = color_palette.get(row1[module_col], '#999999')
                ax.plot([pos1[0], pos2[0]], [pos1[1], pos2[1]],
                       color=c, alpha=alpha, linewidth=0.5)

# Nodos
for _, row in network_genes.iterrows():
    gene_id = row[gene_id_col]
    if gene_id not in positions:
        continue
    pos = positions[gene_id]
    c = color_palette.get(row[module_col], '#999999')
    size = 100 + 300 * abs(row['kME'])
    ax.scatter(pos[0], pos[1], s=size, c=c, alpha=0.8,
              edgecolors='black', linewidth=0.8, zorder=5)
    label = labels.get(gene_id, '')
    if label:
        ax.annotate(label, pos, fontsize=7, fontweight='bold',
                   ha='center', va='bottom',
                   xytext=(0, 8), textcoords='offset points',
                   bbox=dict(boxstyle='round,pad=0.15', facecolor='white', alpha=0.8),
                   zorder=6)

# Etiquetas de modulo
for mod_idx, mod_color in enumerate(mod_list):
    center_x = -3.5 + mod_idx * spacing
    me_name = module_me_map.get(mod_color, mod_color)
    corr = trait_corr_map.get(me_name, 0)
    ax.text(center_x, -3.3, f'{me_name} ({mod_color})\nr={corr:.3f} con muerte',
           fontsize=10, fontweight='bold', ha='center',
           bbox=dict(facecolor=color_palette.get(mod_color, '#DDD'),
                    edgecolor='black', boxstyle='round,pad=0.5', alpha=0.3))

ax.set_xlim(-7, 7)
ax.set_ylim(-4.5, 4)
ax.axis('off')
ax.set_title('Red de Coexpresion: Top 20 Genes Clave por Modulo\n(Tamano del nodo proporcional a kME)',
             fontsize=14, fontweight='bold')

for kme_val, label in [(0.9, 'kME=0.9'), (0.7, 'kME=0.7'), (0.5, 'kME=0.5')]:
    s = 100 + 300 * kme_val
    ax.scatter([], [], s=s, c='gray', alpha=0.6, edgecolors='black', linewidth=0.5, label=label)
ax.legend(title='Module Membership', loc='upper left', fontsize=9, title_fontsize=10)

plt.savefig(os.path.join(RESULTS_DIR, "fig_14_network_coexpression.png"))
plt.savefig(os.path.join(INTEGRATION_DIR, "figures", "04_network_coexpression.png"))
plt.close()
print("[OK] 04_network_coexpression.png")

# --- figura 5: validacion - genes conocidos ---
print("[5/6] Tabla de validacion...")

fig, ax = plt.subplots(figsize=(14, 6))
ax.axis('off')

val_display = validation[['gene_symbol', 'log2FoldChange', 'padj', 'in_DEGs',
                           'in_WGCNA', 'module', 'in_Hub', 'in_Key_Genes']].copy()

val_display['log2FoldChange'] = val_display['log2FoldChange'].round(2)
val_display['padj'] = val_display['padj'].apply(lambda x: f'{x:.2e}' if pd.notna(x) else 'NA')
val_display['in_DEGs'] = val_display['in_DEGs'].map({True: '✓', False: '✗'})
val_display['in_WGCNA'] = val_display['in_WGCNA'].map({True: '✓', False: '✗'})
val_display['in_Hub'] = val_display['in_Hub'].map({True: '✓', False: '✗'})
val_display['in_Key_Genes'] = val_display['in_Key_Genes'].map({True: '✓', False: '✗'})

col_labels = ['Gen', 'log2FC', 'padj', 'DEG', 'WGCNA', 'Modulo', 'Hub', 'Gen Clave']

table = ax.table(
    cellText=val_display.values,
    colLabels=col_labels,
    loc='center',
    cellLoc='center'
)
table.auto_set_font_size(False)
table.set_fontsize(10)
table.auto_set_column_width(col=list(range(len(col_labels))))

for j in range(len(col_labels)):
    table[0, j].set_facecolor('#4472C4')
    table[0, j].set_text_props(color='white', fontweight='bold')

for i in range(len(val_display)):
    gene = validation.iloc[i]
    for j in range(len(col_labels)):
        if gene['in_Key_Genes']:
            table[i+1, j].set_facecolor('#C6EFCE')
        elif gene['in_Hub']:
            table[i+1, j].set_facecolor('#FFEB9C')
        elif gene['in_DEGs']:
            table[i+1, j].set_facecolor('#FCE4D6')

table.scale(1, 1.5)

ax.set_title('Validacion: Genes Conocidos de Cancer de Pancreas\n'
             'Verde = Gen clave | Amarillo = Hub | Naranja = DEG',
             fontsize=14, fontweight='bold', pad=20)

plt.savefig(os.path.join(RESULTS_DIR, "fig_18_validation_known_genes.png"))
plt.savefig(os.path.join(BASE_DIR, "07_validacion", "figures" if os.path.exists(os.path.join(BASE_DIR, "07_validacion", "figures")) else "results", "05_validation_known_genes.png"))
plt.close()
print("[OK] 05_validation_known_genes.png")

# --- figura 6: resumen grafico del pipeline ---
print("[6/6] Resumen grafico del pipeline...")

fig, axes = plt.subplots(2, 2, figsize=(16, 12))

# Panel A: Distribucion log2FC
ax = axes[0, 0]
up_genes = key_genes[key_genes['direction'] == 'UP']
down_genes = key_genes[key_genes['direction'] == 'DOWN']

ax.hist(up_genes['log2FoldChange'], bins=20, alpha=0.7, color='#FF6B6B', label=f'UP (n={len(up_genes)})')
ax.hist(down_genes['log2FoldChange'], bins=20, alpha=0.7, color='#4ECDC4', label=f'DOWN (n={len(down_genes)})')
ax.axvline(x=0, color='black', linestyle='-', linewidth=1)
ax.set_xlabel('log2 Fold Change')
ax.set_ylabel('Numero de genes')
ax.set_title('A) Distribucion de Fold Change\nen genes clave', fontweight='bold')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

# Panel B: kME por modulo
ax = axes[0, 1]
kme_data = []
kme_labels = []
for mod_color, me_name in module_me_map.items():
    mod_kme = key_genes[key_genes[module_col] == mod_color]['kME']
    if len(mod_kme) > 0:
        kme_data.append(mod_kme)
        kme_labels.append(f'{me_name} ({mod_color})\nn={len(mod_kme)}')

bp = ax.boxplot(kme_data,
                labels=kme_labels,
                patch_artist=True,
                widths=0.5)

for i, mod_color in enumerate(list(module_me_map.keys())[:len(kme_data)]):
    hex_color = color_palette.get(mod_color, '#999999')
    bp['boxes'][i].set_facecolor(hex_color + '40')  # semitransparente
    bp['boxes'][i].set_edgecolor(hex_color)

ax.axhline(y=0.5, color='red', linestyle='--', alpha=0.5, label='Umbral kME=0.5')
ax.set_ylabel('Module Membership (kME)')
ax.set_title('B) Distribucion de kME\npor modulo significativo', fontweight='bold')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Panel C: Top 15 genes por score combinado
ax = axes[1, 0]
key_genes_scored = key_genes.copy()
key_genes_scored['score'] = key_genes_scored['kME'].abs() * key_genes_scored['GS_death'].abs() * key_genes_scored['log2FoldChange'].abs()
key_genes_scored = key_genes_scored.sort_values('score', ascending=False)

top15 = key_genes_scored.head(15)
names = top15.apply(
    lambda row: row['gene_symbol'] if pd.notna(row['gene_symbol']) and row['gene_symbol'] != 'NA'
    else row['ensembl_id'][:12], axis=1
)
colors = [color_palette.get(c, '#999999') for c in top15[module_col]]

bars = ax.barh(range(len(top15)), top15['score'], color=colors, alpha=0.8, edgecolor='black', linewidth=0.5)
ax.set_yticks(range(len(top15)))
ax.set_yticklabels(names, fontsize=9, fontweight='bold')
ax.invert_yaxis()
ax.set_xlabel('Score combinado (|kME| x |GS| x |log2FC|)')
ax.set_title('C) Top 15 genes clave\n(score combinado)', fontweight='bold')
ax.grid(True, alpha=0.3, axis='x')

# Panel D: Funnel de seleccion
ax = axes[1, 1]

# valores dinamicos del funnel
n_degs_total = len(degs_all[(degs_all['padj'] < 0.05) & (degs_all['log2FoldChange'].abs() > 1)])
n_wgcna = len(module_assign)
n_hub = len(hub_genes)
n_key = len(key_genes)

categories = ['Total DEGs\n(|log2FC|>1)', 'Genes en\nWGCNA', 'Hub genes\n(kME>0.5, |GS|>0.15)',
              'Genes clave\n(DEG inter Hub)']
values = [n_degs_total, n_wgcna, n_hub, n_key]
colors_bar = ['#4472C4', '#ED7D31', '#A5A5A5', '#FF6B6B']

bars = ax.bar(categories, values, color=colors_bar, edgecolor='black', linewidth=0.5, alpha=0.8)

for bar, val in zip(bars, values):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 200,
            f'{val:,}', ha='center', va='bottom', fontweight='bold', fontsize=11)

ax.set_ylabel('Numero de genes')
ax.set_title('D) Pipeline de filtrado\n(funnel de seleccion)', fontweight='bold')
ax.set_yscale('log')
ax.grid(True, alpha=0.3, axis='y')
ax.set_ylim(50, 100000)

plt.suptitle('Resumen del Analisis de Integracion - TFM Cancer de Pancreas',
             fontsize=16, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(RESULTS_DIR, "fig_00_resumen_pipeline.png"))
plt.savefig(os.path.join(INTEGRATION_DIR, "figures", "06_resumen_pipeline.png"))
plt.close()
print("[OK] 06_resumen_pipeline.png")

# --- resumen ---
print(f"[OK] 6 figuras generadas en {RESULTS_DIR}/")
print("[OK] Visualizaciones finalizadas")
