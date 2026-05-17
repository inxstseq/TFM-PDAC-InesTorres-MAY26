#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
12_enriquecimiento_genes_clave.py
Autor: Ines Torres
Fecha: 17/02/2026
Enriquecimiento GO y KEGG de los genes hub y la firma pronostica LASSO-Cox.
Se ejecuta DESPUES de supervivencia y LASSO para interpretar biologicamente
los genes que han demostrado relevancia clinica.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy as gp
import os
import warnings
warnings.filterwarnings('ignore')

# --- configuracion ---

os.chdir(os.path.expanduser("~/TFM_PDAC"))

OUT_DIR = "06_integracion/enrichment"
FIG_DIR = "06_integracion/enrichment/figures"
RESULTS_DIR = "resultados_TFM_v2"
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

# --- cargo genes ---

# 1) genes clave (DEG intersect Hub) - los que pasaron doble filtro
key_genes = pd.read_csv("06_integracion/results/genes_clave_DEG_intersect_Hub.csv")

# 2) firma LASSO-Cox (los genes seleccionados por el modelo pronostico)
# estos se guardaron en el script 15
lasso_path = "11_tablas_resultados/tables/machine_learning/lasso_cox_signature.csv"
if os.path.exists(lasso_path):
    lasso_genes = pd.read_csv(lasso_path)
    lasso_symbols = lasso_genes['gene_symbol'].dropna().unique().tolist()
    print(f"[OK] Firma LASSO-Cox: {len(lasso_symbols)} genes")
else:
    lasso_symbols = []
    print("[WARN] No se encontro lasso_cox_signature.csv, usando solo hub genes")

# 3) genes con supervivencia significativa (KM p < 0.05)
km_path = "06_integracion/survival/kaplan_meier_summary.csv"
if os.path.exists(km_path):
    km_results = pd.read_csv(km_path)
    km_sig = km_results[km_results['LogRank_p'] < 0.05]
    km_symbols = km_sig['Gene'].dropna().unique().tolist()
    print(f"[OK] Genes KM significativos: {len(km_symbols)}")
else:
    km_symbols = []
    print("[WARN] No se encontro kaplan_meier_summary.csv")

print(f"[OK] Genes clave (DEG+Hub): {len(key_genes)}")

# extraigo simbolos por modulo (dinamico, no hardcodeado)
module_col = 'module_color' if 'module_color' in key_genes.columns else 'Module'
modules_presentes = key_genes[module_col].unique()
modules_presentes = [m for m in modules_presentes if m != 'grey' and pd.notna(m)]

genes_por_modulo = {}
for mod in modules_presentes:
    mod_genes = key_genes[key_genes[module_col] == mod]
    symbols = mod_genes['gene_symbol'].dropna()
    symbols = symbols[symbols != 'NA'].unique().tolist()
    if len(symbols) > 0:
        genes_por_modulo[mod] = symbols

# conjunto final: hub genes + firma LASSO + KM significativos (union)
all_hub_symbols = key_genes['gene_symbol'].dropna()
all_hub_symbols = all_hub_symbols[all_hub_symbols != 'NA'].unique().tolist()

# genes con evidencia clinica: aparecen en KM significativo o firma LASSO
genes_clinicos = list(set(km_symbols + lasso_symbols))
genes_clinicos = [g for g in genes_clinicos if g and g != 'NA']

print(f"\n  Hub genes: {len(all_hub_symbols)}")
for mod, syms in genes_por_modulo.items():
    print(f"  {mod}: {len(syms)} genes")
print(f"  Genes con evidencia clinica: {len(genes_clinicos)}")
if lasso_symbols:
    print(f"  Firma LASSO: {', '.join(lasso_symbols)}")

# --- funcion de enriquecimiento ---

def run_enrichment(gene_list, gene_set_name, output_prefix, description):
    """Ejecuta enriquecimiento con gseapy/Enrichr."""
    print(f"\n  {description} | DB: {gene_set_name} | {len(gene_list)} genes")

    if len(gene_list) < 3:
        print(f"    [WARN] Muy pocos genes ({len(gene_list)}), se omite")
        return pd.DataFrame()

    try:
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=gene_set_name,
            organism='human',
            outdir=None,
            cutoff=0.05
        )

        results = enr.results

        if len(results) == 0:
            print(f"    [WARN] Sin resultados significativos")
            return pd.DataFrame()

        results_filtered = results[results['Adjusted P-value'] < 0.2].copy()
        print(f"    [OK] padj<0.05: {len(results)} | FDR<0.2: {len(results_filtered)}")

        results.to_csv(
            os.path.join(OUT_DIR, f"{output_prefix}_complete.txt"),
            sep="\t", index=False
        )

        if len(results_filtered) > 0:
            results_filtered.to_csv(
                os.path.join(OUT_DIR, f"{output_prefix}_filtered.txt"),
                sep="\t", index=False
            )

        return results_filtered

    except Exception as e:
        print(f"    [ERROR] {str(e)}")
        return pd.DataFrame()


def plot_top_terms(df, title, filename, top_n=20, color='#E64B35'):
    """Barplot horizontal de top terminos enriquecidos."""
    if len(df) == 0:
        return

    df_plot = df.head(top_n).copy()
    df_plot = df_plot.sort_values('Adjusted P-value', ascending=False)

    df_plot['Term_short'] = df_plot['Term'].apply(
        lambda x: x[:55] + '...' if len(x) > 55 else x
    )

    fig, ax = plt.subplots(figsize=(14, max(6, len(df_plot) * 0.45)))

    y_pos = np.arange(len(df_plot))
    ax.barh(
        y_pos,
        -np.log10(df_plot['Adjusted P-value']),
        color=color, alpha=0.8,
        edgecolor='black', linewidth=0.5
    )

    ax.set_yticks(y_pos)
    ax.set_yticklabels(df_plot['Term_short'], fontsize=10)
    ax.set_xlabel('-log10(Adjusted P-value)', fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=13, fontweight='bold', pad=15)
    ax.axvline(x=-np.log10(0.05), color='red', linestyle='--',
               linewidth=2, alpha=0.7, label='padj = 0.05')
    ax.grid(axis='x', alpha=0.3, linestyle='--')
    ax.legend(fontsize=10)

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, filename), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(RESULTS_DIR, f"fig_enrichment_{filename}"), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"    [OK] Figura: {filename}")


# --- enriquecimiento hub genes por modulo (dinamico) ---
# el objetivo es caracterizar funcionalmente los modulos que se asocian a mortalidad

print("\n# --- enriquecimiento por modulo (hub genes) ---")

# paleta de colores para los modulos
module_plot_colors = {
    'turquoise': '#1B9E77', 'blue': '#4472C4', 'brown': '#D95F02',
    'yellow': '#E6AB02', 'green': '#66A61E', 'red': '#E7298A',
    'black': '#333333', 'pink': '#F781BF', 'magenta': '#984EA3',
    'purple': '#9370DB', 'greenyellow': '#A6D854', 'tan': '#D2B48C',
    'salmon': '#FA8072', 'cyan': '#00BFC4', 'midnightblue': '#191970',
    'lightcyan': '#B0E0E6', 'grey60': '#999999', 'lightgreen': '#90EE90'
}

go_bp_results = {}
kegg_results = {}
fig_counter = 1

for mod_name, symbols in genes_por_modulo.items():
    print(f"\n--- MODULO {mod_name.upper()} ({len(symbols)} genes) ---")

    plot_color = module_plot_colors.get(mod_name, '#E64B35')

    go_bp = run_enrichment(symbols, 'GO_Biological_Process_2023',
                           f'{mod_name.upper()}_GO_BP', f'GO BP ({mod_name})')
    kegg = run_enrichment(symbols, 'KEGG_2021_Human',
                          f'{mod_name.upper()}_KEGG', f'KEGG ({mod_name})')

    go_bp_results[mod_name] = go_bp
    kegg_results[mod_name] = kegg

    if len(go_bp) > 0:
        plot_top_terms(go_bp,
                       f'GO BP - Modulo {mod_name}\nHub genes',
                       f'{fig_counter:02d}_GO_BP_{mod_name.upper()}.png',
                       top_n=15, color=plot_color)
        fig_counter += 1
    if len(kegg) > 0:
        plot_top_terms(kegg,
                       f'KEGG - Modulo {mod_name}\nHub genes',
                       f'{fig_counter:02d}_KEGG_{mod_name.upper()}.png',
                       top_n=15, color=plot_color)
        fig_counter += 1

# --- enriquecimiento de genes con evidencia clinica ---
# estos son los que salieron significativos en supervivencia y/o fueron seleccionados por LASSO

go_bp_clin = pd.DataFrame()
kegg_clin = pd.DataFrame()

if len(genes_clinicos) >= 3:
    print("\n--- GENES CON EVIDENCIA CLINICA (KM + LASSO) ---")

    go_bp_clin = run_enrichment(genes_clinicos, 'GO_Biological_Process_2023',
                                 'CLINICAL_GO_BP', 'GO BP (clinicos)')
    kegg_clin = run_enrichment(genes_clinicos, 'KEGG_2021_Human',
                                'CLINICAL_KEGG', 'KEGG (clinicos)')

    if len(go_bp_clin) > 0:
        plot_top_terms(go_bp_clin,
                       'GO BP - Genes con Evidencia Clinica\n(KM significativos + Firma LASSO-Cox)',
                       f'{fig_counter:02d}_GO_BP_CLINICAL.png', top_n=15, color='#2E8B57')
        fig_counter += 1
    if len(kegg_clin) > 0:
        plot_top_terms(kegg_clin,
                       'KEGG - Genes con Evidencia Clinica\n(KM significativos + Firma LASSO-Cox)',
                       f'{fig_counter:02d}_KEGG_CLINICAL.png', top_n=15, color='#2E8B57')
        fig_counter += 1
else:
    print("\n[WARN] Pocos genes clinicos, se omite enriquecimiento especifico")

# --- figura comparativa entre modulos (si hay al menos 2) ---

mods_con_go = {m: df for m, df in go_bp_results.items() if len(df) > 0}

if len(mods_con_go) >= 2:
    n_panels = min(len(mods_con_go), 4)  # maximo 4 paneles
    mod_names = list(mods_con_go.keys())[:n_panels]

    fig, axes = plt.subplots(1, n_panels, figsize=(10 * n_panels, 10))
    if n_panels == 1:
        axes = [axes]

    for idx, mod in enumerate(mod_names):
        top_terms = mods_con_go[mod].head(10).sort_values('Adjusted P-value', ascending=False)
        top_terms['Term_short'] = top_terms['Term'].apply(
            lambda x: x[:45] + '...' if len(x) > 45 else x)
        plot_color = module_plot_colors.get(mod, '#E64B35')

        axes[idx].barh(range(len(top_terms)),
                     -np.log10(top_terms['Adjusted P-value']),
                     color=plot_color, alpha=0.8, edgecolor='black', linewidth=0.5)
        axes[idx].set_yticks(range(len(top_terms)))
        axes[idx].set_yticklabels(top_terms['Term_short'], fontsize=10)
        axes[idx].set_xlabel('-log10(padj)', fontsize=12)
        axes[idx].set_title(f'{mod}',
                           fontsize=12, fontweight='bold')
        axes[idx].axvline(x=-np.log10(0.05), color='red', linestyle='--', alpha=0.5)
        axes[idx].grid(axis='x', alpha=0.3)

    plt.suptitle('Comparacion GO (BP): Modulos con Hub Genes\nDoble Filtro (DEG + WGCNA)',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, f'{fig_counter:02d}_comparison_modules_GO_BP.png'),
                dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(RESULTS_DIR, f'fig_enrichment_{fig_counter:02d}_comparison_modules_GO_BP.png'),
                dpi=300, bbox_inches='tight')
    plt.close()
    print(f"[OK] Figura comparativa: {fig_counter:02d}_comparison_modules_GO_BP.png")

# --- resumen ---

print(f"\n# --- resumen ---")
print(f"  Hub genes analizados: {len(all_hub_symbols)}")
for mod in genes_por_modulo:
    go_n = len(go_bp_results.get(mod, []))
    kegg_n = len(kegg_results.get(mod, []))
    print(f"  {mod}: GO_BP={go_n} KEGG={kegg_n}")
if len(genes_clinicos) >= 3:
    print(f"  Clinicos: GO_BP={len(go_bp_clin)} KEGG={len(kegg_clin)}")

# top terms del primer modulo con resultados
for mod, df in go_bp_results.items():
    if len(df) > 0:
        print(f"\n  Top 5 GO BP ({mod}):")
        for _, row in df.head(5).iterrows():
            print(f"    - {row['Term'][:60]} (padj={row['Adjusted P-value']:.2e})")
        break

print(f"\n  Archivos en: {OUT_DIR}/ | Figuras en: {FIG_DIR}/")
print("[OK] Enriquecimiento completado")
