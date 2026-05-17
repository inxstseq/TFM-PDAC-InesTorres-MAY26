#!/usr/bin/env python3
"""
15_machine_learning.py
Autora: Ines Torres
Fecha: 19/02/2026
Firma pronostica LASSO-Cox para genes clave de PDAC.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
import os
import warnings
warnings.filterwarnings('ignore')

# --- rutas ---
BASE = os.path.expanduser('~/TFM_PDAC')
VST_PATH = f'{BASE}/02_preprocesamiento/outputs_FINAL/matriz_vst_normalizada.txt'
META_PATH = f'{BASE}/02_preprocesamiento/outputs_FINAL/metadata_muestras.txt'
GENES_PATH = f'{BASE}/06_integracion/results/genes_clave_DEG_intersect_Hub.csv'
CLINICAL_PATH = f'{BASE}/01_datos/clinical/TCGA_PAAD_traits_numeric.txt'
FIG_DIR = f'{BASE}/08_figuras/07_machine_learning'
TABLE_DIR = f'{BASE}/11_tablas_resultados/tables/machine_learning'
RESULTS_DIR = f'{BASE}/resultados_TFM_v2'

os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(TABLE_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

np.random.seed(42)

# --- cargar datos ---
print("[OK] Cargando datos...")

vst = pd.read_csv(VST_PATH, sep='\t', index_col=0)
meta = pd.read_csv(META_PATH, sep='\t', index_col=0)
genes_df = pd.read_csv(GENES_PATH)
key_genes = genes_df['ensembl_id'].tolist()
clinical = pd.read_csv(CLINICAL_PATH, sep='\t', index_col=0)

# Quitar version de IDs Ensembl
vst.index = [x.split('.')[0] for x in vst.index]

# Filtrar a genes clave presentes en VST
key_genes_present = [g for g in key_genes if g in vst.index]
print(f"  VST: {vst.shape[0]} genes x {vst.shape[1]} muestras")
print(f"  Genes clave en VST: {len(key_genes_present)}")

# Matriz de features (muestras x genes)
X_full = vst.loc[key_genes_present].T

# Etiquetas (0=Normal, 1=Tumor)
common_samples = X_full.index.intersection(meta.index)
X_full = X_full.loc[common_samples]
meta_sub = meta.loc[common_samples]
y_full = (meta_sub['grupo'] == 'TCGA_Tumor').astype(int).values
X_full_vals = X_full.values.astype(np.float64)

print(f"  Muestras: {len(y_full)} (Tumor: {sum(y_full)}, Normal: {len(y_full)-sum(y_full)})")
gene_names = X_full.columns.tolist()

# Mapeo de simbolos
gene_symbol_map = dict(zip(genes_df['ensembl_id'], genes_df['gene_symbol']))

# Fallback para genes sin simbolo HGNC oficial en la anotacion
symbol_fallback = {
    'ENSG00000000005': 'TNMD',
    'ENSG00000237940': 'LINC01238',
    'ENSG00000274712': 'ENSG00000274712',   # pseudogen BPTF, sin simbolo HGNC
    'ENSG00000260306': 'ENSG00000260306',   # novel transcript, sin simbolo HGNC
    'ENSG00000249993': 'ENSG00000249993',   # lncRNA sin anotacion oficial
}
for ens_id, sym in symbol_fallback.items():
    if ens_id in gene_symbol_map and (pd.isna(gene_symbol_map[ens_id]) or str(gene_symbol_map[ens_id]) in ('NA', 'nan', '')):
        gene_symbol_map[ens_id] = sym

# --- funciones auxiliares ---

def standardize(X_train, X_test):
    """Estandariza (z-score) usando media/std del train."""
    mu = X_train.mean(axis=0)
    sigma = X_train.std(axis=0) + 1e-10
    return (X_train - mu) / sigma, (X_test - mu) / sigma, mu, sigma

def stratified_kfold(y, n_folds=10):
    """Genera indices de k-fold estratificado."""
    idx_0 = np.where(y == 0)[0]
    idx_1 = np.where(y == 1)[0]
    np.random.shuffle(idx_0)
    np.random.shuffle(idx_1)

    folds = [[] for _ in range(n_folds)]
    for i, idx in enumerate(idx_0):
        folds[i % n_folds].append(idx)
    for i, idx in enumerate(idx_1):
        folds[i % n_folds].append(idx)

    for f_idx in range(n_folds):
        test_idx = np.array(folds[f_idx])
        train_idx = np.array([i for j in range(n_folds) if j != f_idx for i in folds[j]])
        yield train_idx, test_idx

def sigmoid(z):
    z = np.clip(z, -500, 500)
    return 1.0 / (1.0 + np.exp(-z))

def concordance_index(times, events, risk_scores):
    """C-index: proporcion de pares concordantes."""
    n = len(times)
    concordant = 0
    discordant = 0
    tied = 0

    for i in range(n):
        if events[i] == 0:
            continue
        for j in range(n):
            if i == j:
                continue
            if times[j] > times[i]:
                if risk_scores[i] > risk_scores[j]:
                    concordant += 1
                elif risk_scores[i] < risk_scores[j]:
                    discordant += 1
                else:
                    tied += 1

    total = concordant + discordant + tied
    if total == 0:
        return 0.5
    return (concordant + 0.5 * tied) / total

def bootstrap_cindex(times, events, risk_scores, n_boot=1000):
    """Bootstrap para IC 95% del C-index."""
    n = len(times)
    c_indices = []

    for _ in range(n_boot):
        idx = np.random.choice(n, size=n, replace=True)
        t_b = times[idx]
        e_b = events[idx]
        r_b = risk_scores[idx]

        if sum(e_b) < 2:
            continue
        c_indices.append(concordance_index(t_b, e_b, r_b))

    c_indices = np.array(c_indices)
    ci_lower = np.percentile(c_indices, 2.5)
    ci_upper = np.percentile(c_indices, 97.5)
    return np.mean(c_indices), ci_lower, ci_upper

# --- firma pronostica LASSO-Cox ---
print("\n[OK] Firma genica pronostica (LASSO-Cox)")

# Datos de supervivencia (solo TCGA)
clinical_surv = clinical[['survival', 'death']].copy()
clinical_surv = clinical_surv.dropna()

tcga_samples = meta[meta['grupo'] == 'TCGA_Tumor'].index.tolist()

common_surv = list(set(X_full.index) & set(clinical_surv.index) & set(tcga_samples))
print(f"  Muestras TCGA con supervivencia: {len(common_surv)}")

X_surv = X_full.loc[common_surv].values.astype(np.float64)
times = clinical_surv.loc[common_surv, 'survival'].values.astype(np.float64)
events = clinical_surv.loc[common_surv, 'death'].values.astype(np.int32)
surv_gene_names = X_full.columns.tolist()

# Quitar muestras con tiempo <= 0
valid_mask = times > 0
X_surv = X_surv[valid_mask]
times = times[valid_mask]
events = events[valid_mask]
print(f"  Muestras validas (t > 0): {len(times)}, eventos: {sum(events)}")

# --- LASSO-Cox ---

class LassoCox:
    """Cox PH con penalizacion LASSO (L1), descenso por coordenadas."""

    def __init__(self, alpha=0.1, max_iter=200, tol=1e-6):
        self.alpha = alpha
        self.max_iter = max_iter
        self.tol = tol
        self.coef_ = None

    def _partial_likelihood_gradient(self, X, times, events, beta):
        """Gradiente de la log-verosimilitud parcial negativa."""
        n, p = X.shape
        risk_scores = X @ beta

        order = np.argsort(-times)
        X_sorted = X[order]
        times_sorted = times[order]
        events_sorted = events[order]
        risk_scores_sorted = risk_scores[order]

        gradient = np.zeros(p)

        exp_scores = np.exp(risk_scores_sorted - np.max(risk_scores_sorted))

        cum_exp = np.cumsum(exp_scores[::-1])[::-1]
        cum_exp_x = np.zeros((n, p))
        for j in range(p):
            cum_exp_x[:, j] = np.cumsum((exp_scores * X_sorted[:, j])[::-1])[::-1]

        for i in range(n):
            if events_sorted[i] == 1:
                gradient -= X_sorted[i]
                if cum_exp[i] > 0:
                    gradient += cum_exp_x[i] / cum_exp[i]

        return gradient / n

    def _soft_threshold(self, z, lam):
        """Soft thresholding para penalizacion L1."""
        if z > lam:
            return z - lam
        elif z < -lam:
            return z + lam
        return 0.0

    def fit(self, X, times, events):
        n, p = X.shape
        self.coef_ = np.zeros(p)

        for iteration in range(self.max_iter):
            old_coef = self.coef_.copy()

            grad = self._partial_likelihood_gradient(X, times, events, self.coef_)

            step_size = 0.1 / (iteration + 1) ** 0.5

            for j in range(p):
                z = self.coef_[j] - step_size * grad[j]
                self.coef_[j] = self._soft_threshold(z, self.alpha * step_size)

            if np.max(np.abs(self.coef_ - old_coef)) < self.tol:
                break

        return self

    def predict_risk(self, X):
        """Risk score (predictor lineal)."""
        return X @ self.coef_

    def get_selected_genes(self):
        """Indices de coeficientes no nulos."""
        return np.where(np.abs(self.coef_) > 1e-8)[0]


# --- buscar lambda optimo con 5-fold CV ---
# Estandarizacion DENTRO de cada fold para evitar data leakage
print("  Buscando lambda optimo con 5-fold CV...")

alphas = [0.5, 0.3, 0.2, 0.15, 0.1, 0.08, 0.05, 0.03, 0.02, 0.01]

alpha_results = []

for alpha in alphas:
    c_indices = []
    n_genes_selected = []

    n = len(times)
    perm = np.random.permutation(n)
    fold_size = n // 5

    for k in range(5):
        test_mask = np.zeros(n, dtype=bool)
        start = k * fold_size
        end = start + fold_size if k < 4 else n
        test_mask[perm[start:end]] = True
        train_mask = ~test_mask

        # Estandarizar usando solo datos de train (evitar data leakage)
        X_tr_raw = X_surv[train_mask]
        X_te_raw = X_surv[test_mask]
        tr_mu = X_tr_raw.mean(axis=0)
        tr_std = X_tr_raw.std(axis=0) + 1e-10
        X_tr = (X_tr_raw - tr_mu) / tr_std
        X_te = (X_te_raw - tr_mu) / tr_std

        t_tr = times[train_mask]
        t_te = times[test_mask]
        e_tr = events[train_mask]
        e_te = events[test_mask]

        model = LassoCox(alpha=alpha, max_iter=200)
        model.fit(X_tr, t_tr, e_tr)

        risk = model.predict_risk(X_te)
        c_idx = concordance_index(t_te, e_te, risk)
        c_indices.append(c_idx)
        n_genes_selected.append(len(model.get_selected_genes()))

    mean_c = np.mean(c_indices)
    std_c = np.std(c_indices)
    mean_genes = np.mean(n_genes_selected)

    alpha_results.append({
        'alpha': alpha,
        'C_index_mean': mean_c,
        'C_index_std': std_c,
        'N_genes_mean': mean_genes
    })
    print(f"    alpha={alpha:.3f}: C-index={mean_c:.4f}+/-{std_c:.4f}, genes={mean_genes:.0f}")

alpha_df = pd.DataFrame(alpha_results)
alpha_df.to_csv(f'{TABLE_DIR}/lasso_cox_alpha_selection.csv', index=False)

best_alpha_row = alpha_df.loc[alpha_df['C_index_mean'].idxmax()]
best_alpha = best_alpha_row['alpha']
print(f"[OK] Mejor alpha: {best_alpha} (C-index={best_alpha_row['C_index_mean']:.4f})")

# --- modelo LASSO-Cox final ---
print(f"  Entrenando modelo final (alpha={best_alpha})...")

# Estandarizar todo el dataset para el modelo final
X_surv_mu = X_surv.mean(axis=0)
X_surv_std = X_surv.std(axis=0) + 1e-10
X_surv_s = (X_surv - X_surv_mu) / X_surv_std

final_model = LassoCox(alpha=best_alpha, max_iter=500)
final_model.fit(X_surv_s, times, events)

selected_idx = final_model.get_selected_genes()
selected_genes_ensembl = [surv_gene_names[i] for i in selected_idx]
selected_genes_symbols = [gene_symbol_map.get(g, g) for g in selected_genes_ensembl]
selected_coefs = final_model.coef_[selected_idx]

print(f"  Genes seleccionados: {len(selected_idx)}")

signature_rows = []
for i in np.argsort(-np.abs(selected_coefs)):
    symbol = selected_genes_symbols[i]
    coef = selected_coefs[i]
    effect = 'Riesgo (+)' if coef > 0 else 'Protector (-)'
    print(f"    {symbol:<15} {coef:>10.4f}   {effect}")
    signature_rows.append({
        'Ensembl_ID': selected_genes_ensembl[i],
        'Gene_Symbol': symbol,
        'LASSO_Coef': coef,
        'Effect': effect,
        'HR': np.exp(coef)
    })

signature_df = pd.DataFrame(signature_rows)
signature_df.to_csv(f'{TABLE_DIR}/lasso_cox_firma_genica.csv', index=False)
signature_df.to_csv(f'{RESULTS_DIR}/tabla_lasso_cox_firma_genica.csv', index=False)
print("[OK] lasso_cox_firma_genica.csv")

# Guardar firma en formato para script 12
signature_for_export = pd.DataFrame({
    'gene_symbol': [r['Gene_Symbol'] for r in signature_rows],
    'lasso_coef': [r['LASSO_Coef'] for r in signature_rows]
})
signature_for_export.to_csv(f'{TABLE_DIR}/lasso_cox_signature.csv', index=False)
print("[OK] lasso_cox_signature.csv")

# --- risk score y estratificacion ---
risk_scores = final_model.predict_risk(X_surv_s)
median_risk = np.median(risk_scores)
risk_groups = np.where(risk_scores >= median_risk, 'High', 'Low')

print(f"  Alto riesgo: {sum(risk_groups == 'High')}, Bajo riesgo: {sum(risk_groups == 'Low')}")

c_index_full = concordance_index(times, events, risk_scores)
print(f"  C-index (modelo completo): {c_index_full:.4f}")

# Bootstrap CI para C-index
c_boot_mean, c_boot_lower, c_boot_upper = bootstrap_cindex(times, events, risk_scores, n_boot=1000)
print(f"  C-index bootstrap: {c_boot_mean:.4f} (IC 95%: {c_boot_lower:.4f}-{c_boot_upper:.4f})")

# --- Kaplan-Meier y log-rank ---

def kaplan_meier(times, events):
    """Curva de Kaplan-Meier."""
    unique_times = np.sort(np.unique(times[events == 1]))
    surv = 1.0
    km_times = [0.0]
    km_surv = [1.0]

    for t in unique_times:
        at_risk = np.sum(times >= t)
        died = np.sum((times == t) & (events == 1))
        if at_risk > 0:
            surv *= (1 - died / at_risk)
        km_times.append(t)
        km_surv.append(surv)

    return np.array(km_times), np.array(km_surv)

def log_rank_test(times1, events1, times2, events2):
    """Test log-rank."""
    all_times = np.sort(np.unique(np.concatenate([
        times1[events1 == 1], times2[events2 == 1]
    ])))

    O1 = 0
    E1 = 0
    V = 0

    for t in all_times:
        d1 = np.sum((times1 == t) & (events1 == 1))
        d2 = np.sum((times2 == t) & (events2 == 1))
        n1 = np.sum(times1 >= t)
        n2 = np.sum(times2 >= t)

        d = d1 + d2
        n = n1 + n2

        if n > 0:
            e1 = d * n1 / n
            O1 += d1
            E1 += e1

            if n > 1:
                V += d * n1 * n2 * (n - d) / (n**2 * (n - 1))

    if V > 0:
        chi2 = (O1 - E1)**2 / V
    else:
        chi2 = 0

    def chi2_cdf_1df(x):
        if x <= 0:
            return 0
        return math.erf(math.sqrt(x / 2))

    p_value = 1 - chi2_cdf_1df(chi2)
    return chi2, p_value

high_mask = risk_groups == 'High'
low_mask = risk_groups == 'Low'

km_high_t, km_high_s = kaplan_meier(times[high_mask], events[high_mask])
km_low_t, km_low_s = kaplan_meier(times[low_mask], events[low_mask])

lr_chi2, lr_pval = log_rank_test(
    times[high_mask], events[high_mask],
    times[low_mask], events[low_mask]
)

print(f"  Log-rank: chi2={lr_chi2:.4f}, p={lr_pval:.2e}")

# Mediana de supervivencia por grupo
def find_median_survival(km_t, km_s):
    for i in range(len(km_s)):
        if km_s[i] <= 0.5:
            return km_t[i]
    return np.nan

med_high = find_median_survival(km_high_t, km_high_s)
med_low = find_median_survival(km_low_t, km_low_s)
print(f"  Mediana alto riesgo: {med_high:.1f}m, bajo riesgo: {med_low:.1f}m")

# --- figura 01: KM high vs low risk ---
fig, ax = plt.subplots(figsize=(9, 7))
ax.step(km_high_t, km_high_s, where='post', color='#E53935', linewidth=2.5,
        label=f'High Risk (n={sum(high_mask)}, median={med_high:.1f}m)')
ax.step(km_low_t, km_low_s, where='post', color='#1E88E5', linewidth=2.5,
        label=f'Low Risk (n={sum(low_mask)}, median={med_low:.1f}m)')

for mask, color, km_t, km_s in [
    (high_mask, '#E53935', km_high_t, km_high_s),
    (low_mask, '#1E88E5', km_low_t, km_low_s)
]:
    censored_times = times[mask][events[mask] == 0]
    for ct in censored_times:
        idx = np.searchsorted(km_t, ct, side='right') - 1
        if idx >= 0 and idx < len(km_s):
            ax.plot(ct, km_s[idx], '|', color=color, markersize=8, alpha=0.5)

ax.set_xlabel('Time (months)', fontsize=12)
ax.set_ylabel('Survival Probability', fontsize=12)
n_sig_genes = len(selected_idx)
ax.set_title(f'LASSO-Cox Prognostic Signature ({n_sig_genes}-Gene Model)\n'
             f'Log-rank p = {lr_pval:.2e} | C-index = {c_index_full:.4f}', fontsize=13)
ax.legend(loc='upper right', fontsize=11)
ax.set_ylim(-0.05, 1.05)
ax.set_xlim(0, max(times) + 2)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f'{FIG_DIR}/01_KM_risk_groups_LASSO.png', dpi=300, bbox_inches='tight')
plt.savefig(f'{RESULTS_DIR}/fig_21_KM_risk_groups_LASSO.png', dpi=300, bbox_inches='tight')
plt.close()
print("[OK] 01_KM_risk_groups_LASSO.png")

# --- figura 02: alpha selection ---
fig, ax1 = plt.subplots(figsize=(9, 6))

ax1.errorbar(alpha_df['alpha'], alpha_df['C_index_mean'],
            yerr=alpha_df['C_index_std'], color='#1E88E5', marker='o',
            linewidth=2, capsize=4, label='C-index')
ax1.set_xlabel('LASSO Penalty (alpha)', fontsize=12)
ax1.set_ylabel('C-index (5-fold CV)', fontsize=12, color='#1E88E5')
ax1.tick_params(axis='y', labelcolor='#1E88E5')
ax1.axvline(x=best_alpha, color='red', linestyle='--', alpha=0.7, label=f'Best alpha = {best_alpha}')

ax2 = ax1.twinx()
ax2.plot(alpha_df['alpha'], alpha_df['N_genes_mean'], color='#4CAF50',
         marker='s', linewidth=2, label='N genes')
ax2.set_ylabel('Number of Selected Genes', fontsize=12, color='#4CAF50')
ax2.tick_params(axis='y', labelcolor='#4CAF50')

ax1.set_title('LASSO-Cox: Lambda Selection via Cross-Validation', fontsize=13)
ax1.legend(loc='upper left', fontsize=10)
ax2.legend(loc='upper right', fontsize=10)
ax1.set_xscale('log')
ax1.invert_xaxis()
ax1.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(f'{FIG_DIR}/02_lasso_alpha_selection.png', dpi=300, bbox_inches='tight')
plt.savefig(f'{RESULTS_DIR}/fig_22_lasso_alpha_selection.png', dpi=300, bbox_inches='tight')
plt.close()
print("[OK] 02_lasso_alpha_selection.png")

# --- figura 03: LASSO coefficients ---
if len(selected_idx) > 0:
    fig, ax = plt.subplots(figsize=(10, max(5, len(selected_idx) * 0.4)))

    sort_idx = np.argsort(np.abs(selected_coefs))
    sorted_symbols = [selected_genes_symbols[i] for i in sort_idx]
    sorted_coefs = selected_coefs[sort_idx]

    colors_coef = ['#E53935' if c > 0 else '#1E88E5' for c in sorted_coefs]

    y_pos = np.arange(len(sorted_symbols))
    ax.barh(y_pos, sorted_coefs, color=colors_coef, alpha=0.8, edgecolor='gray')
    ax.set_yticks(y_pos)
    ax.set_yticklabels(sorted_symbols, fontsize=10)
    ax.axvline(x=0, color='black', linewidth=0.8)
    ax.set_xlabel('LASSO-Cox Coefficient', fontsize=12)
    ax.set_title(f'LASSO-Cox Gene Signature Coefficients\n'
                 f'(alpha={best_alpha}, {len(selected_idx)} genes selected)', fontsize=13)

    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='#E53935', label='Risk Factor (beta > 0)'),
                       Patch(facecolor='#1E88E5', label='Protective (beta < 0)')]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=10)
    ax.grid(True, alpha=0.3, axis='x')

    plt.tight_layout()
    plt.savefig(f'{FIG_DIR}/03_lasso_coefficients.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{RESULTS_DIR}/fig_23_lasso_coefficients.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("[OK] 03_lasso_coefficients.png")

# --- figura 04: distribucion de risk scores ---
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

ax = axes[0]
ax.hist(risk_scores[high_mask], bins=20, alpha=0.7, color='#E53935', label='High Risk', edgecolor='gray')
ax.hist(risk_scores[low_mask], bins=20, alpha=0.7, color='#1E88E5', label='Low Risk', edgecolor='gray')
ax.axvline(x=median_risk, color='black', linestyle='--', linewidth=1.5, label=f'Median = {median_risk:.3f}')
ax.set_xlabel('Risk Score', fontsize=11)
ax.set_ylabel('Count', fontsize=11)
ax.set_title('Risk Score Distribution', fontsize=12)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

ax = axes[1]
sorted_risk_idx = np.argsort(risk_scores)
sorted_risks = risk_scores[sorted_risk_idx]
sorted_events_plot = events[sorted_risk_idx]
colors_waterfall = ['#E53935' if r >= median_risk else '#1E88E5' for r in sorted_risks]

ax.bar(range(len(sorted_risks)), sorted_risks, color=colors_waterfall, width=1.0, edgecolor='none')
ax.axhline(y=median_risk, color='black', linestyle='--', linewidth=1, alpha=0.7)
ax.set_xlabel('Patients (ranked)', fontsize=11)
ax.set_ylabel('Risk Score', fontsize=11)
ax.set_title('Waterfall Plot - Risk Score per Patient', fontsize=12)
ax.grid(True, alpha=0.3, axis='y')

plt.suptitle('LASSO-Cox Risk Score Analysis', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(f'{FIG_DIR}/04_risk_score_distribution.png', dpi=300, bbox_inches='tight')
plt.savefig(f'{RESULTS_DIR}/fig_24_risk_score_distribution.png', dpi=300, bbox_inches='tight')
plt.close()
print("[OK] 04_risk_score_distribution.png")

# --- guardar resultados finales ---

# Risk scores por paciente
risk_df = pd.DataFrame({
    'Sample': np.array(common_surv)[valid_mask],
    'Risk_Score': risk_scores,
    'Risk_Group': risk_groups,
    'Survival_Months': times,
    'Event': events
})
risk_df = risk_df.sort_values('Risk_Score', ascending=False)
risk_df.to_csv(f'{TABLE_DIR}/risk_scores_pacientes.csv', index=False)

# Resumen LASSO-Cox (incluye bootstrap CI)
lasso_summary = pd.DataFrame({
    'Parameter': ['Best_alpha', 'N_genes_selected', 'C_index_full',
                   'C_index_bootstrap_mean', 'C_index_CI95_lower', 'C_index_CI95_upper',
                   'LogRank_chi2', 'LogRank_p', 'Median_High_Risk',
                   'Median_Low_Risk', 'N_High_Risk', 'N_Low_Risk'],
    'Value': [best_alpha, len(selected_idx), c_index_full,
              c_boot_mean, c_boot_lower, c_boot_upper,
              lr_chi2, lr_pval, med_high, med_low,
              sum(high_mask), sum(low_mask)]
})
lasso_summary.to_csv(f'{TABLE_DIR}/lasso_cox_resumen.csv', index=False)
lasso_summary.to_csv(f'{RESULTS_DIR}/tabla_lasso_cox_resumen.csv', index=False)

print("[OK] Tablas finales guardadas")

# --- resumen ---
print(f"\n[OK] LASSO-Cox: {len(selected_idx)} genes, C-index={c_index_full:.4f} "
      f"(IC 95%: {c_boot_lower:.4f}-{c_boot_upper:.4f}), log-rank p={lr_pval:.2e}")
print("[OK] Analisis de firma pronostica finalizado")
