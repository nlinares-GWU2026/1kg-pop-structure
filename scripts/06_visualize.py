import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import numpy as np
import os

# ============================================================
# CONFIGURATION
# ============================================================

# Superpopulation colors — these are the standard colors used
# in most 1000 Genomes publications so readers recognize them
SUPERPOP_COLORS = {
    "AFR": "#E41A1C",  # red
    "EUR": "#377EB8",  # blue
    "EAS": "#4DAF4A",  # green
    "SAS": "#984EA3",  # purple
    "AMR": "#FF7F00"   # orange
}

# Output directory
os.makedirs("results/plots", exist_ok=True)

# ============================================================
# STEP 1 — LOAD DATA
# ============================================================

print("Loading data...")

# Load population panel — maps sample ID to population and superpopulation
panel = pd.read_csv(
    "data/integrated_call_samples_v3.20130502.ALL.panel",
    sep="\t",
    usecols=["sample", "pop", "super_pop"]
)
print(f"  Panel: {len(panel)} samples")

# Load PCA eigenvectors — PC scores for each individual
eigenvec = pd.read_csv(
    "results/pca/chr22_pca.eigenvec",
    sep=r"\s+",       # whitespace separated
    header=None       # no column names in file
)
# Name the columns — first two are IDs, rest are PC scores
eigenvec.columns = ["FID", "IID"] + [f"PC{i}" for i in range(1, 11)]
print(f"  Eigenvec: {len(eigenvec)} individuals, 10 PCs")

# Load PCA eigenvalues — variance explained per PC
eigenval = pd.read_csv(
    "results/pca/chr22_pca.eigenval",
    header=None,
    names=["eigenvalue"]
)
# Calculate percent variance explained
eigenval["pct_variance"] = (
    eigenval["eigenvalue"] / eigenval["eigenvalue"].sum() * 100
)
eigenval["PC"] = range(1, 11)
print(f"  Eigenvalues loaded, PC1 explains {eigenval['pct_variance'][0]:.1f}%")

# Load sNMF Q matrix for best K=5
# LEA saves Q files as: chr22.snmf/K5/run{n}/chr22_r{n}.5.Q
# We need the run with lowest cross-entropy — run the best one
# Check which run files exist
import glob
q_files = sorted(glob.glob(
    "results/admixture/chr22.snmf/K5/run*/chr22_r*.5.Q"
))
print(f"  Found {len(q_files)} Q files for K=5: {q_files}")

# Load the first run's Q file (we'll use run1 — LEA picks best internally)
q_matrix = pd.read_csv(
    q_files[0],
    sep=r"\s+",
    header=None
)
q_matrix.columns = [f"Anc{i+1}" for i in range(5)]
print(f"  Q matrix: {q_matrix.shape[0]} individuals x {q_matrix.shape[1]} ancestral populations")

# ============================================================
# STEP 2 — MERGE DATA
# ============================================================

# Merge eigenvec with panel to get population labels
pca_data = eigenvec.merge(
    panel,
    left_on="IID",
    right_on="sample",
    how="left"
)

# Add Q matrix to the merged data
# Both are ordered the same way (same individuals, same order)
for col in q_matrix.columns:
    pca_data[col] = q_matrix[col].values

print(f"\nMerged dataset: {len(pca_data)} individuals")
print(f"Superpopulations: {pca_data['super_pop'].value_counts().to_dict()}")

# ============================================================
# PLOT 1 — PCA SCATTER PLOT (PC1 vs PC2)
# ============================================================

print("\nGenerating PCA scatter plot...")

fig, ax = plt.subplots(figsize=(9, 7))

# Plot each superpopulation separately so they get individual legend entries
for spop, group in pca_data.groupby("super_pop"):
    ax.scatter(
        group["PC1"],
        group["PC2"],
        c=SUPERPOP_COLORS[spop],
        label=spop,
        alpha=0.6,       # slight transparency so overlapping points visible
        s=18,            # point size
        linewidths=0     # no edge lines on points
    )

# Axis labels include variance explained — standard in population genetics
pc1_var = eigenval.loc[eigenval["PC"] == 1, "pct_variance"].values[0]
pc2_var = eigenval.loc[eigenval["PC"] == 2, "pct_variance"].values[0]
ax.set_xlabel(f"PC1 ({pc1_var:.1f}% variance explained)", fontsize=12)
ax.set_ylabel(f"PC2 ({pc2_var:.1f}% variance explained)", fontsize=12)
ax.set_title(
    "Population Structure — 1000 Genomes Phase 3\nPCA (chr22, 7,751 SNPs)",
    fontsize=13,
    pad=12
)

# Legend
ax.legend(
    title="Superpopulation",
    title_fontsize=10,
    fontsize=9,
    framealpha=0.8,
    markerscale=1.5
)

ax.grid(True, alpha=0.2, linewidth=0.5)
plt.tight_layout()
plt.savefig("results/plots/pca_pc1_pc2.png", dpi=150, bbox_inches="tight")
plt.close()
print("  Saved: results/plots/pca_pc1_pc2.png")

# ============================================================
# PLOT 2 — SCREE PLOT
# ============================================================

print("Generating scree plot...")

fig, ax = plt.subplots(figsize=(7, 4))

ax.plot(
    eigenval["PC"],
    eigenval["pct_variance"],
    marker="o",
    color="#377EB8",
    linewidth=2,
    markersize=7,
    markerfacecolor="white",
    markeredgewidth=2
)

# Annotate PC1-PC2 dominance
ax.annotate(
    "Major structure\n(PC1+PC2 = 72.8%)",
    xy=(2, eigenval.loc[eigenval["PC"]==2, "pct_variance"].values[0]),
    xytext=(3.5, 20),
    fontsize=9,
    arrowprops=dict(arrowstyle="->", color="gray"),
    color="gray"
)

# Annotate the elbow — the drop between PC4 and PC5
ax.annotate(
    "Elbow\n(PC4→PC5)",
    xy=(5, eigenval.loc[eigenval["PC"]==5, "pct_variance"].values[0]),
    xytext=(6, 5),
    fontsize=9,
    arrowprops=dict(arrowstyle="->", color="gray"),
    color="gray"
)

# Annotate plateau beginning
ax.annotate(
    "Plateau\n(noise floor)",
    xy=(5, eigenval.loc[eigenval["PC"]==5, "pct_variance"].values[0]),
    xytext=(6.5, 8),
    fontsize=9,
    arrowprops=dict(arrowstyle="->", color="gray"),
    color="gray"
)

ax.set_xlabel("Principal Component", fontsize=12)
ax.set_ylabel("Variance Explained (%)", fontsize=12)
ax.set_title(
    "PCA Scree Plot — chr22\nVariance explained per PC",
    fontsize=13,
    pad=12
)
ax.set_xticks(range(1, 11))
ax.grid(True, alpha=0.3, linewidth=0.5)
plt.tight_layout()
plt.savefig("results/plots/pca_scree.png", dpi=150, bbox_inches="tight")
plt.close()
print("  Saved: results/plots/pca_scree.png")

# ============================================================
# PLOT 3 — ADMIXTURE BAR CHART (K=5)
# ============================================================

print("Generating ADMIXTURE bar chart...")

# Sort individuals by superpopulation then population
# so bars group cleanly by ancestry
pca_data_sorted = pca_data.sort_values(["super_pop", "pop"]).reset_index(drop=True)

anc_cols = [f"Anc{i+1}" for i in range(5)]

# Assign colors using the Hungarian algorithm for optimal 1-to-1 mapping
# between ancestral components and superpopulations.
# The greedy approach (take each superpop's dominant component) fails
# when two superpopulations share the same dominant component.
from scipy.optimize import linear_sum_assignment

SPOP_ORDER = ["AFR", "EUR", "EAS", "SAS", "AMR"]

# Build a 5x5 matrix of mean ancestry proportions
# rows = superpopulations, columns = ancestral components
mean_q = np.zeros((5, 5))
for i, spop in enumerate(SPOP_ORDER):
    subset = pca_data_sorted[
        pca_data_sorted["super_pop"] == spop
    ][anc_cols]
    mean_q[i] = subset.mean().values

print("  Mean Q matrix (superpops x components):")
for i, spop in enumerate(SPOP_ORDER):
    print(f"    {spop}: {[f'{v:.3f}' for v in mean_q[i]]}")

# Hungarian algorithm finds the optimal 1-to-1 assignment
# linear_sum_assignment minimizes cost, so we negate to maximize
row_ind, col_ind = linear_sum_assignment(-mean_q)
# col_ind[i] = which ancestral component belongs to superpop i

print("  Component assignments:")
anc_colors = ["#aaaaaa"] * 5  # default gray for unassigned
for i, spop in enumerate(SPOP_ORDER):
    component_idx = col_ind[i]
    anc_colors[component_idx] = SUPERPOP_COLORS[spop]
    print(f"    {spop} → Anc{component_idx+1} "
          f"(mean proportion: {mean_q[i, component_idx]:.3f})")

fig, ax = plt.subplots(figsize=(16, 3.5))

# Draw stacked bars — one per individual
bar_width = 1.0
bottom = np.zeros(len(pca_data_sorted))

for i, col in enumerate(anc_cols):
    ax.bar(
        range(len(pca_data_sorted)),
        pca_data_sorted[col].values,
        bottom=bottom,
        color=anc_colors[i],
        width=bar_width,
        linewidth=0      # no gaps between bars
    )
    bottom += pca_data_sorted[col].values

# Add superpopulation labels below the chart
spop_order = pca_data_sorted["super_pop"].values
boundaries = []
current = spop_order[0]
start = 0
for i, s in enumerate(spop_order):
    if s != current or i == len(spop_order) - 1:
        boundaries.append((start, i, current))
        start = i
        current = s

for start, end, spop in boundaries:
    mid = (start + end) / 2
    ax.text(
        mid, -0.15, spop,
        ha="center", va="top",
        fontsize=9,
        color=SUPERPOP_COLORS.get(spop, "black"),
        fontweight="bold",
        transform=ax.get_xaxis_transform()
    )
    # Draw vertical divider between superpopulations
    if end < len(pca_data_sorted) - 1:
        ax.axvline(x=end, color="white", linewidth=1.5)

ax.set_xlim(-0.5, len(pca_data_sorted) - 0.5)
ax.set_ylim(0, 1)
ax.set_xticks([])
ax.set_ylabel("Ancestry proportion", fontsize=11)
ax.set_title(
    "Ancestry Estimation — sNMF K=5\n1000 Genomes Phase 3 (chr22, 7,751 SNPs)",
    fontsize=13,
    pad=12
)

plt.tight_layout()
plt.savefig(
    "results/plots/admixture_K5.png",
    dpi=150,
    bbox_inches="tight"
)
plt.close()
print("  Saved: results/plots/admixture_K5.png")

# ============================================================
# DONE
# ============================================================

print("\nAll plots saved to results/plots/")
print("Files:")
for f in sorted(os.listdir("results/plots/")):
    print(f"  {f}")
