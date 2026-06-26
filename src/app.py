import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import streamlit as st

# 1. PAGE CONFIGURATION
st.set_page_config(
    page_title="ALS Microglia Explorer", page_icon="🧬", layout="wide"
)

st.title("🧬 Mapping Microglial Heterogeneity in the ALS Spinal Cord")
st.markdown(
    "Demonstrating production-grade pipeline deployment for the Ajami Lab position."
)
st.markdown("---")


# 2. OPTIMIZED DATA LOADING
@st.cache_data
def load_data():
    return pd.read_csv("data/processed/dashboard_data.csv", index_col=0)


df = load_data()

# 3. SIDEBAR PATIENT FILTERS
st.sidebar.title("Dashboard Controls")
all_donors = ["All Patients"] + list(df["Donor_ID"].unique())
selected_donor = st.sidebar.selectbox("Filter Landscape by Patient:", all_donors)

# Dynamically slice the dataframe based on sidebar selection
if selected_donor != "All Patients":
    plot_df = df[df["Donor_ID"] == selected_donor]
else:
    plot_df = df

# 4. SPLIT SCREEN LAYOUT (Columns)
col_left, col_right = st.columns([1, 1])

# LEFT COLUMN: INTERACTIVE MANIFOLD MAP
with col_left:
    st.subheader("High-Dimensional Manifold Projection")

    gene_list = [
        c
        for c in df.columns
        if c not in ["UMAP_1", "UMAP_2", "Cell_Type", "Donor_ID"]
    ]
    selected_gene = st.selectbox("Select a target biomarker:", gene_list)

    fig1, ax1 = plt.subplots(figsize=(6, 4.5))
    points = ax1.scatter(
        x=plot_df["UMAP_1"],
        y=plot_df["UMAP_2"],
        c=plot_df[selected_gene],
        cmap="viridis",
        s=4,
        alpha=0.8,
    )
    ax1.set_axis_off()

    # Colorbar configuration
    cbar = fig1.colorbar(points, ax=ax1, shrink=0.7, pad=0.02)
    cbar.set_label(
        f"Expression Vector Intensity: {selected_gene}",
        fontsize=8,
        weight="bold",
    )
    cbar.ax.tick_params(labelsize=7)

    st.pyplot(fig1)

# RIGHT COLUMN: DYNAMIC COMPOSITIONAL MODELING
with col_right:
    st.subheader("Immune Ecosystem Abundance Profiles")
    st.write("")  # Whitespace buffer to align dropdown alignment

    # Programmatically calculate percentages on the fly based on the current data slice
    count_table = pd.crosstab(plot_df["Donor_ID"], plot_df["Cell_Type"])
    percentage_table = count_table.div(count_table.sum(axis=1), axis=0) * 100

    fig2, ax2 = plt.subplots(figsize=(6, 4))
    percentage_table.plot(kind="bar", stacked=True, ax=ax2, colormap="tab20")

    plt.title("Relative Abundance per Patient Cohort", fontsize=10, weight="bold")
    plt.ylabel("Relative Abundance (%)", fontsize=8)
    plt.xlabel("Patient Donor ID", fontsize=8)
    plt.xticks(rotation=0, fontsize=8)
    plt.yticks(fontsize=8)

    # Place legend cleanly outside the plot canvas
    plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7)
    plt.tight_layout()

    st.pyplot(fig2)