import logging
from pathlib import Path
import pandas as pd
import anndata as ad

logger = logging.getLogger(__name__)

class ClinicalAnalyticsPipeline:
    """
    Handles Phase 6 of the pipeline: Programmatic cell-type annotation tracking,
    compositional cross-tabulations, severity stratification, and dashboard dataframe extraction.
    """
    def __init__(self, config: dict):
        # Establish explicit biological annotation directory matching Phase 5 validation profiles
        self.cluster_annotations = {
            '0': 'Homeostatic Microglia',
            '1': 'Homeostatic Microglia',
            '2': 'Pre-inflammatory Microglia',
            '3': 'Stress-Response Microglia',
            '4': 'Infiltrating T-Cells',
            '5': 'Homeostatic Microglia',
            '6': 'Erythrocyte Contamination',
            '7': 'XIST+ Microglia Subset',
            '8': 'Disease-Associated Microglia (DAM)',
            '9': 'Infiltrating T-Cells',
            '10': 'Synaptic Debris/Astrocytes',
            '11': 'Activated Microglia',
            '12': 'Infiltrating Monocytes'
        }
        # Centralized targeting of verified biological biomarkers used in web dashboard
        self.target_markers = [
            "NEAT1", "XIST", "SORL1", "SPP1", "RGCC", "LPL", 
            "C1QC", "FTL", "CD3E", "CCL5", "CD2", "HBB", "GFAP"
        ]

    def annotate_cell_populations(self, adata: ad.AnnData) -> ad.AnnData:
        """Maps computed categorical cluster integers into explicit functional names."""
        logger.info("Mapping biological lexicon definitions onto cell cluster indices...")
        adata.obs['cell_type'] = adata.obs['leiden'].map(self.cluster_annotations)
        logger.info("✅ Programmatic cell-type annotation mapping successful.")
        return adata

    def analyze_composition_and_stratify(self, adata: ad.AnnData):
        """Constructs donor-level relative abundance percentages to stratify patient groups."""
        logger.info("Generating multi-patient compositional cross-tabulation arrays...")
        count_table = pd.crosstab(adata.obs['donor_id'], adata.obs['cell_type'])
        percentage_table = count_table.div(count_table.sum(axis=1), axis=0) * 100
        
        # Execute clinical stratification metric outputs inside internal logs
        if 'Disease-Associated Microglia (DAM)' in percentage_table.columns:
            dam_burden = percentage_table['Disease-Associated Microglia (DAM)'].sort_values(ascending=False)
            logger.info("📋 CLINICAL SEVERITY PROFILE STATUS - RANKED BY DAM BURDEN (%):")
            for rank, (donor, value) in enumerate(dam_burden.items(), 1):
                logger.info(f"   Rank {rank}: Patient {donor} -> {value:.2f}% DAM Abundance Burden")
                
        return percentage_table

    def extract_dashboard_data_layer(self, adata: ad.AnnData, output_path: str = "data/processed/dashboard_data.csv"):
        """Extracts and serializes a high-performance 4MB matrix for the interactive Streamlit Cloud tier."""
        logger.info("📦 Initializing data extraction layer for web app deployment...")
        
        # Track 2D manifold coordinates
        dash_df = pd.DataFrame(
            adata.obsm["X_umap"], columns=["UMAP_1", "UMAP_2"], index=adata.obs_names
        )
        dash_df["Cell_Type"] = adata.obs["cell_type"].astype(str)
        dash_df["Donor_ID"] = adata.obs["donor_id"].astype(str)

        # Vectorize selected expression tracks
        for gene in self.target_markers:
            if gene in adata.var_names:
                dash_df[gene] = adata[:, gene].X.toarray().flatten()

        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        dash_df.to_csv(output_path)
        logger.info(f"🎉 Production web application layer safely generated at: {output_path}")

    def run(self, adata: ad.AnnData) -> ad.AnnData:
        """Master execution loop for Phase 6 clinical analytics."""
        adata = self.annotate_cell_populations(adata)
        _ = self.analyze_composition_and_stratify(adata)
        self.extract_dashboard_data_layer(adata)
        return adata
