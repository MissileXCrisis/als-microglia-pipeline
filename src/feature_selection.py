import logging
from pathlib import Path
import anndata as ad
import pandas as pd
import scanpy as sc

logger = logging.getLogger(__name__)

class FeatureSelectionPipeline:
    """
    Handles Phase 3 of the pipeline: Size normalization, variance stabilization,
    Highly Variable Gene (HVG) masking, unit-variance scaling, and PCA loading analysis.
    """
    def __init__(self, config: dict):
        self.n_top_genes = config['feature_selection']['n_top_genes']
        self.n_pcs = config['feature_selection']['n_pcs']

    def stabilize_variance(self, adata: ad.AnnData) -> ad.AnnData:
        """Saves a raw copy of counts, log-normalizes matrix to uniform library depth."""
        logger.info("Freezing deep raw count layer into .raw matrix slot for downstream DE tests...")
        adata.raw = adata.copy()

        logger.info("Library-size normalizing cells to a uniform target sum of 10,000 counts...")
        sc.pp.normalize_total(adata, target_sum=1e4)

        logger.info("Stabilizing matrix variance using natural log-transformation [log(X + 1)]...")
        sc.pp.log1p(adata)
        return adata

    def extract_highly_variable_features(self, adata: ad.AnnData) -> ad.AnnData:
        """Isolates biological variance from technical background noise using a batch-aware mask."""
        logger.info(f"Isolating top {self.n_top_genes} High-Variance features using batch-aware Seurat strategy...")
        sc.pp.highly_variable_genes(
            adata, 
            n_top_genes=self.n_top_genes, 
            flavor='seurat', 
            batch_key='batch'
        )
        return adata

    def compute_pca_and_loadings(self, adata: ad.AnnData) -> ad.AnnData:
        """Scales gene variance, clips outliers, and calculates 50 orthogonal axes of variation."""
        logger.info("Scaling highly variable genes to unit-variance (clipping extreme outliers at 10)...")
        sc.pp.scale(adata, max_value=10)

        logger.info(f"Compressing data dimensions into {self.n_comps if hasattr(self, 'n_comps') else 50} PCA coordinates...")
        sc.tl.pca(adata, n_comps=self.n_pcs, use_highly_variable=True)
        
        # Log structural loadings validation info directly to developer terminal
        pca_loadings = pd.DataFrame(adata.varm['PCs'], index=adata.var_names)
        logger.info("✅ Mathematical dimension compression complete.")
        logger.info(f"Top 3 genes driving positive PC1 axis: {list(pca_loadings[0].sort_values(ascending=False).head(3).index)}")
        logger.info(f"Top 3 genes driving positive PC2 axis: {list(pca_loadings[1].sort_values(ascending=False).head(3).index)}")
        return adata

    def save_integrated_workspace(self, adata: ad.AnnData, output_path: str = "data/processed/als_spinal_cord_phase3_pca.h5ad"):
        """Saves the complete integrated state to disk using high-efficiency gzip compression."""
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        logger.info(f"Writing finalized compressed integrated workspace state to: {output_path}")
        
        # Save state with gzip compression to maintain a clean disk footprint
        adata.write_h5ad(output_path, compression="gzip")
        logger.info("🎉 Phase 3 Complete! Workspace serialized to disk and insulated against memory leaks.")

    def run(self, adata: ad.AnnData) -> ad.AnnData:
        """Master execution loop for Phase 3 feature engineering."""
        adata = self.stabilize_variance(adata)
        adata = self.extract_highly_variable_features(adata)
        adata = self.compute_pca_and_loadings(adata)
        self.save_integrated_workspace(adata)
        return adata
