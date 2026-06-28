import logging
import scanpy as sc
import anndata as ad

logger = logging.getLogger(__name__)

class ClusterEnginePipeline:
    """
    Handles Phase 4 & Phase 5 of the production pipeline: Neighborhood graph engineering,
    Leiden community detection, UMAP coordinate projection, and non-parametric Wilcoxon DE testing.
    """
    def __init__(self, config: dict):
        self.neighbors_k = config['clustering']['neighbors_k']
        self.n_pcs = config['feature_selection']['n_pcs']
        self.resolution = config['clustering']['leiden_resolution']

    def construct_topology_and_cluster(self, adata: ad.AnnData) -> ad.AnnData:
        """Assembles the cell-to-cell neighborhood network graph and computes community partitions."""
        logger.info(f"Engineering topological k-NN graph network with k={self.neighbors_k} utilizing {self.n_pcs} PCs...")
        sc.pp.neighbors(adata, n_neighbors=self.neighbors_k, n_pcs=self.n_pcs)
        
        logger.info(f"Executing unsupervised Leiden community partitioning at resolution granularity: {self.resolution}...")
        sc.tl.leiden(adata, resolution=self.resolution)
        logger.info(f"✅ Leiden clustering complete. Total unique cell-states discovered: {adata.obs['leiden'].nunique()}")
        
        logger.info("Projecting network topology manifold down to continuous 2D UMAP spaces...")
        sc.tl.umap(adata)
        return adata

    def run_differential_expression(self, adata: ad.AnnData) -> ad.AnnData:
        """Executes global non-parametric Wilcoxon rank-sum testing with Benjamini-Hochberg FDR adjustments."""
        logger.info("Extracting clean log-normalized workspace matrix from .raw storage layer for DE alignment...")
        adata_de = adata.raw.to_adata()
        sc.pp.normalize_total(adata_de, target_sum=1e4)
        sc.pp.log1p(adata_de)

        logger.info("Running global non-parametric Wilcoxon rank-sum marker genes group testing...")
        sc.tl.rank_genes_groups(
            adata_de, 
            groupby='leiden', 
            method='wilcoxon', 
            corr_method='benjamini-hochberg',
            key_added='wilcoxon_global'
        )
        
        # Merge the computed differential expression metadata dictionary slots back to main workspace
        adata.uns['wilcoxon_global'] = adata_de.uns['wilcoxon_global']
        logger.info("✅ Global marker discovery complete. Statistical enrichment parameters stored in .uns workspace.")
        return adata

    def run(self, adata: ad.AnnData) -> ad.AnnData:
        """Master execution loop for the clustering engine layer."""
        adata = self.construct_topology_and_cluster(adata)
        adata = self.run_differential_expression(adata)
        return adata
