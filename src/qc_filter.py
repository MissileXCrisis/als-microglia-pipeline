import os
import sqlite3
import logging
from pathlib import Path
import anndata as ad
import pandas as pd
import scanpy as sc

# Initialize production logging configuration
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

class QualityControlPipeline:
    """
    Handles Phase 2 of the pipeline: Ingestion via SQL, structural cell/gene filtering,
    mitochondrial degradation checks, and batch-aware doublet removal loops.
    """
    def __init__(self, config: dict):
        self.db_path = config['data_infrastructure']['db_path'] if 'data_infrastructure' in config else "data/als_microglia.db"
        self.min_genes = config['quality_control']['min_genes_per_cell']
        self.min_cells = config['quality_control']['min_cells_per_gene']
        self.max_mito = config['quality_control']['max_mitochondrial_pct']
        self.max_genes = config['quality_control'].get('max_genes_per_cell', 6000)

    def extract_cohort_from_db(self) -> pd.DataFrame:
        """Queries local SQLite instance to isolate target spinal cord cohorts."""
        logger.info(f"Querying relational database at {self.db_path} for spinal cord samples...")
        conn = sqlite3.connect(self.db_path)
        query = """
            SELECT gsm_id, donor_id, tissue, h5_filename 
            FROM cohort_metadata 
            WHERE tissue = 'SC' AND batch_protocol = 'Single-Run';
        """
        cohort_df = pd.read_sql_query(query, conn)
        conn.close()
        logger.info(f"Successfully identified {len(cohort_df)} cohort runs matching criteria.")
        return cohort_df

    def ingest_and_concatenate(self, cohort_df: pd.DataFrame) -> ad.AnnData:
        """Streams independent HDF5 matrices and merges them into a unified CSR matrix workspace."""
        adata_list = []
        for _, row in cohort_df.iterrows():
            h5_path = f"data/raw/{row['h5_filename']}"
            logger.info(f"Loading matrix {row['gsm_id']} ({row['donor_id']}) from {h5_path}...")
            
            if not os.path.exists(h5_path):
                raise FileNotFoundError(f"Missing mandatory input matrix: {h5_path}")
                
            sample_adata = sc.read_10x_h5(h5_path)
            sample_adata.var_names_make_unique()
            
            # Inject relational clinical metadata fields into cell observation space
            sample_adata.obs['gsm_id'] = row['gsm_id']
            sample_adata.obs['donor_id'] = row['donor_id']
            sample_adata.obs['tissue'] = row['tissue']
            
            adata_list.append(sample_adata)

        logger.info("Concatenating sample arrays into a single unified workspace matrix...")
        adata = ad.concat(adata_list, join="outer", label="batch", merge="unique")
        adata.X = adata.X.tocsr() # Enforce optimized Compressed Sparse Row layout
        logger.info(f"Unified Workspace Formed: {adata.shape[0]} cells x {adata.shape[1]} raw genes.")
        return adata

    def execute_quality_filters(self, adata: ad.AnnData) -> ad.AnnData:
        """Filters out zero-variance genes, empty droplets, and lysed/degraded cells."""
        logger.info(f"Pruning low-abundance transcripts appearing in fewer than {self.min_cells} cells...")
        sc.pp.filter_genes(adata, min_cells=self.min_cells)

        logger.info("Computing mitochondrial read metrics to flag cell lysis points...")
        adata.var['mt'] = adata.var_names.str.startswith('MT-')
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

        logger.info(f"Applying hard boundaries: Genes ({self.min_genes}-{self.max_genes}), Mito (<{self.max_mito*100}%)...")
        adata = adata[
            (adata.obs['pct_counts_mt'] < (self.max_mito * 100)) & 
            (adata.obs['n_genes_by_counts'] > self.min_genes) & 
            (adata.obs['n_genes_by_counts'] < self.max_genes), :
        ].copy()
        
        logger.info(f"Post-filtering background cleanup: {adata.shape[0]} cells remain.")
        return adata

    def remove_doublets(self, adata: ad.AnnData) -> ad.AnnData:
        """Executes a batch-aware doublet removal loop to permanently prune cellular collisions."""
        logger.info("Initiating batch-aware doublet simulation loops via Scrublet...")
        filtered_batches = []

        for batch_id in adata.obs['batch'].unique():
            batch_data = adata[adata.obs['batch'] == batch_id].copy()
            # Run scanpy wrapper over Scrublet engine to establish topological graph simulations
            sc.pp.scrublet(batch_data, verbose=False)
            filtered_batches.append(batch_data)

        # Recombine cleared cell blocks and drop predicted multiplets
        adata = ad.concat(filtered_batches, merge="unique")
        adata = adata[adata.obs['predicted_doublet'] == False, :].copy()
        
        logger.info(f"🎉 Phase 2 Complete! Cleaned Workspace Dimensions: {adata.shape[0]} cells x {adata.shape[1]} genes.")
        return adata

    def run(self) -> ad.AnnData:
        """Master execution loop for Phase 2 data cleaning."""
        cohort_df = self.extract_cohort_from_db()
        adata = self.ingest_and_concatenate(cohort_df)
        adata = self.execute_quality_filters(adata)
        adata = self.remove_doublets(adata)
        return adata
