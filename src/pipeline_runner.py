import yaml
import logging
import scanpy as sc
from src.db_manager import DatabaseManager
from src.qc_filter import QualityControlPipeline
from src.feature_selection import FeatureSelectionPipeline
from src.cluster_engine import ClusterEnginePipeline
from src.clinical_analytics import ClinicalAnalyticsPipeline

# Global production logging configurations
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

def main():
    logger.info("==========================================================================")
    logger.info("🧬 RUNNING END-TO-END PRODUCTION ALS SINGLE-CELL PIPELINE RUNNER")
    logger.info("==========================================================================")

    # 1. Load centralized configuration architecture profiles
    try:
        with open("config/pipeline_config.yaml", "r") as f:
            config = yaml.safe_load(f)
        logger.info("✅ Configuration profiles successfully imported.")
    except Exception as e:
        logger.error(f"Failed to ingest centralized configuration file: {e}")
        return

    # 2. Phase 1: Database Construction & Data Streaming
    logger.info("--- STARTING PHASE 1: RELATIONAL INGESTION ---")
    db_engine = DatabaseManager(db_path="data/als_microglia.db")
    db_engine.parse_and_ingest_metadata()

    # 3. Phase 2: Quality Control Sanitization 
    logger.info("--- STARTING PHASE 2: DATA PROCESSING & SANITIZATION ---")
    qc_engine = QualityControlPipeline(config)
    adata_clean = qc_engine.run()

    # 4. Phase 3: Mathematical Dimension Compression 
    logger.info("--- STARTING PHASE 3: FEATURE SELECTION & REDUCTION ---")
    feature_engine = FeatureSelectionPipeline(config)
    adata_pca = feature_engine.run(adata_clean)

    # 5. Phase 4 & 5: Unsupervised Ecosystem Mapping & DE Testing
    logger.info("--- STARTING PHASE 4 & 5: CLUSTERING & DIFFERENTIAL EXPRESSION ---")
    cluster_engine = ClusterEnginePipeline(config)
    adata_clustered = cluster_engine.run(adata_pca)

    # 6. Phase 6: Clinical Abundance Analytics & Deployed Extractor
    logger.info("--- STARTING PHASE 6: CLINICAL ANALYTICS & PRODUCTION EXTRACTION ---")
    analytics_engine = ClinicalAnalyticsPipeline(config)
    _ = analytics_engine.run(adata_clustered)

    logger.info("==========================================================================")
    logger.info("🎉 SUCCESS: END-TO-END PIPELINE SHUTDOWN CLEANLY WITH ZERO ERRORS")
    logger.info("==========================================================================")

if __name__ == "__main__":
    main()