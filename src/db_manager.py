import os
import gzip
import sqlite3
import logging
from pathlib import Path
import pandas as pd

# Set up clean logging infrastructure
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

class DatabaseManager:
    """
    Production-grade ETL pipeline manager that transforms non-standardized public 
    genomics metadata from NCBI GEO into a structured relational SQLite database.
    """
    def __init__(self, db_path: str = "data/als_microglia.db", 
                 matrix_file: str = "data/metadata/GSE204704_series_matrix.txt.gz"):
        self.db_path = db_path
        self.matrix_file = matrix_file
        
        # Ensure directories exist
        Path(self.db_path).parent.mkdir(parents=True, exist_ok=True)
        Path(self.matrix_file).parent.mkdir(parents=True, exist_ok=True)

    def setup_database(self, df: pd.DataFrame):
        """Populates the clean metadata dataframe into local SQLite instance."""
        try:
            logger.info(f"Connecting to relational database at: {self.db_path}")
            conn = sqlite3.connect(self.db_path)
            
            # Write or replace the table
            df.to_sql("cohort_metadata", conn, if_exists="replace", index=False)
            logger.info("✅ Table 'cohort_metadata' successfully written and indexed.")
            
            # Verify table contents
            verification_df = pd.read_sql_query("SELECT * FROM cohort_metadata", conn)
            logger.info(f"📊 Database Verification: Successfully retrieved {len(verification_df)} rows from disk.")
            
            conn.close()
            return True
        except Exception as e:
            logger.error(f"❌ Failed to populate database: {e}")
            raise e

    def parse_and_ingest_metadata(self) -> pd.DataFrame:
        """Parses the NCBI SOFT-formatted compressed series matrix file and feature-engineers metadata fields."""
        if not os.path.exists(self.matrix_file):
            logger.error(f"❌ Series matrix file not found at {self.matrix_file}. Run data download sequence first.")
            raise FileNotFoundError(f"Missing source file: {self.matrix_file}")

        logger.info(f"Opening and parsing compressed matrix file: {self.matrix_file}")
        
        # Initialize explicit target storage matrices
        gsm_ids = []
        titles = []
        tissues = []
        antibodies = []
        h5_filenames = []

        try:
            with gzip.open(self.matrix_file, 'rt') as f:
                for line in f:
                    # 1. Isolate core identifiers
                    if line.startswith("!Sample_title"):
                        titles = [x.strip('"') for x in line.split('\t')[1:]]
                    elif line.startswith("!Sample_geo_accession"):
                        gsm_ids = [x.strip('"') for x in line.split('\t')[1:]]
                        
                    # 2. Extract dynamic target metadata tracks
                    elif line.startswith("!Sample_characteristics_ch1"):
                        elements = [x.strip('"') for x in line.split('\t')[1:]]
                        
                        # Isolate physical tissue anatomical regions
                        if ":" in elements[0] and "tissue" in elements[0]:
                            tissues = [x.split(":")[-1].strip() for x in elements]
                            
                        # Isolate cellular hashing antibody/hashing batch controls
                        elif any("antibodies" in x for x in elements):
                            antibodies = [x.split(":")[-1].strip() if "antibodies" in x else "Single-Run" for x in elements]
                            
                    # 3. Pull supplementary HDF5 URL file endpoints
                    elif line.startswith("!Sample_supplementary_file_1"):
                        urls = [x.strip('"') for x in line.split('\t')[1:]]
                        h5_filenames = [os.path.basename(url) for url in urls]

            # In case the specific multiplexed row wasn't populated for all columns, default to Single-Run
            if not antibodies:
                antibodies = ["Single-Run"] * len(gsm_ids)

            # Assemble raw arrays into a unified DataFrame workspace
            df_clean = pd.DataFrame({
                'gsm_id': gsm_ids,
                'sample_title': titles,
                'tissue': tissues,
                'batch_protocol': antibodies,
                'h5_filename': h5_filenames
            })

            # 4. Feature Engineering Core Layouts:
            # Extract sample-level donor_id tokens directly from metadata labels (e.g. "ALS4 GM MA001" -> "ALS4")
            df_clean['donor_id'] = df_clean['sample_title'].apply(
                lambda x: x.split()[0] if len(x.split()) > 0 else "Unknown"
            )

            logger.info(f"Parsing complete. Unified shape: {df_clean.shape[0]} samples x {df_clean.shape[1]} columns.")
            
            # Synchronize records to disk
            self.setup_database(df_clean)
            return df_clean

        except Exception as e:
            logger.error(f"❌ Error occurred while parsing matrix metadata: {e}")
            raise e

    def fetch_cohort_by_tissue(self, tissue_name: str) -> pd.DataFrame:
        """Helper method to test programmatic SQL fetching of filtered cohorts."""
        conn = sqlite3.connect(self.db_path)
        query = "SELECT * FROM cohort_metadata WHERE tissue = ?"
        df_cohort = pd.read_sql_query(query, conn, params=(tissue_name,))
        conn.close()
        logger.info(f"Fetched {len(df_cohort)} records matching tissue cohort: '{tissue_name}'")
        return df_cohort


if __name__ == "__main__":
    # Self-test code path for local validation
    db_manager = DatabaseManager(
        db_path="data/als_microglia.db",
        matrix_file="data/metadata/GSE204704_series_matrix.txt.gz"
    )
    db_manager.parse_and_ingest_metadata()
