# Mapping Microglial Heterogeneity in Amyotrophic Lateral Sclerosis (ALS) Using Single-Cell Transcriptomics

🚀 **Live Production Dashboard Portal:** https://als-microglia-pipeline-aappnfbukwnn9f2fpayeifw.streamlit.app/

---

## 🔬 Project Executive Summary
This production-grade single-cell pipeline maps microglial cell-state transitions and immune ecosystems in human post-mortem spinal cord tissue across the public genomics dataset **GSE204704**. By integrating a relational database framework (SQLite) with batch-aware quality control, high-dimensional manifold learning, and automated cell-type annotation, this architecture isolates structural neuronal contamination and tracks peripheral T-cell infiltration into the central nervous system.

## 🏗️ Repository Architecture
Adhering to strict data science infrastructure and reproducibility standards, this production stack is modularized into discrete functional components. 

```text
als-microglia-pipeline/
├── .gitignore                  # Prevents committing raw data and heavy DB files to Git
├── README.md                   # Complete pipeline documentation and abstract
├── requirements.txt            # Explicit environment dependency freeze file
├── config/
│   └── pipeline_config.yaml    # Centralized control file for hard QC and PCA thresholds
├── data/                       # Kept local / Git-ignored
│   ├── metadata/               # Raw series matrix files
│   ├── raw/                    # Raw unzipped sample .h5 matrices 
│   ├── processed/              # Stores checkpoint data layer (dashboard_data.csv)
│   └── als_microglia.db        # Phase 1 SQLite relational database instance
├── notebooks/                  # Preserved strictly for exploratory R&D and visualization
│   ├── Phase1_Data.ipynb
│   ├── Phase2&3_Processing.ipynb
│   └── Phase_4_5_6.ipynb
└── src/                        # Core production source code directory
    ├── __init__.py             # Recognizes src as a python package
    ├── db_manager.py           # Phase 1: SQL Database pipeline construction
    ├── qc_filter.py            # Phase 2: Transcriptomic sanitization and Doublet loops
    ├── feature_selection.py    # Phase 3: Highly Variable Gene (HVG) selection & PCA loading math
    ├── cluster_engine.py       # Phase 4 & 5: Leiden clustering and Wilcoxon DE marker discovery
    ├── clinical_analytics.py   # Phase 6: Annotation mapping and compositional modeling
    └── app.py                  # Streamlit Web Application entry point
	
🧬 Core Biological Conclusions
1. High-Dimensional Manifold & Cluster Topologies (Phase 3 & 4)
Mathematical Gradients of Active Phagocytosis (PC1): Principal Component 1 (PC1) isolated strong positive weights for structural neuronal/axonal markers (PCDH9, NTM, NOVA1, NCAM1), mathematically capturing ambient tissue contamination alongside active microglial phagocytosis of degenerating motor neurons.

Peripheral Immune Infiltration (PC2): PC2 isolated a strict immunological axis dominated by T-cell receptor components and lymphocytic cytokines (CCL5, CD2, IL32, CD3E, CD3D), drawing a hard spatial boundary between central resident microglia and infiltrating peripheral T-lymphocytes crossing the disrupted blood-brain barrier.

2. Differential Expression Unmasks the DAM Phenotype (Phase 5)
Global non-parametric Wilcoxon rank-sum testing paired with Benjamini-Hochberg FDR corrections identified Cluster 8 as the highly active Disease-Associated Microglia (DAM) state, showing robust, exclusive transcript enrichment of hallmarks SPP1, RGCC, and LPL.

3. Patient Stratification & Clinical Abundance Modeling (Phase 6)
Cross-tabulation and clinical stratification across the donor cohorts successfully captured extensive patient-to-patient heterogeneity from autopsy specimens:

Donor ALS3 presented a highly localized neurodegenerative climax, leading the cohort with a 7.22% DAM phenotype burden.

Donor ALS4 exhibited a severe structural breakdown of tissue protection, marked by a massive 34.92% peripheral T-cell influx.

Sex-Linked QC Verification: Tracking sex-specific biomarker expression safely profiled donor ALS4 as biologically female due to a distinct 21.89% XIST+ microglial footprint, verifying internal dataset controls.

💻 Quickstart & Reproducibility Guide (Debian/Linux)
Follow these terminal instructions to instantiate a clean-room environment replication, stream the source repositories, process the raw sequencing arrays, and deploy the application layer locally.
sudo apt update
sudo apt install git -y 

git clone https://github.com/MissileXCrisis/als-microglia-pipeline.git 

# Create a fresh Python virtual environment
sudo apt update
sudo apt install python3 python3-venv python3-pip -y
python3 -m venv venv_test

# Activate the isolated testing environment
source venv_test/bin/activate 

# Upgrade the core package manager
pip install --upgrade pip 
# Install the exact software stack from your requirements file 
pip install -r requirements.txt 

# Create the local data architecture layouts 
mkdir -p data/raw data/metadata 

#Direct download stream if testing completely from scratch:
wget -O data/metadata/GSE204704_series_matrix.txt.gz "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE204nnn/GSE204704/matrix/GSE204704_series_matrix.txt.gz"
wget -O GSE204704_RAW.tar "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE204704&format=file"
tar -xvf GSE204704_RAW.tar -C ./data/raw/

#Execute the master Orchestrator
python -m src.pipeline_runner

#Execture the stream lit  dashboard from the local data
streamlit run src/app.py