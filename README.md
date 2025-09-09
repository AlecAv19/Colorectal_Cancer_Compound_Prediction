# Colorectal_Cancer_Compound_Prediction

## Virtual Screening Pipeline for Colorectal Cancer Drug Discovery
A comprehensive machine learning pipeline for virtual screening of natural products and small molecules against colorectal cancer cell lines.

### Overview
This pipeline implements a complete virtual screening workflow including SMILES processing, data balancing, machine learning model training, virtual screening, and toxicity assessment. The system is specifically designed for colorectal cancer drug discovery using validated cell line models.

### Features
#### - Multi-step ML Pipeline: From raw SMILES to validated predictions
#### - Data Balancing: Automated handling of imbalanced datasets
#### - Multiple Algorithms: Support for RF, DTREE, KNN, BNB, GNB
#### - Ensemble Modeling: Advanced genetic algorithm-based model optimization
#### - Virtual Screening: Large-scale compound screening capabilities
#### - Toxicity Assessment: Safety evaluation of hit compounds
#### - Natural Product Focus: COCONUT, LOTUS, and FOODB database support

### Step 1: SMILES Processing and Data Balancing

```
python Step1SmilesBalancingPartition.py -i input_file.txt -p Step1ParamFile.txt
```
#### - Sanitizes and validates SMILES structures
#### - Computes molecular descriptors
#### - Balances datasets using clustering
#### - Creates train/test/external splits

### Step 2: Data Preprocessing

```
python Step2DataProcessing1.py -i training_file.txt -p Step2ParamFile.txt
```
#### - Data normalization (StandardScaler)
#### - Feature selection and reduction
#### - Saves preprocessing parameters for consistency



