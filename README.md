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

### Step 3: Model Training

```
python Step3ADataModelling.py -p Step3AParamFile.txt
```
#### - Genetic algorithm-based feature selection
#### - Multiple ML algorithm training
#### - Ensemble model optimization
#### - Application domain definition
#### - Comprehensive model evaluation

### Step 4: Model Evaluation

```
python Step4ModelEvaluation.py
```
#### - Performance comparison across models
#### - Model ranking and selection
#### - Deployment readiness assessment

### Step 5: DECOY Generation

```
python Step5Decoy_Generation.py
```
#### - Generates property-matched decoy compounds
#### - Creates realistic virtual screening scenarios
#### - Supports high active:inactive ratios

### Step 6: Virtual Screening

```
python Step6Virtual_Screening.py
```
#### - Automated model and data detection
#### - Large-scale compound screening
#### - Confidence scoring and ranking

### Step 7: Multi-Model Comparison

```
python Step7Screening_Comparison.py
```
#### - Compares results across different models
#### - Consensus hit identification
#### - Performance analytics

### Step 8: Database Screening

```
python Step8Database_Screening.py
```
#### - Natural product database screening
#### - COCONUT, LOTUS, FOODB support
#### - Cross-database consensus analysis

### Step 9: Toxicity Assessment

```
python Step9Toxicity_Screening.py
```
#### - Safety evaluation of top hits
#### - Multiple toxicity endpoints
#### - Risk prioritization

## Installation

### Requirements
```
python >= 3.7
pandas >= 1.3.0
numpy >= 1.20.0
scikit-learn >= 1.0.0
rdkit-pypi >= 2022.3.1
openpyxl >= 3.0.0
```

### Dependencies
```
pip install pandas numpy scikit-learn rdkit-pypi openpyxl
```

## Quick Start

### 1. Prepare Your Data
#### Create input file with format:
```
ID	SMILES	Class
Compound1	CCO	1
Compound2	c1ccccc1	0
```

### 2. Configure Parameters
#### Edit parameter files for each step:
#### - Step1ParamFile.txt - SMILES processing options
#### - Step2ParamFile.txt - Preprocessing parameters
#### - Step3AParamFile.txt - Model training settings

### 3. Run Complete Pipeline
```
# Step 1: Data processing
python Step1SmilesBalancingPartition.py -i your_data.txt -p Step1ParamFile.txt

# Step 2: Preprocessing
python Step2DataProcessing1.py -i your_data_Curated_Balanced_Training.txt -p Step2ParamFile.txt

# Step 3A: Model training
python Step3ADataModelling.py -p Step3AParamFile.txt

# Step 4: Model evaluation
python Step4ModelEvaluation.py

# Continue with virtual screening steps...
```

## Parameter Configuration
### Step 1 Parameters (Step1ParamFile.txt)

```
Fix_Salts = True
Include_Isomers = False
Ignore_Charge = True
Descriptors_Type = 1
Partition_Ratio = 0.3
Metric_Cluster = euclidean
N_Cluster_Min = 5
Split_Ratio = 70-15-15
```

### Step 2 Parameters (Step2ParamFile.txt)

```
Data_Normalization = True
Simple_Variable_Reduction = True
Var_Curoff = 0.1
```

### Step 3A Parameters (Step3AParamFile.txt)

```
Training_File = your_training_file.txt
Test_File = your_test_file.txt
External_File = your_external_file.txt
Normalization_File = your_normalization_file.pkl
Fitness_Functions = RF, DTREE, KNN
Fitness_Metric = BCR
population_size = 50
num_generations = 100
```

## Output Files
### Model Files

#### *_Model.pkl - Trained models for deployment
#### *_Scaler.pkl - Preprocessing parameters
#### ApplicationDomain_Conditions.pkl - Application domain definition

### Results Files

#### VS_Results_*.xlsx - Virtual screening results
#### Model_Performance_Summary.csv - Model comparison
#### Top50_Focused_Toxicity_Results_*.xlsx - Safety assessment

### Log Files

#### *_logFile.log - Detailed processing logs for each step

## Cell Line Models
### The pipeline supports multiple colorectal cancer cell line models:

#### HCT-15: Best overall performer (Quality: 54.6/100)
#### LoVo: Strong consensus model (Quality: 51.6/100)
#### DLD-1: Conservative, low false positives (Quality: 46.1/100)

## Database Support
### Natural Product Databases

#### - COCONUT: >400k natural products
#### - LOTUS: Traditional medicine compounds
#### - FOODB: Food-derived compounds

### DECOY Databases

#### - Property-matched inactive compounds
#### - Configurable active:inactive ratios

## Performance Metrics
### The pipeline uses multiple metrics for model evaluation:

#### - BCR: Balanced Classification Rate
#### - Accuracy: Overall prediction accuracy
#### - F1-Score: Harmonic mean of precision/recall
#### - Sensitivity/Specificity: Individual class performance

## Troubleshooting
### Common Issues

###   1. SMILES Processing Errors

####   - Check SMILES validity with RDKit
####   - Verify input file format


###   2. Memory Issues

####   - Reduce dataset size for testing
####   - Use data chunking for large files


###   3. Model Loading Errors

####   - Ensure all preprocessing files are available
####   - Check file paths in parameter files

## Citation
### If you use this pipeline in your research, please cite:
```
@software{virtual_screening_pipeline,
  title={Virtual Screening Pipeline for Colorectal Cancer Drug Discovery},
  author={[Alec Avila]},
  year={2024},
  url={https://github.com/[AlecAv19]/virtual-screening-pipeline}
}
```

## License
### This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments
### RDKit for molecular descriptor calculation
### scikit-learn for machine learning algorithms
### Natural product databases: COCONUT, LOTUS, FOODB

