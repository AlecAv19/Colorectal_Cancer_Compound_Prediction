# -*- coding: utf-8 -*-
"""
Created on Tue Aug 19 16:23:44 2025

@author: ASUS TUF F15
"""
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
import warnings
warnings.filterwarnings('ignore')

def integrate_virtual_screening_data(excel_file_path):
    """
    Integrate virtual screening data from NCI, CHEMBL, and GDSC sources
    with prioritized labels for virtual screening applications.
    
    Args:
        excel_file_path (str): Path to the VC_Data.xlsx file
    
    Returns:
        pd.DataFrame: Integrated dataset with optimized VS features
    """
    
    # Read data from all three sources
    print("Reading data from three sources...")
    nci_data = pd.read_excel(excel_file_path, sheet_name='NCI')
    chembl_data = pd.read_excel(excel_file_path, sheet_name='CHEMBL')
    gdsc_data = pd.read_excel(excel_file_path, sheet_name='GDSC')
    
    print(f"NCI: {len(nci_data)} compounds")
    print(f"CHEMBL: {len(chembl_data)} compounds") 
    print(f"GDSC: {len(gdsc_data)} compounds")
    
    # Clean and standardize SMILES
    def clean_smiles(smiles):
        if pd.isna(smiles) or smiles == 'null' or smiles == '':
            return None
        return str(smiles).strip()
    
    # Prepare NCI data with VS-relevant features
    nci_processed = nci_data.copy()
    nci_processed['canonical_smiles'] = nci_processed['canonical_smiles'].apply(clean_smiles)
    nci_processed['source'] = 'NCI'
    nci_processed['data_priority'] = 1  # Highest priority for VS
    
    # Select most relevant NCI features for VS
    nci_features = [
        'canonical_smiles', 'compound_name', 'source', 'data_priority',
        # Core VS features
        'consensus_activity', 'vs_score', 'activity_ratio',
        # Potency features  
        'mean_gi50', 'std_gi50', 'min_gi50', 'max_gi50',
        # Reliability features
        'confidence', 'gi50_consistency', 'has_conflict',
        # Selectivity features
        'selectivity', 'num_cell_lines_tested', 'active_cell_lines', 'inactive_cell_lines',
        # ML-ready features
        'tier', 'tier_description', 'ml_weight', 'training_priority',
        # Additional context
        'cell_lines_tested'
    ]
    
    nci_vs = nci_processed[nci_features].copy()
    
    # Prepare CHEMBL data
    chembl_processed = chembl_data.copy()
    chembl_processed['canonical_smiles'] = chembl_processed['canonical_smiles'].apply(clean_smiles)
    chembl_processed['source'] = 'CHEMBL'
    chembl_processed['data_priority'] = 2
    
    # Map CHEMBL features to common VS schema
    chembl_vs = pd.DataFrame({
        'canonical_smiles': chembl_processed['canonical_smiles'],
        'compound_name': chembl_processed['molregno_chembl_id'],  # Use ChEMBL ID as name
        'source': chembl_processed['source'],
        'data_priority': chembl_processed['data_priority'],
        'consensus_activity': chembl_processed['class'],  # Map class to consensus_activity
        'vs_score': None,  # Not available in CHEMBL
        'activity_ratio': None,
        'mean_gi50': chembl_processed['standard_value'],  # Map to potency measure
        'std_gi50': None,
        'min_gi50': None,
        'max_gi50': None,
        'confidence': None,
        'gi50_consistency': None,
        'has_conflict': None,
        'selectivity': None,
        'num_cell_lines_tested': None,
        'active_cell_lines': None,
        'inactive_cell_lines': None,
        'tier': None,
        'tier_description': None,
        'ml_weight': None,
        'training_priority': None,
        'cell_lines_tested': chembl_processed['All_CellLines'],
        # CHEMBL-specific features
        'chembl_id': chembl_processed['molregno_chembl_id'],
        'num_assays': chembl_processed['num_assays'],
        'cell_chembl_id': chembl_processed['cell_chembl_id']
    })
    
    # Prepare GDSC data
    gdsc_processed = gdsc_data.copy()
    gdsc_processed['canonical_smiles'] = gdsc_processed['canonical_smiles'].apply(clean_smiles)
    gdsc_processed['source'] = 'GDSC'
    gdsc_processed['data_priority'] = 3
    
    # Map GDSC features to common VS schema
    gdsc_vs = pd.DataFrame({
        'canonical_smiles': gdsc_processed['canonical_smiles'],
        'compound_name': gdsc_processed['drug_name'],
        'source': gdsc_processed['source'],
        'data_priority': gdsc_processed['data_priority'],
        'consensus_activity': gdsc_processed['activity'],
        'vs_score': None,
        'activity_ratio': gdsc_processed['consistency_ratio'],
        'mean_gi50': gdsc_processed['mean_z_score'],  # Z-score as potency measure
        'std_gi50': gdsc_processed['std_z_score'],
        'min_gi50': gdsc_processed['min_z_score'],
        'max_gi50': gdsc_processed['max_z_score'],
        'confidence': None,
        'gi50_consistency': gdsc_processed['consistency_ratio'],
        'has_conflict': None,
        'selectivity': None,
        'num_cell_lines_tested': gdsc_processed['num_cell_lines'],
        'active_cell_lines': None,
        'inactive_cell_lines': None,
        'tier': None,
        'tier_description': None,
        'ml_weight': None,
        'training_priority': None,
        'cell_lines_tested': gdsc_processed['cell_lines_tested'],
        # GDSC-specific features
        'targets': gdsc_processed['targets'],
        'z_score_range': gdsc_processed['z_score_range'],
        'num_measurements': gdsc_processed['num_measurements']
    })
    
    # Combine all datasets
    print("\nCombining datasets...")
    combined_data = pd.concat([nci_vs, chembl_vs, gdsc_vs], ignore_index=True, sort=False)
    
    # Remove rows without SMILES
    combined_data = combined_data.dropna(subset=['canonical_smiles'])
    
    print(f"Combined dataset: {len(combined_data)} entries")
    
    return combined_data


def create_integrated_vs_dataset(combined_data):
    """
    Create a consolidated virtual screening dataset by merging overlapping compounds
    and prioritizing the most reliable data sources.
    """
    
    print("\nCreating integrated VS dataset...")
    
    # Group by SMILES to handle overlapping compounds
    grouped = combined_data.groupby('canonical_smiles')
    
    integrated_records = []
    
    for smiles, group in grouped:
        if len(group) == 1:
            # Single source - use as is but add metadata
            record = group.iloc[0].to_dict()
            record['num_data_sources'] = 1
            record['data_sources'] = record['source']
            record['composite_confidence'] = 0.33  # Base confidence for single source
            integrated_records.append(record)
        else:
            # Multiple sources - create integrated record
            print(f"Integrating {len(group)} sources for compound: {group.iloc[0]['compound_name']}")
            
            # Sort by data priority (NCI=1, CHEMBL=2, GDSC=3)
            group_sorted = group.sort_values('data_priority')
            
            # Start with highest priority source
            integrated_record = group_sorted.iloc[0].to_dict()
            
            # Enhance with data from other sources
            for _, row in group_sorted.iloc[1:].iterrows():
                # Fill missing values with data from other sources
                for col in integrated_record:
                    if pd.isna(integrated_record[col]) and not pd.isna(row[col]):
                        integrated_record[col] = row[col]
                
                # Add source-specific information
                if row['source'] == 'CHEMBL':
                    integrated_record['chembl_validation'] = True
                    integrated_record['chembl_num_assays'] = row['num_assays']
                elif row['source'] == 'GDSC':
                    integrated_record['gdsc_validation'] = True
                    integrated_record['known_targets'] = row.get('targets', None)
            
            # Create composite confidence score
            sources = group['source'].tolist()
            confidence_score = len(sources) / 3.0  # Normalize by max possible sources
            
            if 'NCI' in sources:
                nci_conf = group[group['source'] == 'NCI']['confidence'].iloc[0]
                if nci_conf == 'high':
                    confidence_score += 0.2
            
            integrated_record['composite_confidence'] = confidence_score
            integrated_record['num_data_sources'] = len(sources)
            integrated_record['data_sources'] = ','.join(sources)
            
            integrated_records.append(integrated_record)
    
    integrated_df = pd.DataFrame(integrated_records)
    
    print(f"Integrated dataset: {len(integrated_df)} unique compounds")
    print(f"Multi-source compounds: {len([r for r in integrated_records if r['num_data_sources'] > 1])}")
    
    return integrated_df


def add_molecular_descriptors(df, smiles_col='canonical_smiles'):
    """Add key molecular descriptors for virtual screening."""
    
    print("\nCalculating molecular descriptors...")
    
    descriptors = []
    failed_smiles = []
    
    for idx, smiles in enumerate(df[smiles_col]):
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                desc_dict = {
                    'molecular_weight': Descriptors.MolWt(mol),
                    'logp': Descriptors.MolLogP(mol),
                    'hbd': Descriptors.NumHDonors(mol),
                    'hba': Descriptors.NumHAcceptors(mol),
                    'tpsa': Descriptors.TPSA(mol),
                    'rotatable_bonds': Descriptors.NumRotatableBonds(mol),
                    'aromatic_rings': Descriptors.NumAromaticRings(mol),
                    'lipinski_violations': sum([
                        Descriptors.MolWt(mol) > 500,
                        Descriptors.MolLogP(mol) > 5,
                        Descriptors.NumHDonors(mol) > 5,
                        Descriptors.NumHAcceptors(mol) > 10
                    ])
                }
                descriptors.append(desc_dict)
            else:
                failed_smiles.append(smiles)
                descriptors.append({k: None for k in ['molecular_weight', 'logp', 'hbd', 'hba', 'tpsa', 'rotatable_bonds', 'aromatic_rings', 'lipinski_violations']})
        except Exception as e:
            failed_smiles.append(smiles)
            descriptors.append({k: None for k in ['molecular_weight', 'logp', 'hbd', 'hba', 'tpsa', 'rotatable_bonds', 'aromatic_rings', 'lipinski_violations']})
    
    if failed_smiles:
        print(f"Warning: Failed to process {len(failed_smiles)} SMILES")
    
    # Add descriptors to dataframe
    desc_df = pd.DataFrame(descriptors)
    result_df = pd.concat([df.reset_index(drop=True), desc_df], axis=1)
    
    return result_df


def create_vs_labels(df):
    """Create standardized labels for virtual screening."""
    
    print("\nCreating VS labels...")
    
    # Standardize activity labels
    def standardize_activity(row):
        activity = row.get('consensus_activity')
        if not pd.isna(activity):
            if activity in [1, '1', 'active', 'Active']:
                return 1
            elif activity in [0, '0', 'inactive', 'Inactive']:
                return 0
        return None
    
    df['vs_activity_label'] = df.apply(standardize_activity, axis=1)
    
    # Create confidence categories
    def categorize_confidence(row):
        composite_conf = row.get('composite_confidence', 0)
        if pd.isna(composite_conf):
            composite_conf = 0
        
        if composite_conf >= 0.8:
            return 'high'
        elif composite_conf >= 0.5:
            return 'medium'
        else:
            return 'low'
    
    df['vs_confidence_category'] = df.apply(categorize_confidence, axis=1)
    
    # Create training set recommendations
    def recommend_for_training(row):
        score = 0
        
        # Clear activity label
        if not pd.isna(row.get('vs_activity_label')):
            score += 3
        
        # Confidence level
        conf_cat = row.get('vs_confidence_category', 'low')
        if conf_cat == 'high':
            score += 2
        elif conf_cat == 'medium':
            score += 1
            
        # Multiple data sources
        num_sources = row.get('num_data_sources', 1)
        if not pd.isna(num_sources) and num_sources > 1:
            score += 2
            
        # Drug-like properties
        lipinski_viol = row.get('lipinski_violations', 5)
        if not pd.isna(lipinski_viol) and lipinski_viol <= 1:
            score += 1
            
        return score >= 5
    
    df['recommend_for_training'] = df.apply(recommend_for_training, axis=1)
    
    # Create VS priority tiers
    def assign_vs_tier(row):
        training_rec = row.get('recommend_for_training', False)
        conf_cat = row.get('vs_confidence_category', 'low')
        activity_label = row.get('vs_activity_label')
        
        if training_rec and conf_cat == 'high':
            return 'Tier_1_Premium'
        elif activity_label is not None and conf_cat in ['high', 'medium']:
            return 'Tier_2_Reliable'
        elif activity_label is not None:
            return 'Tier_3_Usable'
        else:
            return 'Tier_4_Exploratory'
    
    df['vs_tier'] = df.apply(assign_vs_tier, axis=1)
    
    return df


def generate_summary_report(df):
    """Generate a summary report of the integrated dataset."""
    
    print("\n" + "="*60)
    print("VIRTUAL SCREENING DATASET INTEGRATION SUMMARY")
    print("="*60)
    
    print(f"\n📊 DATASET OVERVIEW:")
    print(f"   Total compounds: {len(df)}")
    print(f"   Unique SMILES: {df['canonical_smiles'].nunique()}")
    
    print(f"\n⚡ ACTIVITY LABELS:")
    if 'vs_activity_label' in df.columns:
        activity_counts = df['vs_activity_label'].value_counts()
        print(f"   Active compounds: {activity_counts.get(1, 0)}")
        print(f"   Inactive compounds: {activity_counts.get(0, 0)}")
        print(f"   Unlabeled compounds: {df['vs_activity_label'].isna().sum()}")
    
    print(f"\n🔬 CONFIDENCE DISTRIBUTION:")
    if 'vs_confidence_category' in df.columns:
        conf_counts = df['vs_confidence_category'].value_counts()
        for conf, count in conf_counts.items():
            print(f"   {conf.capitalize()}: {count}")
    
    print(f"\n🎯 VS TIER DISTRIBUTION:")
    if 'vs_tier' in df.columns:
        tier_counts = df['vs_tier'].value_counts()
        for tier, count in tier_counts.items():
            print(f"   {tier}: {count}")
    
    print(f"\n🧪 TRAINING SET RECOMMENDATIONS:")
    if 'recommend_for_training' in df.columns:
        training_ready = df['recommend_for_training'].sum()
        print(f"   Recommended for training: {training_ready}")
        print(f"   Training set percentage: {training_ready/len(df)*100:.1f}%")
    
    print(f"\n💊 MOLECULAR PROPERTIES:")
    if 'lipinski_violations' in df.columns:
        lipinski_compliant = (df['lipinski_violations'] <= 1).sum()
        print(f"   Lipinski compliant: {lipinski_compliant}/{len(df)} ({lipinski_compliant/len(df)*100:.1f}%)")
    
    print(f"\n🔗 DATA SOURCE INTEGRATION:")
    if 'num_data_sources' in df.columns:
        multi_source = (df['num_data_sources'] > 1).sum()
        print(f"   Multi-source compounds: {multi_source}")
        print(f"   Single-source compounds: {len(df) - multi_source}")
    
    # Top compounds for VS
    print(f"\n🌟 TOP COMPOUNDS FOR VIRTUAL SCREENING:")
    if 'vs_tier' in df.columns:
        top_compounds = df[df['vs_tier'] == 'Tier_1_Premium'].head(5)
        if len(top_compounds) == 0:
            top_compounds = df[df['vs_tier'] == 'Tier_2_Reliable'].head(5)
        
        for idx, row in top_compounds.iterrows():
            data_sources = row.get('data_sources', 'Unknown')
            activity = row.get('vs_activity_label', 'Unknown')
            compound_name = row.get('compound_name', 'Unknown')
            print(f"   • {compound_name} - Activity: {activity}, Sources: {data_sources}")


def save_vs_strategy_summary(df, output_file):
    """Create VS strategy comparison sheet."""
    
    print("\nCreating VS strategy summary...")
    
    strategies = {
        'All_Compounds': df,
        'High_Confidence_Only': df[df.get('vs_confidence_category', '') == 'high'],
        'Multi_Source_Only': df[df.get('num_data_sources', 1) > 1],
        'Training_Ready': df[df.get('recommend_for_training', False) == True],
        'NCI_Priority': df[df.get('source', '') == 'NCI']
    }
    
    strategy_summary = []
    
    for strategy_name, subset in strategies.items():
        if len(subset) > 0:
            summary = {
                'Strategy': strategy_name,
                'Total_Compounds': len(subset),
                'Active_Compounds': subset[subset.get('vs_activity_label', None) == 1].shape[0] if 'vs_activity_label' in subset.columns else 0,
                'Inactive_Compounds': subset[subset.get('vs_activity_label', None) == 0].shape[0] if 'vs_activity_label' in subset.columns else 0,
                'High_Confidence': subset[subset.get('vs_confidence_category', '') == 'high'].shape[0],
                'Multi_Source': subset[subset.get('num_data_sources', 1) > 1].shape[0],
                'Lipinski_Compliant': subset[subset.get('lipinski_violations', 5) <= 1].shape[0] if 'lipinski_violations' in subset.columns else 0
            }
        else:
            summary = {
                'Strategy': strategy_name,
                'Total_Compounds': 0,
                'Active_Compounds': 0,
                'Inactive_Compounds': 0,
                'High_Confidence': 0,
                'Multi_Source': 0,
                'Lipinski_Compliant': 0
            }
        
        strategy_summary.append(summary)
    
    strategy_df = pd.DataFrame(strategy_summary)
    
    # Save with multiple sheets
    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        df.to_excel(writer, sheet_name='Integrated_Dataset', index=False)
        strategy_df.to_excel(writer, sheet_name='VS_Strategy_Summary', index=False)
        
        # Individual strategy sheets
        for strategy_name, subset in strategies.items():
            if len(subset) > 0:
                subset.to_excel(writer, sheet_name=f'VS_{strategy_name}', index=False)
    
    print(f"✅ VS strategy summary saved with {len(strategies)} different approaches")
    return strategy_df


# Main execution function
def main():
    """Main function to execute the virtual screening data integration."""
    
    # File path - update this to your actual file location
    excel_file = 'VC_Data.xlsx'
    
    try:
        # Step 1: Load and process raw data
        print("🚀 Starting Virtual Screening Data Integration...")
        combined_data = integrate_virtual_screening_data(excel_file)
        
        # Step 2: Create integrated dataset
        integrated_data = create_integrated_vs_dataset(combined_data)
        
        # Step 3: Add molecular descriptors
        try:
            enriched_data = add_molecular_descriptors(integrated_data)
            print("✅ Molecular descriptors calculated successfully")
        except Exception as e:
            print(f"⚠️ Warning: Could not calculate molecular descriptors: {e}")
            print("   Continuing without RDKit descriptors...")
            enriched_data = integrated_data.copy()
            # Add placeholder columns
            enriched_data['molecular_weight'] = None
            enriched_data['logp'] = None
            enriched_data['lipinski_violations'] = None
        
        # Step 4: Create VS labels and tiers
        final_data = create_vs_labels(enriched_data)
        
        # Step 5: Generate summary report
        generate_summary_report(final_data)
        
        # Step 6: Save the integrated dataset with strategy comparison
        output_file = 'Integrated_VS_Dataset_Complete.xlsx'
        strategy_summary = save_vs_strategy_summary(final_data, output_file)
        
        # Also save a simple CSV version
        csv_output = 'Integrated_VS_Dataset.csv'
        final_data.to_csv(csv_output, index=False)
        
        print(f"\n🎉 SUCCESS! Integration completed!")
        print(f"📁 Main results: {output_file}")
        print(f"📁 CSV version: {csv_output}")
        print(f"💡 Next steps:")
        print(f"  1. Review 'VS_Strategy_Summary' sheet to choose approach")
        print(f"  2. Use individual VS sheets for different screening strategies")
        print(f"  3. Apply compound tiers and confidence scores for prioritization")
        
        return final_data, strategy_summary
        
    except Exception as e:
        print(f"❌ Error during integration: {str(e)}")
        import traceback
        traceback.print_exc()
        return None, None


if __name__ == "__main__":
    # Run the integration
    integrated_dataset, strategy_summary = main()
    
    # Display sample results if successful
    if integrated_dataset is not None:
        print("\n" + "="*60)
        print("SAMPLE OF INTEGRATED DATASET (Key Columns)")
        print("="*60)
        
        # Show key columns that should always exist
        display_columns = ['compound_name', 'source', 'vs_activity_label', 'vs_confidence_category', 'vs_tier']
        available_columns = [col for col in display_columns if col in integrated_dataset.columns]
        
        if available_columns:
            print(integrated_dataset[available_columns].head(10).to_string(index=False))
        else:
            print("Sample data:")
            print(integrated_dataset.head(5).to_string(index=False))