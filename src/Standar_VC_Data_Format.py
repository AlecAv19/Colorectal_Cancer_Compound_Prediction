# -*- coding: utf-8 -*-
"""
Created on Tue Aug 19 16:41:30 2025

@author: ASUS TUF F15
"""

"""
BRIDGE SCRIPT: VS Data to Decoy Generator Format Converter

This script converts your integrated VS dataset from the 3 sources (NCI, CHEMBL, GDSC)
into the format expected by the decoy generation script.

Input: Integrated_VS_Dataset_Complete.xlsx (from VS integration script)
Output: compounds_across_celllines.xlsx (for decoy generation script)
"""

import pandas as pd
import numpy as np
import os
from datetime import datetime

def convert_vs_data_for_decoy_generation(vs_input_file='Integrated_VS_Dataset_Complete.xlsx'):
    """
    Convert integrated VS dataset to format expected by decoy generator
    """
    
    print("🔄 CONVERTING VS DATA FOR DECOY GENERATION")
    print("=" * 50)
    
    try:
        # Check if input file exists
        if not os.path.exists(vs_input_file):
            print(f"❌ Input file not found: {vs_input_file}")
            print("Please run the VS integration script first to create this file.")
            return False
        
        print(f"📁 Loading VS data from: {vs_input_file}")
        
        # Load the integrated dataset
        vs_data = pd.read_excel(vs_input_file, sheet_name='Integrated_Dataset')
        
        print(f"✅ Loaded {len(vs_data)} compounds from VS dataset")
        
        # Check required columns
        required_columns = ['canonical_smiles', 'vs_activity_label']
        missing_columns = [col for col in required_columns if col not in vs_data.columns]
        
        if missing_columns:
            print(f"❌ Missing required columns: {missing_columns}")
            print("Available columns:", list(vs_data.columns))
            return False
        
        # Clean and prepare the data
        print("\n🧹 Cleaning and preparing data...")
        
        # Remove rows without SMILES or activity labels
        clean_data = vs_data.dropna(subset=['canonical_smiles', 'vs_activity_label']).copy()
        
        print(f"   Compounds with valid SMILES and activity: {len(clean_data)}")
        
        # Create the format expected by decoy generator
        print("\n🔄 Converting to decoy generator format...")
        
        # Map compound names - use the best available identifier
        def get_compound_id(row):
            if pd.notna(row.get('compound_name')) and row.get('compound_name') != '':
                return str(row['compound_name'])
            elif pd.notna(row.get('chembl_id')) and row.get('chembl_id') != '':
                return str(row['chembl_id'])
            elif pd.notna(row.get('nsc_id')) and row.get('nsc_id') != '':
                return str(row['nsc_id'])
            else:
                return f"COMPOUND_{row.name}"
        
        # Create standardized dataset
        standardized_data = pd.DataFrame({
            'molregno_chembl_id': clean_data.apply(get_compound_id, axis=1),
            'canonical_smiles': clean_data['canonical_smiles'],
            'class': clean_data['vs_activity_label'].astype(int),  # Ensure integer format
            'source': clean_data.get('source', 'VS_Integrated'),
            'vs_confidence': clean_data.get('vs_confidence_category', 'unknown'),
            'vs_tier': clean_data.get('vs_tier', 'unknown'),
            'num_data_sources': clean_data.get('num_data_sources', 1),
            'data_sources': clean_data.get('data_sources', 'VS_Integrated')
        })
        
        # Add any additional useful columns from the VS data
        optional_columns = [
            'vs_score', 'activity_ratio', 'mean_gi50', 'confidence', 
            'selectivity', 'tier', 'ml_weight', 'training_priority',
            'molecular_weight', 'logp', 'hbd', 'hba', 'tpsa'
        ]
        
        for col in optional_columns:
            if col in clean_data.columns:
                standardized_data[col] = clean_data[col]
        
        # Separate active and inactive compounds
        active_compounds = standardized_data[standardized_data['class'] == 1].copy()
        inactive_compounds = standardized_data[standardized_data['class'] == 0].copy()
        
        print(f"   Active compounds: {len(active_compounds)}")
        print(f"   Inactive compounds: {len(inactive_compounds)}")
        
        if len(active_compounds) == 0:
            print("❌ No active compounds found in dataset!")
            return False
        
        if len(inactive_compounds) == 0:
            print("⚠️ No inactive compounds found - decoy generator will create all decoys")
        
        # Create output file in the format expected by decoy generator
        output_file = 'compounds_across_celllines.xlsx'
        
        print(f"\n💾 Saving to decoy generator format: {output_file}")
        
        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            # Main detailed data (what the decoy generator loads first)
            standardized_data.to_excel(writer, sheet_name='Detailed_Data', index=False)
            
            # Separated active and inactive sheets (required by decoy generator)
            active_compounds.to_excel(writer, sheet_name='Active', index=False)
            inactive_compounds.to_excel(writer, sheet_name='Inactive', index=False)
            
            # Additional analysis sheets
            # Data source breakdown
            source_analysis = standardized_data.groupby(['source', 'class']).size().unstack(fill_value=0)
            source_analysis.columns = ['Inactive', 'Active']
            source_analysis['Total'] = source_analysis.sum(axis=1)
            source_analysis.reset_index().to_excel(writer, sheet_name='Source_Analysis', index=False)
            
            # Confidence analysis
            if 'vs_confidence' in standardized_data.columns:
                conf_analysis = standardized_data.groupby(['vs_confidence', 'class']).size().unstack(fill_value=0)
                conf_analysis.columns = ['Inactive', 'Active']
                conf_analysis['Total'] = conf_analysis.sum(axis=1)
                conf_analysis.reset_index().to_excel(writer, sheet_name='Confidence_Analysis', index=False)
            
            # Processing metadata
            metadata = pd.DataFrame({
                'Parameter': [
                    'Input_File',
                    'Processing_Date',
                    'Total_Compounds',
                    'Active_Compounds',
                    'Inactive_Compounds',
                    'Data_Sources',
                    'Activity_Ratio',
                    'Ready_for_Decoy_Generation'
                ],
                'Value': [
                    vs_input_file,
                    datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                    len(standardized_data),
                    len(active_compounds),
                    len(inactive_compounds),
                    ', '.join(standardized_data['source'].unique()),
                    f"1:{len(inactive_compounds)/len(active_compounds):.2f}" if len(active_compounds) > 0 else "N/A",
                    "YES"
                ]
            })
            metadata.to_excel(writer, sheet_name='Processing_Metadata', index=False)
        
        print(f"✅ Successfully converted VS data to decoy generator format!")
        print(f"📁 Output file: {output_file}")
        
        print(f"\n📊 CONVERSION SUMMARY:")
        print(f"   ✓ Total compounds processed: {len(standardized_data)}")
        print(f"   ✓ Active compounds: {len(active_compounds)}")
        print(f"   ✓ Inactive compounds: {len(inactive_compounds)}")
        print(f"   ✓ Data sources: {', '.join(standardized_data['source'].unique())}")
        
        print(f"\n🎯 NEXT STEPS:")
        print(f"   1. File '{output_file}' is now ready for the decoy generator")
        print(f"   2. Run the decoy generation script as normal")
        print(f"   3. The decoy generator will use your VS compounds as the basis")
        print(f"   4. Additional property-matched decoys will be generated")
        
        return True
        
    except Exception as e:
        print(f"❌ Error during conversion: {e}")
        import traceback
        traceback.print_exc()
        return False


def analyze_vs_data_for_decoy_suitability(vs_input_file='Integrated_VS_Dataset_Complete.xlsx'):
    """
    Analyze the VS data to assess its suitability for decoy generation
    """
    
    print("\n🔍 ANALYZING VS DATA FOR DECOY GENERATION SUITABILITY")
    print("=" * 60)
    
    try:
        vs_data = pd.read_excel(vs_input_file, sheet_name='Integrated_Dataset')
        
        print(f"📊 VS DATASET ANALYSIS:")
        print(f"   Total compounds: {len(vs_data)}")
        
        # Check activity labels
        if 'vs_activity_label' in vs_data.columns:
            activity_counts = vs_data['vs_activity_label'].value_counts()
            print(f"   Activity distribution:")
            print(f"     Active (1): {activity_counts.get(1, 0)}")
            print(f"     Inactive (0): {activity_counts.get(0, 0)}")
            print(f"     Unlabeled: {vs_data['vs_activity_label'].isna().sum()}")
        
        # Check SMILES availability
        valid_smiles = vs_data['canonical_smiles'].notna().sum()
        print(f"   Valid SMILES: {valid_smiles}/{len(vs_data)} ({valid_smiles/len(vs_data)*100:.1f}%)")
        
        # Check confidence levels
        if 'vs_confidence_category' in vs_data.columns:
            conf_counts = vs_data['vs_confidence_category'].value_counts()
            print(f"   Confidence levels:")
            for conf, count in conf_counts.items():
                print(f"     {conf}: {count}")
        
        # Check data sources
        if 'source' in vs_data.columns:
            source_counts = vs_data['source'].value_counts()
            print(f"   Data sources:")
            for source, count in source_counts.items():
                print(f"     {source}: {count}")
        
        # Assess suitability
        print(f"\n🎯 SUITABILITY ASSESSMENT:")
        
        suitable_compounds = vs_data[
            vs_data['canonical_smiles'].notna() & 
            vs_data['vs_activity_label'].notna()
        ]
        
        active_suitable = suitable_compounds[suitable_compounds['vs_activity_label'] == 1]
        inactive_suitable = suitable_compounds[suitable_compounds['vs_activity_label'] == 0]
        
        print(f"   Compounds suitable for decoy generation: {len(suitable_compounds)}")
        print(f"   Active compounds: {len(active_suitable)}")
        print(f"   Inactive compounds: {len(inactive_suitable)}")
        
        if len(active_suitable) >= 10:
            print("   ✅ Sufficient active compounds for decoy generation")
        else:
            print("   ⚠️ Limited active compounds - consider including more data")
        
        if len(inactive_suitable) >= 10:
            print("   ✅ Good baseline of inactive compounds")
        else:
            print("   ℹ️ Few inactive compounds - decoy generator will create additional ones")
        
        # Recommendations
        print(f"\n💡 RECOMMENDATIONS:")
        
        if len(suitable_compounds) >= 20:
            print("   ✅ Dataset is well-suited for decoy generation")
            print("   → Proceed with conversion and decoy generation")
        else:
            print("   ⚠️ Small dataset - consider:")
            print("   → Including more compounds from additional sources")
            print("   → Using a higher number of generated decoys")
        
        if len(active_suitable) > 0 and len(inactive_suitable) == 0:
            print("   ℹ️ Only active compounds available")
            print("   → Decoy generator will create all inactive compounds")
        
        return True
        
    except Exception as e:
        print(f"❌ Error during analysis: {e}")
        return False


def main():
    """
    Main function to bridge VS data to decoy generator format
    """
    
    print("🌉 VS DATA TO DECOY GENERATOR BRIDGE")
    print("=" * 40)
    print("This script converts your integrated VS dataset")
    print("into the format required by the decoy generator.")
    print("=" * 40)
    
    # Step 1: Analyze the VS data
    analysis_success = analyze_vs_data_for_decoy_suitability()
    
    if not analysis_success:
        print("\n❌ Could not analyze VS data")
        return False
    
    # Step 2: Convert the data
    conversion_success = convert_vs_data_for_decoy_generation()
    
    if conversion_success:
        print(f"\n🎉 SUCCESS! VS data converted for decoy generation")
        print(f"\n🚀 YOU CAN NOW RUN THE DECOY GENERATOR SCRIPT")
        print(f"   The file 'compounds_across_celllines.xlsx' is ready")
        print(f"   Simply run your decoy generation script as normal")
        
        print(f"\n📋 WHAT HAPPENS NEXT:")
        print(f"   1. Decoy generator will load your converted VS compounds")
        print(f"   2. It will calculate molecular properties")
        print(f"   3. It will generate additional property-matched decoys")
        print(f"   4. You'll get a balanced dataset for ML training")
        
        return True
    else:
        print(f"\n❌ Conversion failed")
        print(f"   Check error messages above for details")
        return False


if __name__ == "__main__":
    success = main()
    
    if success:
        print(f"\n✨ Bridge conversion completed successfully!")
        print(f"   Ready for decoy generation! 🚀")
    else:
        print(f"\n💥 Bridge conversion failed")
        print(f"   Check your input files and try again")