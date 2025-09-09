# -*- coding: utf-8 -*-
"""
Created on Tue Aug 19 16:45:30 2025

@author: ASUS TUF F15
"""

"""
DIAGNOSTIC DECOY GENERATOR FIX
Identifies and fixes SMILES processing issues
"""

import pandas as pd
import numpy as np
import os
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
import warnings
warnings.filterwarnings('ignore')

def diagnose_smiles_issues(compounds_file='compounds_across_celllines.xlsx'):
    """
    Diagnose SMILES processing issues in your data
    """
    print("🔍 DIAGNOSING SMILES PROCESSING ISSUES")
    print("=" * 45)
    
    try:
        # Load the data
        excel_file = pd.ExcelFile(compounds_file)
        print(f"Available sheets: {excel_file.sheet_names}")
        
        # Load active compounds
        active_data = pd.read_excel(compounds_file, sheet_name='Active')
        print(f"Active compounds loaded: {len(active_data)}")
        print(f"Columns: {list(active_data.columns)}")
        
        # Check SMILES column
        if 'canonical_smiles' not in active_data.columns:
            print("❌ No 'canonical_smiles' column found")
            return False
        
        # Analyze SMILES strings
        smiles_series = active_data['canonical_smiles']
        print(f"\n📊 SMILES Analysis:")
        print(f"   Total SMILES: {len(smiles_series)}")
        print(f"   Non-null SMILES: {smiles_series.notna().sum()}")
        print(f"   Empty strings: {(smiles_series == '').sum()}")
        print(f"   'null' strings: {(smiles_series == 'null').sum()}")
        
        # Test individual SMILES
        valid_smiles = []
        invalid_smiles = []
        
        print(f"\n🧪 Testing SMILES with RDKit:")
        for idx, smiles in enumerate(smiles_series.dropna()):
            if smiles and smiles != 'null' and smiles != '':
                try:
                    mol = Chem.MolFromSmiles(str(smiles))
                    if mol is not None:
                        valid_smiles.append(smiles)
                        print(f"   ✅ {idx+1}: Valid - {smiles[:50]}...")
                    else:
                        invalid_smiles.append(smiles)
                        print(f"   ❌ {idx+1}: Invalid - {smiles[:50]}...")
                except Exception as e:
                    invalid_smiles.append(smiles)
                    print(f"   ❌ {idx+1}: Error - {smiles[:50]}... ({e})")
            else:
                print(f"   ⚠️ {idx+1}: Empty/null SMILES")
        
        print(f"\n📈 SMILES Validation Results:")
        print(f"   Valid SMILES: {len(valid_smiles)}")
        print(f"   Invalid SMILES: {len(invalid_smiles)}")
        
        if len(invalid_smiles) > 0:
            print(f"\n❌ Invalid SMILES found:")
            for i, smiles in enumerate(invalid_smiles[:5]):  # Show first 5
                print(f"     {i+1}: {smiles}")
        
        # Test molecular property calculation on valid SMILES
        if len(valid_smiles) > 0:
            print(f"\n🧪 Testing molecular property calculation:")
            test_smiles = valid_smiles[0]
            try:
                mol = Chem.MolFromSmiles(test_smiles)
                mw = Descriptors.MolWt(mol)
                logp = Descriptors.MolLogP(mol)
                print(f"   ✅ Successfully calculated properties for: {test_smiles[:30]}...")
                print(f"      MW: {mw:.2f}, LogP: {logp:.2f}")
                return True
            except Exception as e:
                print(f"   ❌ Property calculation failed: {e}")
                return False
        else:
            print(f"   ❌ No valid SMILES found for testing")
            return False
            
    except Exception as e:
        print(f"❌ Error during diagnosis: {e}")
        return False

def fix_smiles_data(compounds_file='compounds_across_celllines.xlsx'):
    """
    Fix SMILES data issues and create a cleaned version
    """
    print(f"\n🔧 FIXING SMILES DATA")
    print("=" * 25)
    
    try:
        # Load all sheets
        excel_file = pd.ExcelFile(compounds_file)
        
        fixed_sheets = {}
        
        for sheet_name in ['Active', 'Inactive', 'Detailed_Data']:
            if sheet_name in excel_file.sheet_names:
                print(f"Processing {sheet_name} sheet...")
                
                df = pd.read_excel(compounds_file, sheet_name=sheet_name)
                original_count = len(df)
                
                # Clean SMILES
                if 'canonical_smiles' in df.columns:
                    # Remove null, empty, and 'null' strings
                    df = df[df['canonical_smiles'].notna()]
                    df = df[df['canonical_smiles'] != '']
                    df = df[df['canonical_smiles'] != 'null']
                    
                    # Validate remaining SMILES
                    valid_rows = []
                    for idx, row in df.iterrows():
                        smiles = row['canonical_smiles']
                        try:
                            mol = Chem.MolFromSmiles(str(smiles))
                            if mol is not None:
                                valid_rows.append(row)
                        except:
                            continue
                    
                    if valid_rows:
                        cleaned_df = pd.DataFrame(valid_rows)
                        fixed_sheets[sheet_name] = cleaned_df
                        print(f"   {sheet_name}: {original_count} → {len(cleaned_df)} compounds")
                    else:
                        print(f"   ❌ No valid compounds in {sheet_name}")
                else:
                    print(f"   ⚠️ No canonical_smiles column in {sheet_name}")
        
        # Save fixed data
        if fixed_sheets:
            output_file = 'compounds_across_celllines_FIXED.xlsx'
            with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
                for sheet_name, df in fixed_sheets.items():
                    df.to_excel(writer, sheet_name=sheet_name, index=False)
            
            print(f"✅ Fixed data saved to: {output_file}")
            
            # Summary
            print(f"\n📊 FIXED DATA SUMMARY:")
            for sheet_name, df in fixed_sheets.items():
                print(f"   {sheet_name}: {len(df)} valid compounds")
            
            return output_file
        else:
            print(f"❌ No valid data to save")
            return None
            
    except Exception as e:
        print(f"❌ Error fixing data: {e}")
        return None

def create_simple_decoy_generator(fixed_file):
    """
    Simple decoy generator that works with the fixed data
    """
    print(f"\n🏗️ CREATING SIMPLE DECOY GENERATOR")
    print("=" * 35)
    
    try:
        # Load fixed data
        active_data = pd.read_excel(fixed_file, sheet_name='Active')
        
        if 'Inactive' in pd.ExcelFile(fixed_file).sheet_names:
            inactive_data = pd.read_excel(fixed_file, sheet_name='Inactive')
        else:
            inactive_data = pd.DataFrame()
        
        print(f"Active compounds: {len(active_data)}")
        print(f"Inactive compounds: {len(inactive_data)}")
        
        # Calculate properties for active compounds
        print(f"\nCalculating properties for active compounds...")
        
        active_properties = []
        for idx, row in active_data.iterrows():
            try:
                mol = Chem.MolFromSmiles(row['canonical_smiles'])
                if mol is not None:
                    props = {
                        'compound_id': row['molregno_chembl_id'],
                        'smiles': row['canonical_smiles'],
                        'activity': 1,
                        'MW': Descriptors.MolWt(mol),
                        'LogP': Descriptors.MolLogP(mol),
                        'HBD': Descriptors.NumHDonors(mol),
                        'HBA': Descriptors.NumHAcceptors(mol),
                        'TPSA': Descriptors.TPSA(mol),
                        'RotBonds': Descriptors.NumRotatableBonds(mol)
                    }
                    active_properties.append(props)
            except Exception as e:
                print(f"   Error with compound {idx}: {e}")
                continue
        
        print(f"   Successfully processed: {len(active_properties)} active compounds")
        
        # Process inactive compounds if available
        inactive_properties = []
        if len(inactive_data) > 0:
            print(f"\nCalculating properties for inactive compounds...")
            for idx, row in inactive_data.iterrows():
                try:
                    mol = Chem.MolFromSmiles(row['canonical_smiles'])
                    if mol is not None:
                        props = {
                            'compound_id': row['molregno_chembl_id'],
                            'smiles': row['canonical_smiles'],
                            'activity': 0,
                            'MW': Descriptors.MolWt(mol),
                            'LogP': Descriptors.MolLogP(mol),
                            'HBD': Descriptors.NumHDonors(mol),
                            'HBA': Descriptors.NumHAcceptors(mol),
                            'TPSA': Descriptors.TPSA(mol),
                            'RotBonds': Descriptors.NumRotatableBonds(mol)
                        }
                        inactive_properties.append(props)
                except Exception as e:
                    continue
            
            print(f"   Successfully processed: {len(inactive_properties)} inactive compounds")
        
        # Generate simple decoys
        print(f"\nGenerating simple decoys...")
        
        if len(active_properties) > 0:
            # Calculate property ranges from active compounds
            active_df = pd.DataFrame(active_properties)
            
            mw_range = (active_df['MW'].min(), active_df['MW'].max())
            logp_range = (active_df['LogP'].min(), active_df['LogP'].max())
            
            print(f"   MW range: {mw_range[0]:.1f} - {mw_range[1]:.1f}")
            print(f"   LogP range: {logp_range[0]:.1f} - {logp_range[1]:.1f}")
            
            # Generate simple decoys
            simple_decoys = []
            decoy_smiles = [
                'CCCCCCCC', 'CCCCCCCCC', 'CCCCCCCCCC',  # Alkanes
                'c1ccccc1', 'c1ccc2ccccc2c1',  # Aromatics
                'CCO', 'CCCO', 'CCCCO',  # Alcohols
                'CCC(=O)O', 'CCCC(=O)O',  # Acids
                'c1ccc(cc1)C', 'c1ccc(cc1)CC'  # Substituted benzenes
            ]
            
            for i, smiles in enumerate(decoy_smiles):
                try:
                    mol = Chem.MolFromSmiles(smiles)
                    if mol is not None:
                        mw = Descriptors.MolWt(mol)
                        logp = Descriptors.MolLogP(mol)
                        
                        # Check if within reasonable range
                        if (mw_range[0] * 0.5 <= mw <= mw_range[1] * 1.5 and
                            logp_range[0] - 2 <= logp <= logp_range[1] + 2):
                            
                            decoy = {
                                'compound_id': f'DECOY_{i+1:03d}',
                                'smiles': smiles,
                                'activity': 0,
                                'MW': mw,
                                'LogP': logp,
                                'HBD': Descriptors.NumHDonors(mol),
                                'HBA': Descriptors.NumHAcceptors(mol),
                                'TPSA': Descriptors.TPSA(mol),
                                'RotBonds': Descriptors.NumRotatableBonds(mol)
                            }
                            simple_decoys.append(decoy)
                except:
                    continue
            
            print(f"   Generated: {len(simple_decoys)} simple decoys")
        
        # Combine all data
        all_data = active_properties + inactive_properties + simple_decoys
        final_df = pd.DataFrame(all_data)
        
        # Save results
        output_file = 'Simple_Balanced_Dataset.xlsx'
        
        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            final_df.to_excel(writer, sheet_name='Complete_Dataset', index=False)
            
            # Separate sheets
            active_df = final_df[final_df['activity'] == 1]
            inactive_df = final_df[final_df['activity'] == 0]
            
            active_df.to_excel(writer, sheet_name='Active_Compounds', index=False)
            inactive_df.to_excel(writer, sheet_name='Inactive_Compounds', index=False)
            
            # Summary
            summary = pd.DataFrame({
                'Metric': ['Total_Compounds', 'Active_Compounds', 'Inactive_Compounds', 'Ratio_A:I'],
                'Value': [len(final_df), len(active_df), len(inactive_df), 
                         f"1:{len(inactive_df)/len(active_df):.1f}" if len(active_df) > 0 else "N/A"]
            })
            summary.to_excel(writer, sheet_name='Summary', index=False)
        
        print(f"\n✅ Simple balanced dataset created: {output_file}")
        print(f"   Total compounds: {len(final_df)}")
        print(f"   Active: {len(active_df)}")
        print(f"   Inactive: {len(inactive_df)}")
        
        return output_file
        
    except Exception as e:
        print(f"❌ Error creating simple dataset: {e}")
        import traceback
        traceback.print_exc()
        return None

def main():
    """
    Main diagnostic and fix workflow
    """
    print("🔧 DIAGNOSTIC DECOY GENERATOR FIX")
    print("=" * 35)
    print("This script will:")
    print("1. Diagnose SMILES processing issues")
    print("2. Fix invalid SMILES data")
    print("3. Create a simple balanced dataset")
    print("=" * 35)
    
    # Step 1: Diagnose issues
    compounds_file = 'compounds_across_celllines.xlsx'
    
    print(f"\n📋 Step 1: Diagnosing issues...")
    diagnosis_success = diagnose_smiles_issues(compounds_file)
    
    if not diagnosis_success:
        print(f"❌ Diagnosis failed - check your input file")
        return False
    
    # Step 2: Fix data
    print(f"\n📋 Step 2: Fixing data...")
    fixed_file = fix_smiles_data(compounds_file)
    
    if fixed_file is None:
        print(f"❌ Could not fix data")
        return False
    
    # Step 3: Create simple dataset
    print(f"\n📋 Step 3: Creating simple balanced dataset...")
    output_file = create_simple_decoy_generator(fixed_file)
    
    if output_file:
        print(f"\n🎉 SUCCESS!")
        print(f"📁 Fixed data: {fixed_file}")
        print(f"📁 Balanced dataset: {output_file}")
        print(f"\n💡 You now have a working balanced dataset for ML training!")
        return True
    else:
        print(f"❌ Failed to create dataset")
        return False

if __name__ == "__main__":
    success = main()
    
    if success:
        print(f"\n✨ Diagnostic fix completed successfully!")
    else:
        print(f"\n💥 Diagnostic fix failed")