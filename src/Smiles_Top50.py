# -*- coding: utf-8 -*-
"""
Created on Mon Sep  8 00:10:02 2025

@author: ASUS TUF F15
"""

# -*- coding: utf-8 -*-
"""
SMILES Data Adder for Top50 Sheets
Adds SMILES columns directly to COCONUT_Top50, FOODB_Top50, and LOTUS_Top50 sheets
"""

import pandas as pd
import numpy as np
import glob
import os
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

class SMILESAdder:
    """Add SMILES data to Top50 sheets"""
    
    def __init__(self):
        self.smiles_databases = {}
        
    def load_original_databases(self):
        """Load SMILES data from original database files"""
        print("🔍 LOADING ORIGINAL DATABASE FILES")
        print("=" * 40)
        
        databases = ['COCONUT', 'LOTUS', 'FOODB']
        
        for db_name in databases:
            print(f"\n📊 Loading {db_name} database...")
            
            try:
                # Find database files
                if db_name == 'COCONUT':
                    files = glob.glob("*COCONUT*.csv")
                elif db_name == 'LOTUS':
                    files = glob.glob("*LOTUS*.csv")
                elif db_name == 'FOODB':
                    files = glob.glob("*FOODB*.csv")
                
                if not files:
                    print(f"   ❌ No {db_name} database file found")
                    continue
                
                # Use the largest file (likely the most complete)
                db_file = max(files, key=os.path.getsize)
                print(f"   📁 Using file: {db_file}")
                
                # Try different encodings
                df = None
                for encoding in ['utf-8', 'latin-1', 'cp1252', 'iso-8859-1']:
                    try:
                        df = pd.read_csv(db_file, encoding=encoding, nrows=50000)  # Limit for performance
                        print(f"   ✅ Loaded with {encoding} encoding")
                        break
                    except:
                        continue
                
                if df is None:
                    print(f"   ❌ Could not load {db_file} with any encoding")
                    continue
                
                print(f"   📏 Shape: {df.shape[0]:,} rows × {df.shape[1]:,} columns")
                print(f"   📋 Columns: {list(df.columns[:5])}...")
                
                # Find SMILES and ID columns
                smiles_cols = [col for col in df.columns if 'smiles' in col.lower()]
                id_cols = [col for col in df.columns if any(x in col.lower() for x in ['id', 'identifier'])]
                
                print(f"   🧪 SMILES columns found: {smiles_cols}")
                print(f"   🆔 ID columns found: {id_cols}")
                
                if not smiles_cols or not id_cols:
                    print(f"   ❌ Missing required columns for {db_name}")
                    continue
                
                # Extract ID-SMILES mapping
                smiles_col = smiles_cols[0]
                id_col = id_cols[0]
                
                id_smiles_map = df[[id_col, smiles_col]].copy()
                id_smiles_map.columns = ['ID', 'SMILES']
                id_smiles_map = id_smiles_map.dropna()
                
                # Clean up IDs (remove .0 suffixes if they exist)
                id_smiles_map['ID'] = id_smiles_map['ID'].astype(str).str.replace(r'\.0$', '', regex=True)
                
                self.smiles_databases[db_name] = id_smiles_map
                print(f"   ✅ Stored {len(id_smiles_map):,} ID-SMILES pairs for {db_name}")
                
            except Exception as e:
                print(f"   ❌ Error loading {db_name}: {e}")
        
        print(f"\n📊 SUMMARY: Loaded {len(self.smiles_databases)} databases")
        for db, data in self.smiles_databases.items():
            print(f"   {db}: {len(data):,} compounds")
    
    def add_smiles_to_top50_file(self, input_file='Top_50_Results.xlsx'):
        """Add SMILES columns to the Top50 Excel file"""
        print(f"\n🔬 ADDING SMILES TO TOP50 SHEETS")
        print("=" * 35)
        
        # Check if input file exists
        if not os.path.exists(input_file):
            print(f"❌ File not found: {input_file}")
            return None
        
        # Load the Excel file
        try:
            excel_data = pd.read_excel(input_file, sheet_name=None)  # Load all sheets
            print(f"✅ Loaded {input_file}")
            print(f"📋 Sheets found: {list(excel_data.keys())}")
        except Exception as e:
            print(f"❌ Error loading Excel file: {e}")
            return None
        
        # Process each sheet
        updated_sheets = {}
        smiles_added_summary = {}
        
        for sheet_name, df in excel_data.items():
            print(f"\n📊 Processing sheet: {sheet_name}")
            print(f"   📏 Original size: {len(df):,} rows × {len(df.columns):,} columns")
            
            # Determine database from sheet name
            database = None
            for db in ['COCONUT', 'LOTUS', 'FOODB']:
                if db in sheet_name.upper():
                    database = db
                    break
            
            if database is None:
                print(f"   ⚠️ Cannot determine database for {sheet_name}, skipping")
                updated_sheets[sheet_name] = df
                continue
            
            print(f"   🗄️ Database identified: {database}")
            
            # Check if we have SMILES data for this database
            if database not in self.smiles_databases:
                print(f"   ❌ No SMILES data available for {database}")
                updated_sheets[sheet_name] = df
                smiles_added_summary[sheet_name] = 0
                continue
            
            # Find ID column in the sheet
            id_columns = [col for col in df.columns if any(x in str(col).lower() for x in ['id', 'identifier', 'compound_id'])]
            
            if not id_columns:
                print(f"   ❌ No ID column found in {sheet_name}")
                print(f"   📋 Available columns: {list(df.columns)}")
                updated_sheets[sheet_name] = df
                smiles_added_summary[sheet_name] = 0
                continue
            
            id_col = id_columns[0]
            print(f"   🆔 Using ID column: {id_col}")
            
            # Check if SMILES column already exists
            existing_smiles = [col for col in df.columns if 'smiles' in str(col).lower()]
            if existing_smiles:
                print(f"   ✅ SMILES column already exists: {existing_smiles[0]}")
                updated_sheets[sheet_name] = df
                smiles_added_summary[sheet_name] = df[existing_smiles[0]].notna().sum()
                continue
            
            # Prepare compound IDs for matching
            df_copy = df.copy()
            df_copy['Clean_ID'] = df_copy[id_col].astype(str).str.replace(r'\.0$', '', regex=True)
            
            # Get SMILES database for this database
            smiles_db = self.smiles_databases[database]
            
            # Create lookup dictionary for faster matching
            smiles_lookup = dict(zip(smiles_db['ID'], smiles_db['SMILES']))
            
            # Add SMILES column
            df_copy['SMILES'] = df_copy['Clean_ID'].map(smiles_lookup)
            
            # Drop the temporary Clean_ID column
            df_copy = df_copy.drop('Clean_ID', axis=1)
            
            # Count successful matches
            smiles_found = df_copy['SMILES'].notna().sum()
            smiles_added_summary[sheet_name] = smiles_found
            
            print(f"   ✅ Added SMILES for {smiles_found:,}/{len(df):,} compounds ({smiles_found/len(df)*100:.1f}%)")
            
            # Show some examples
            if smiles_found > 0:
                examples = df_copy[df_copy['SMILES'].notna()].head(3)
                print(f"   📝 Examples:")
                for _, row in examples.iterrows():
                    compound_name = row.get('Compound_Name', 'Unknown')
                    smiles = row['SMILES']
                    print(f"      {compound_name}: {smiles[:50]}{'...' if len(smiles) > 50 else ''}")
            
            updated_sheets[sheet_name] = df_copy
        
        # Save updated file
        output_file = input_file.replace('.xlsx', '_with_SMILES.xlsx')
        
        try:
            with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
                for sheet_name, df in updated_sheets.items():
                    df.to_excel(writer, sheet_name=sheet_name, index=False)
            
            print(f"\n💾 SAVED UPDATED FILE")
            print("=" * 25)
            print(f"✅ Output file: {output_file}")
            
            # Summary
            total_original = sum(len(df) for df in excel_data.values())
            total_smiles = sum(smiles_added_summary.values())
            
            print(f"\n📊 SMILES ADDITION SUMMARY:")
            print(f"   Total compounds: {total_original:,}")
            print(f"   SMILES added: {total_smiles:,} ({total_smiles/total_original*100:.1f}%)")
            
            for sheet_name, count in smiles_added_summary.items():
                original_count = len(excel_data[sheet_name])
                print(f"   {sheet_name}: {count:,}/{original_count:,} ({count/original_count*100:.1f}%)")
            
            return output_file
            
        except Exception as e:
            print(f"❌ Error saving file: {e}")
            return None
    
    def create_smiles_lookup_files(self):
        """Create separate CSV files with ID-SMILES mappings for reference"""
        print(f"\n📁 CREATING SMILES LOOKUP FILES")
        print("=" * 35)
        
        for db_name, smiles_data in self.smiles_databases.items():
            output_file = f"{db_name}_ID_SMILES_Lookup.csv"
            
            try:
                smiles_data.to_csv(output_file, index=False)
                print(f"✅ Created {output_file} ({len(smiles_data):,} entries)")
            except Exception as e:
                print(f"❌ Error creating {output_file}: {e}")

def main():
    """Main function to add SMILES to Top50 sheets"""
    print("🧪 SMILES DATA ADDER FOR TOP50 SHEETS")
    print("=" * 40)
    print("This script will add SMILES columns to your Top50 sheets")
    print()
    
    # Initialize SMILES adder
    smiles_adder = SMILESAdder()
    
    # Load original databases
    smiles_adder.load_original_databases()
    
    if not smiles_adder.smiles_databases:
        print("❌ No database files found! Please ensure you have:")
        print("   • COCONUT database CSV file")
        print("   • LOTUS database CSV file") 
        print("   • FOODB database CSV file")
        return
    
    # Add SMILES to Top50 file
    output_file = smiles_adder.add_smiles_to_top50_file()
    
    if output_file:
        print(f"\n🎉 SUCCESS!")
        print("=" * 15)
        print(f"✅ SMILES data added to Top50 sheets")
        print(f"✅ Updated file: {output_file}")
        print("✅ Ready for toxicity screening!")
        
        # Create lookup files for reference
        smiles_adder.create_smiles_lookup_files()
        
        print(f"\n🎯 NEXT STEPS:")
        print("1. Use the new file with SMILES data for toxicity screening")
        print("2. Run the toxicity screening script with the updated file")
        print("3. Expect much faster processing with ~150 compounds total")
    
    else:
        print("\n❌ SMILES addition failed!")
        print("Please check error messages above and try again")

if __name__ == "__main__":
    main()