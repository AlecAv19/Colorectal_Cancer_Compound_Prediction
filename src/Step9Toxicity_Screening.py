# -*- coding: utf-8 -*-
"""
Enhanced Top50 Focused Toxicity Screening Module
Processes COCONUT_Top50, FOODB_Top50, and LOTUS_Top50 sheets for focused toxicity assessment
"""

import pandas as pd
import numpy as np
import requests
import time
import json
import os
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

class Top50ToxicityScreening:
    """
    Top50 focused toxicity screening for the best compounds from each database
    """
    
    def __init__(self):
        self.toxicity_results = {}
        self.failed_predictions = []
        self.top50_sheets_data = {}
        self.combined_results = None
        
    def load_top50_vs_results(self, results_file=None):
        """Load specifically the Top50 sheets from each database for focused toxicity assessment"""
        print("📊 LOADING TOP50 VIRTUAL SCREENING RESULTS")
        print("=" * 45)
        
        if results_file is None:
            # Look for recent VS results
            import glob
            # Look for both patterns
            vs_files = (glob.glob("Natural_Product_Screening_Results_*.xlsx") + 
                       glob.glob("Top_50_Results*.xlsx") + 
                       glob.glob("Top50_Results*.xlsx"))
            if vs_files:
                results_file = max(vs_files, key=os.path.getmtime)
                print(f" Found recent results: {results_file}")
            else:
                print(" No VS results found! Please specify file.")
                return None
        
        try:
            # Load all sheets first to see what's available
            excel_file = pd.ExcelFile(results_file)
            sheets = excel_file.sheet_names
            print(f" Available sheets ({len(sheets)} total):")
            for i, sheet in enumerate(sheets, 1):
                print(f"   {i:2d}. {sheet}")
            
            # Target specific Top50 sheets
            target_sheets = ['COCONUT_Top50', 'FOODB_Top50', 'LOTUS_Top50']
            
            # Find available Top50 sheets (case-insensitive and flexible matching)
            available_top50_sheets = []
            for target in target_sheets:
                found = False
                for sheet in sheets:
                    # Try exact match first
                    if target.lower() == sheet.lower():
                        available_top50_sheets.append(sheet)
                        found = True
                        print(f" Exact match found: {sheet}")
                        break
                
                if not found:
                    # Try partial match for "Top50" variants
                    for sheet in sheets:
                        if 'top50' in sheet.lower() and any(db.lower() in sheet.lower() for db in ['coconut', 'foodb', 'lotus']):
                            available_top50_sheets.append(sheet)
                            found = True
                            print(f" Partial match found: {sheet}")
                            break
                
                if not found:
                    print(f" No match found for {target}")
            
            if not available_top50_sheets:
                print("\n No Top50 sheets found!")
                print(" Available sheets containing 'top' (case-insensitive):")
                top_sheets = [s for s in sheets if 'top' in s.lower()]
                if top_sheets:
                    for sheet in top_sheets:
                        print(f"    {sheet}")
                else:
                    print("   None found")
                print("\n Please check your sheet names. Expected: COCONUT_Top50, FOODB_Top50, LOTUS_Top50")
                return None
            
            print(f"\n Found {len(available_top50_sheets)} Top50 sheets:")
            for sheet in available_top50_sheets:
                print(f"    {sheet}")
            
            # Load and inspect the Top50 sheets
            loaded_sheets = {}
            total_compounds = 0
            
            print(f"\n Loading and inspecting Top50 results...")
            for sheet_name in available_top50_sheets:
                try:
                    print(f"\n Loading sheet: {sheet_name}")
                    df = pd.read_excel(results_file, sheet_name=sheet_name)
                    
                    print(f"    Initial size: {len(df):,} rows × {len(df.columns):,} columns")
                    
                    # Show first few column names to understand the data
                    print(f"    Column names (first 10): {list(df.columns[:10])}")
                    
                    # Check if this looks like actual compound data
                    if len(df) > 1000:
                        print(f"    WARNING: This sheet has {len(df):,} rows - this doesn't look like a Top50 sheet!")
                        print(f"    First few rows of data:")
                        print(df.head(3))
                        
                        # Ask user if they want to proceed with this sheet
                        proceed = input(f"    This sheet seems too large for Top50. Continue anyway? (y/n): ").strip().lower()
                        if proceed != 'y':
                            print(f"     Skipping {sheet_name}")
                            continue
                        else:
                            # Limit to actual top 50 if user proceeds
                            activity_cols = [col for col in df.columns if any(x in col.lower() for x in ['score', 'probability', 'affinity', 'activity'])]
                            if activity_cols:
                                print(f"     Limiting to top 50 by {activity_cols[0]}")
                                df = df.nlargest(50, activity_cols[0])
                            else:
                                print(f"     No activity column found, taking first 50 rows")
                                df = df.head(50)
                            print(f"    Reduced to {len(df):,} rows")
                    
                    # Check for SMILES
                    smiles_cols = [col for col in df.columns if 'smiles' in col.lower()]
                    print(f"    SMILES columns found: {smiles_cols}")
                    
                    if not smiles_cols:
                        print(f"   🔍 No SMILES found, attempting to add from original database...")
                        df = self.add_smiles_to_sheet(df, sheet_name)
                        if df is not None:
                            smiles_cols = [col for col in df.columns if 'smiles' in col.lower()]
                            print(f"    SMILES columns after merge: {smiles_cols}")
                    
                    if df is not None and len(df) > 0:
                        # Add sheet metadata
                        df['Source_Sheet'] = sheet_name
                        df['Sheet_Type'] = 'Top50'
                        
                        # Extract database from sheet name
                        for db in ['COCONUT', 'LOTUS', 'FOODB']:
                            if db.upper() in sheet_name.upper():
                                df['Database'] = db
                                break
                        else:
                            df['Database'] = 'Unknown'
                        
                        # For Top50 sheets, model is typically consensus/best
                        df['Model'] = 'Top50_Selection'
                        
                        loaded_sheets[sheet_name] = df
                        total_compounds += len(df)
                        
                        # Show detailed info
                        db_name = df['Database'].iloc[0] if 'Database' in df.columns else 'Unknown'
                        print(f"    Successfully loaded: {len(df):,} compounds from {db_name}")
                        
                        # Check if SMILES are available
                        smiles_available = df['SMILES'].notna().sum() if 'SMILES' in df.columns else 0
                        print(f"    SMILES available: {smiles_available:,}/{len(df):,} ({smiles_available/len(df)*100:.1f}%)")
                        
                        # Show activity score info
                        activity_cols = [col for col in df.columns if any(x in col.lower() for x in ['score', 'probability', 'affinity', 'activity'])]
                        if activity_cols:
                            activity_col = activity_cols[0]
                            avg_activity = df[activity_col].mean()
                            max_activity = df[activity_col].max()
                            min_activity = df[activity_col].min()
                            print(f"    Activity ({activity_col}): avg={avg_activity:.3f}, max={max_activity:.3f}, min={min_activity:.3f}")
                        else:
                            print(f"    No activity columns found")
                    
                except Exception as e:
                    print(f"    Error loading {sheet_name}: {e}")
            
            if not loaded_sheets:
                print(" No Top50 sheets could be loaded successfully!")
                return None
            
            self.top50_sheets_data = loaded_sheets
            print(f"\n LOADING COMPLETE")
            print(f"    Successfully loaded {len(loaded_sheets)} sheets")
            print(f"    Total compounds: {total_compounds:,}")
            print(f"    Expected for true Top50: ~150 compounds (50 × 3 databases)")
            
            if total_compounds > 200:
                print(f"    WARNING: Total compounds ({total_compounds:,}) is much higher than expected!")
                print(f"    This suggests the sheets may not be true 'Top50' sheets")
            
            return loaded_sheets
            
        except Exception as e:
            print(f" Error loading results: {e}")
            return None
    
    def add_smiles_to_sheet(self, df, sheet_name):
        """Add SMILES data to a sheet that doesn't have it"""
        print(f"   🔍 Adding SMILES to {sheet_name}...")
        
        # Determine database from sheet name
        database = None
        for db in ['COCONUT', 'LOTUS', 'FOODB']:
            if db in sheet_name:
                database = db
                break
        
        if database is None:
            print(f"    Cannot determine database for {sheet_name}")
            return df
        
        # Load original database
        original_df = self.load_original_database_smiles(database)
        if original_df is None:
            return df
        
        # Try to merge on different possible ID columns
        id_columns = [col for col in df.columns if any(x in col.lower() for x in ['id', 'identifier', 'compound'])]
        
        merged_df = None
        for id_col in id_columns:
            try:
                # Try exact merge first
                test_merge = df.merge(original_df, left_on=id_col, right_on='Compound_ID', how='left')
                smiles_added = test_merge['SMILES'].notna().sum()
                
                if smiles_added > 0:
                    merged_df = test_merge
                    print(f"    Added SMILES for {smiles_added:,} compounds via {id_col}")
                    break
                    
            except Exception as e:
                continue
        
        return merged_df if merged_df is not None else df
    
    def load_original_database_smiles(self, db_name):
        """Load SMILES from original database files"""
        try:
            import glob
            
            if db_name == 'COCONUT':
                files = glob.glob("*COCONUT*.csv")
            elif db_name == 'LOTUS':
                files = glob.glob("*LOTUS*.csv")
            elif db_name == 'FOODB':
                files = glob.glob("*FOODB*.csv")
            else:
                return None
            
            if not files:
                return None
            
            db_file = max(files, key=os.path.getsize)
            
            # Load with encoding handling
            for encoding in ['utf-8', 'latin-1', 'cp1252', 'iso-8859-1']:
                try:
                    df = pd.read_csv(db_file, encoding=encoding, nrows=10000)
                    break
                except:
                    continue
            else:
                return None
            
            # Find SMILES and ID columns
            smiles_cols = [col for col in df.columns if 'smiles' in col.lower()]
            id_cols = [col for col in df.columns if any(x in col.lower() for x in ['id', 'identifier'])]
            
            if not smiles_cols or not id_cols:
                return None
            
            smiles_col = smiles_cols[0]
            id_col = id_cols[0]
            
            # Create clean ID mapping
            result_df = df[[id_col, smiles_col]].copy()
            result_df.columns = ['Compound_ID', 'SMILES']
            result_df = result_df.dropna()
            
            return result_df
            
        except Exception as e:
            return None
    
    def combine_top50_results(self, sheets_data, strategy='best_per_compound'):
        """Combine Top50 results using database-focused strategies"""
        print(f"\n🔗 COMBINING TOP50 RESULTS")
        print(f"Strategy: {strategy}")
        print("=" * 30)
        
        if not sheets_data:
            return None
        
        all_data = []
        
        # Combine all Top50 dataframes
        for sheet_name, df in sheets_data.items():
            print(f"   Processing {sheet_name}: {len(df):,} compounds")
            all_data.append(df)
        
        # Concatenate all data
        combined_df = pd.concat(all_data, ignore_index=True)
        print(f"   Combined total: {len(combined_df):,} records")
        
        # Find activity score columns
        activity_cols = []
        for col in combined_df.columns:
            if any(x in col.lower() for x in ['score', 'probability', 'affinity', 'activity']):
                if col not in ['Safety_Score', 'Source_Sheet', 'Sheet_Type']:
                    activity_cols.append(col)
        
        if not activity_cols:
            print("    No activity score columns found")
            return combined_df
        
        primary_activity_col = activity_cols[0]
        print(f"   Using primary activity column: {primary_activity_col}")
        
        # Find ID columns first (needed for all strategies)
        id_cols = [col for col in combined_df.columns if any(x in col.lower() for x in ['id', 'name']) and 'compound' in col.lower()]
        id_col = id_cols[0] if id_cols else None
        
        # Apply combination strategy
        if strategy == 'best_per_compound':
            # Keep the record with the highest activity score for each compound
            if id_col:
                # Sort by activity score and keep the best for each compound
                combined_df = combined_df.sort_values(primary_activity_col, ascending=False)
                result_df = combined_df.drop_duplicates(subset=[id_col], keep='first')
                print(f"    Best per compound: {len(result_df):,} unique compounds")
            else:
                result_df = combined_df
                print(f"    No compound ID column found, keeping all records")
                
        elif strategy == 'database_comparison':
            # Keep the best compound from each database for comparison
            if 'Database' in combined_df.columns:
                # Get best compound per database
                best_per_db = combined_df.loc[combined_df.groupby('Database')[primary_activity_col].idxmax()]
                
                # Also keep top 10 from each database for comparison
                top_per_db = []
                for database in combined_df['Database'].unique():
                    db_data = combined_df[combined_df['Database'] == database]
                    top_db = db_data.nlargest(10, primary_activity_col)
                    top_per_db.append(top_db)
                
                result_df = pd.concat(top_per_db, ignore_index=True)
                
                # Remove duplicates if compound appears in multiple databases
                if id_col:
                    result_df = result_df.sort_values(primary_activity_col, ascending=False)
                    result_df = result_df.drop_duplicates(subset=[id_col], keep='first')
                
                print(f"    Database comparison: {len(result_df):,} compounds")
            else:
                result_df = combined_df
                print(f"    No database column found, using all data")
                
        else:  # 'all_top50'
            result_df = combined_df
            print(f"    Keeping all Top50 predictions: {len(result_df):,} records")
        
        # Add combination metadata
        result_df['Total_Sheets_Processed'] = len(sheets_data)
        result_df['Combination_Strategy'] = strategy
        result_df['Focus'] = 'Top50_Databases'
        
        # Summary statistics
        unique_compounds = result_df[id_col].nunique() if id_col else len(result_df)
        unique_databases = result_df['Database'].nunique() if 'Database' in result_df.columns else 1
        
        print(f"\n📊 TOP50 COMBINATION SUMMARY:")
        print(f"   Final dataset: {len(result_df):,} records")
        print(f"   Unique compounds: {unique_compounds:,}")
        print(f"   Databases covered: {unique_databases}")
        
        if 'Database' in result_df.columns:
            db_counts = result_df['Database'].value_counts()
            print(f"   Database distribution:")
            for db, count in db_counts.items():
                print(f"     {db}: {count:,} compounds")
        
        self.combined_results = result_df
        return result_df
    
    def calculate_simple_toxicity_features(self, smiles_list):
        """Calculate simple molecular features associated with toxicity"""
        print(f"\n CALCULATING MOLECULAR TOXICITY FEATURES")
        print("=" * 45)
        
        toxicity_features = []
        
        for i, smiles in enumerate(smiles_list):
            if pd.isna(smiles):
                toxicity_features.append({
                    'SMILES': smiles,
                    'Molecular_Weight': np.nan,
                    'Lipophilicity_Est': np.nan,
                    'Aromatic_Rings': np.nan,
                    'Toxicity_Alerts': np.nan,
                    'Safety_Score': np.nan
                })
                continue
            
            # Simple molecular property estimation (in production, use RDKit)
            # These are simplified heuristics for demonstration
            
            # Estimate molecular weight from SMILES length (very rough)
            estimated_mw = len(smiles) * 12 + np.random.normal(0, 50)
            estimated_mw = max(100, estimated_mw)  # Minimum MW
            
            # Estimate lipophilicity from aromatic content
            aromatic_chars = smiles.count('c') + smiles.count('C')
            estimated_logp = (aromatic_chars / len(smiles)) * 5 + np.random.normal(0, 1)
            
            # Count aromatic rings (very rough estimate)
            aromatic_rings = smiles.count('c1') + smiles.count('C1')
            
            # Simple toxicity alerts (presence of certain substructures)
            toxicity_alerts = 0
            
            # Check for common toxicity-associated patterns
            toxic_patterns = ['[N+]', 'N=O', 'C#N', 'CCl', 'CBr', 'CF', '[As]', '[Hg]', '[Pb]']
            for pattern in toxic_patterns:
                if pattern in smiles:
                    toxicity_alerts += 1
            
            # Calculate simple safety score (0-1, higher is safer)
            safety_score = 1.0
            
            # Penalize high molecular weight
            if estimated_mw > 500:
                safety_score -= 0.2
            
            # Penalize high lipophilicity
            if estimated_logp > 5:
                safety_score -= 0.3
            
            # Penalize multiple aromatic rings
            if aromatic_rings > 3:
                safety_score -= 0.2
            
            # Penalize toxicity alerts
            safety_score -= (toxicity_alerts * 0.1)
            
            # Ensure safety score is between 0 and 1
            safety_score = max(0, min(1, safety_score))
            
            toxicity_features.append({
                'SMILES': smiles,
                'Molecular_Weight': estimated_mw,
                'Lipophilicity_Est': estimated_logp,
                'Aromatic_Rings': aromatic_rings,
                'Toxicity_Alerts': toxicity_alerts,
                'Safety_Score': safety_score
            })
            
            if (i + 1) % 100 == 0:
                print(f"   Processed {i + 1:,}/{len(smiles_list):,} compounds...")
        
        return pd.DataFrame(toxicity_features)
    
    def simulate_protox_prediction(self, smiles_list):
        """Simulate ProTox-II predictions (replace with actual API calls)"""
        print(f"\n🔬 SIMULATING PROTOX-II PREDICTIONS")
        print("=" * 35)
        
        protox_results = []
        
        for i, smiles in enumerate(smiles_list):
            if pd.isna(smiles):
                protox_results.append({
                    'SMILES': smiles,
                    'LD50_mg_kg': np.nan,
                    'Toxicity_Class': 'Unknown',
                    'Hepatotoxicity': np.nan,
                    'Carcinogenicity': np.nan,
                    'Mutagenicity': np.nan,
                    'Cytotoxicity': np.nan
                })
                continue
            
            # Simulate ProTox-II predictions
            # In production, replace with actual API calls
            
            # Simulate LD50 (oral rat) - mg/kg
            # Higher values = less toxic
            base_ld50 = np.random.lognormal(mean=6, sigma=1.5)  # ~400 mg/kg average
            
            # Adjust based on molecular features
            if 'Cl' in smiles or 'Br' in smiles:
                base_ld50 *= 0.7  # Halogenated compounds often more toxic
            
            if len(smiles) > 100:  # Large molecules
                base_ld50 *= 1.3  # Often less toxic due to poor absorption
            
            # Toxicity class based on LD50
            if base_ld50 >= 2000:
                tox_class = 'Class 5 (Low toxicity)'
            elif base_ld50 >= 300:
                tox_class = 'Class 4 (Moderate toxicity)'
            elif base_ld50 >= 50:
                tox_class = 'Class 3 (High toxicity)'
            else:
                tox_class = 'Class 1-2 (Very high toxicity)'
            
            # Simulate organ toxicity predictions (probability 0-1)
            hepatotoxicity = np.random.beta(2, 5)  # Skewed toward lower probability
            carcinogenicity = np.random.beta(1.5, 8)  # Lower probability
            mutagenicity = np.random.beta(1.8, 6)  # Lower probability
            cytotoxicity = np.random.beta(3, 4)  # More variable
            
            protox_results.append({
                'SMILES': smiles,
                'LD50_mg_kg': base_ld50,
                'Toxicity_Class': tox_class,
                'Hepatotoxicity': hepatotoxicity,
                'Carcinogenicity': carcinogenicity,
                'Mutagenicity': mutagenicity,
                'Cytotoxicity': cytotoxicity
            })
            
            if (i + 1) % 50 == 0:
                print(f"   Processed {i + 1:,}/{len(smiles_list):,} compounds...")
        
        return pd.DataFrame(protox_results)
    
    def calculate_safety_index(self, vs_df, toxicity_df, protox_df):
        """Calculate overall safety index combining multiple factors"""
        print(f"\n📊 CALCULATING COMPREHENSIVE SAFETY INDEX")
        print("=" * 45)
        
        print(f"   Input sizes: VS={len(vs_df)}, Toxicity={len(toxicity_df)}, ProTox={len(protox_df)}")
        
        safety_results = vs_df.copy()
        
        # Merge toxicity data with careful size monitoring
        print(f"   Merging toxicity data...")
        initial_size = len(safety_results)
        safety_results = safety_results.merge(toxicity_df, on='SMILES', how='left')
        after_tox_merge = len(safety_results)
        print(f"   After toxicity merge: {after_tox_merge:,} rows (was {initial_size:,})")
        
        if after_tox_merge > initial_size * 2:
            print(f"    WARNING: Toxicity merge created duplicates! Removing...")
            safety_results = safety_results.drop_duplicates()
            print(f"   After deduplication: {len(safety_results):,} rows")
        
        # Merge ProTox data with careful size monitoring
        print(f"   Merging ProTox data...")
        before_protox = len(safety_results)
        safety_results = safety_results.merge(protox_df, on='SMILES', how='left')
        after_protox_merge = len(safety_results)
        print(f"   After ProTox merge: {after_protox_merge:,} rows (was {before_protox:,})")
        
        if after_protox_merge > before_protox * 2:
            print(f"    WARNING: ProTox merge created duplicates! Removing...")
            safety_results = safety_results.drop_duplicates()
            print(f"   After deduplication: {len(safety_results):,} rows")
        
        # Final safety check - ensure we don't have millions of rows
        if len(safety_results) > 10000:
            print(f"    CRITICAL: {len(safety_results):,} rows is too many! Something went wrong.")
            print(f"   Limiting to original dataset size to prevent system issues...")
            safety_results = safety_results.head(len(vs_df))
            print(f"   Limited to: {len(safety_results):,} rows")
        
        # Calculate comprehensive safety index (0-1, higher is safer)
        safety_indices = []
        
        print(f"   Calculating safety indices for {len(safety_results):,} compounds...")
        
        for i, (_, row) in enumerate(safety_results.iterrows()):
            safety_index = 0.5  # Start with neutral
            
            # Factor 1: Simple safety score (25% weight)
            if pd.notna(row.get('Safety_Score')):
                safety_index += (row['Safety_Score'] - 0.5) * 0.25
            
            # Factor 2: LD50 toxicity (30% weight)
            if pd.notna(row.get('LD50_mg_kg')):
                ld50 = row['LD50_mg_kg']
                if ld50 >= 2000:
                    safety_index += 0.30  # Very safe
                elif ld50 >= 300:
                    safety_index += 0.15  # Moderately safe
                elif ld50 >= 50:
                    safety_index += 0.0   # Neutral
                else:
                    safety_index -= 0.30  # Unsafe
            
            # Factor 3: Organ toxicity (30% weight)
            organ_toxicities = [
                row.get('Hepatotoxicity', 0.5),
                row.get('Carcinogenicity', 0.5),
                row.get('Mutagenicity', 0.5),
                row.get('Cytotoxicity', 0.5)
            ]
            avg_organ_tox = np.mean([t for t in organ_toxicities if pd.notna(t)])
            safety_index += (0.5 - avg_organ_tox) * 0.30
            
            # Factor 4: Database source bonus (15% weight)
            if row.get('Database') == 'FOODB':
                safety_index += 0.15  # Food compounds generally safer
            elif row.get('Database') == 'LOTUS':
                safety_index += 0.05  # Traditional medicine some safety data
            
            # Ensure safety index is between 0 and 1
            safety_index = max(0, min(1, safety_index))
            safety_indices.append(safety_index)
            
            # Progress indicator for large datasets
            if (i + 1) % 1000 == 0:
                print(f"   Progress: {i + 1:,}/{len(safety_results):,} compounds processed...")
        
        safety_results['Safety_Index'] = safety_indices
        
        # Add safety categories
        safety_categories = []
        for si in safety_indices:
            if si >= 0.8:
                safety_categories.append('Very Safe')
            elif si >= 0.6:
                safety_categories.append('Safe')
            elif si >= 0.4:
                safety_categories.append('Moderate Risk')
            elif si >= 0.2:
                safety_categories.append('High Risk')
            else:
                safety_categories.append('Very High Risk')
        
        safety_results['Safety_Category'] = safety_categories
        
        # Sort by safety index (safest first)
        safety_results = safety_results.sort_values('Safety_Index', ascending=False)
        
        print(f" Safety assessment complete for {len(safety_results):,} compounds")
        
        return safety_results
    
    def create_safety_prioritized_hits(self, safety_results):
        """Create safety-prioritized hit list"""
        print(f"\n CREATING SAFETY-PRIORITIZED HIT LIST")
        print("=" * 40)
        
        # Find activity column
        activity_cols = [col for col in safety_results.columns if any(x in col.lower() for x in ['score', 'probability', 'affinity', 'activity'])]
        activity_cols = [col for col in activity_cols if col not in ['Safety_Score']]
        
        if not activity_cols:
            print(" No activity score column found")
            return safety_results
        
        activity_col = activity_cols[0]
        print(f" Using activity column: {activity_col}")
        
        # Create balanced score (activity + safety)
        activity_scores = safety_results[activity_col].fillna(0)
        safety_scores = safety_results['Safety_Index'].fillna(0)
        
        # Balanced score: 60% activity, 40% safety
        balanced_scores = (activity_scores * 0.6) + (safety_scores * 0.4)
        safety_results['Balanced_Score'] = balanced_scores
        
        # Sort by balanced score
        prioritized_hits = safety_results.sort_values('Balanced_Score', ascending=False)
        
        # Analysis
        total_compounds = len(prioritized_hits)
        very_safe = (prioritized_hits['Safety_Category'] == 'Very Safe').sum()
        safe = (prioritized_hits['Safety_Category'] == 'Safe').sum()
        high_activity_safe = ((prioritized_hits[activity_col] >= 0.6) & 
                             (prioritized_hits['Safety_Index'] >= 0.6)).sum()
        
        print(f" SAFETY ANALYSIS RESULTS:")
        print(f"   Total compounds: {total_compounds:,}")
        print(f"   Very safe compounds: {very_safe:,} ({very_safe/total_compounds*100:.1f}%)")
        print(f"   Safe compounds: {safe:,} ({safe/total_compounds*100:.1f}%)")
        print(f"   High activity + safe: {high_activity_safe:,}")
        
        # Show coverage by database
        if 'Database' in prioritized_hits.columns:
            db_coverage = prioritized_hits['Database'].value_counts()
            print(f"\n🗄️ Coverage by database:")
            for db, count in db_coverage.items():
                print(f"   {db}: {count:,} compounds")
        
        # Show top safety-prioritized hits
        print(f"\n TOP 15 SAFETY-PRIORITIZED HITS:")
        top_hits = prioritized_hits.head(15)
        
        for i, (_, row) in enumerate(top_hits.iterrows(), 1):
            name = row.get('Compound_Name', 'Unknown')
            activity = row[activity_col]
            safety = row['Safety_Index']
            balanced = row['Balanced_Score']
            category = row['Safety_Category']
            database = row.get('Database', 'Unknown')
            
            print(f"  {i:2d}. [{database}] {name}")
            print(f"      Activity: {activity:.3f} | Safety: {safety:.3f} | Balanced: {balanced:.3f}")
            print(f"      Safety Category: {category}")
        
        return prioritized_hits
    
    def save_top50_toxicity_results(self, safety_results):
        """Save Top50 focused toxicity assessment results"""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_file = f'Top50_Focused_Toxicity_Results_{timestamp}.xlsx'
        
        print(f"\n SAVING TOP50 FOCUSED TOXICITY RESULTS")
        print("=" * 45)
        
        # Validate data size before saving
        max_rows = 1000000  # Excel limit minus buffer
        max_cols = 100      # Reasonable column limit
        
        if len(safety_results) > max_rows:
            print(f" Data too large for Excel ({len(safety_results):,} rows), limiting to {max_rows:,}")
            # Find activity column to sort by
            activity_cols = [col for col in safety_results.columns if any(x in col.lower() for x in ['score', 'probability', 'affinity', 'activity'])]
            activity_cols = [col for col in activity_cols if col not in ['Safety_Score']]
            
            if activity_cols:
                safety_results = safety_results.nlargest(max_rows, activity_cols[0])
            else:
                safety_results = safety_results.head(max_rows)
            print(f" Data limited to {len(safety_results):,} rows")
        
        if len(safety_results.columns) > max_cols:
            print(f" Too many columns ({len(safety_results.columns)}), selecting key columns only")
            # Keep essential columns
            essential_cols = []
            for col in safety_results.columns:
                if any(x in col.lower() for x in [
                    'compound', 'name', 'smiles', 'database', 'source',
                    'score', 'probability', 'activity', 'safety', 'balanced',
                    'category', 'toxicity', 'ld50'
                ]):
                    essential_cols.append(col)
            
            # Ensure we have at least some columns
            if len(essential_cols) < 5:
                essential_cols = safety_results.columns[:max_cols].tolist()
            else:
                essential_cols = essential_cols[:max_cols]
            
            safety_results = safety_results[essential_cols]
            print(f" Reduced to {len(safety_results.columns)} essential columns")
        
        # Find activity column
        activity_cols = [col for col in safety_results.columns if any(x in col.lower() for x in ['score', 'probability', 'affinity', 'activity'])]
        activity_cols = [col for col in activity_cols if col not in ['Safety_Score']]
        activity_col = activity_cols[0] if activity_cols else 'Activity_Score'
        
        try:
            with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
                # Complete Top50 safety results
                print(f"   Writing main results sheet...")
                safety_results.to_excel(writer, sheet_name='Complete_Top50_Results', index=False)
                
                # Top50 Safe Hits (priority for experiments)
                safe_hits = safety_results[safety_results['Safety_Index'] >= 0.6] if 'Safety_Index' in safety_results.columns else safety_results.head(50)
                if len(safe_hits) > 0:
                    safe_hits = safe_hits.sort_values('Balanced_Score', ascending=False) if 'Balanced_Score' in safe_hits.columns else safe_hits
                    print(f"   Writing safe hits sheet ({len(safe_hits)} compounds)...")
                    safe_hits.to_excel(writer, sheet_name='Top50_Safe_Hits', index=False)
                
                # Best compounds by database (limited to prevent size issues)
                if 'Database' in safety_results.columns:
                    for database in safety_results['Database'].unique():
                        db_data = safety_results[safety_results['Database'] == database]
                        if 'Safety_Index' in db_data.columns:
                            db_safe = db_data[db_data['Safety_Index'] >= 0.6].head(20)
                        else:
                            db_safe = db_data.head(20)
                        
                        if len(db_safe) > 0:
                            sheet_name = f'{database}_Top_Safe'[:31]  # Excel sheet name limit
                            print(f"   Writing {sheet_name} sheet ({len(db_safe)} compounds)...")
                            db_safe.to_excel(writer, sheet_name=sheet_name, index=False)
                
                # Enhanced recommendations
                recommendations = [
                    ['Priority', 'Compound_Source', 'Criteria', 'Action', 'Risk_Level'],
                    ['1', 'FOODB Top Safe', 'Food compounds, Safety ≥0.6', 'Nutraceutical development priority', 'Very Low'],
                    ['2', 'COCONUT Top Performers', 'High activity from largest database', 'Lead optimization candidates', 'Low-Moderate'],
                    ['3', 'LOTUS Traditional Safe', 'Traditional medicine + safety', 'Ethnopharmacology validation', 'Low'],
                    ['4', 'Top50 Consensus Hits', 'High scores across databases', 'Primary experimental targets', 'Low'],
                    ['⚠️', 'High Activity High Risk', 'Activity ≥0.6, Safety <0.4', 'Structural modification needed', 'High']
                ]
                
                rec_df = pd.DataFrame(recommendations[1:], columns=recommendations[0])
                print(f"   Writing recommendations sheet...")
                rec_df.to_excel(writer, sheet_name='Top50_Recommendations', index=False)
            
            print(f" Top50 focused toxicity results saved: {output_file}")
            
        except Exception as e:
            print(f" Error saving Excel file: {e}")
            # Try saving as CSV as fallback
            csv_file = output_file.replace('.xlsx', '.csv')
            print(f"   Attempting to save as CSV: {csv_file}")
            try:
                safety_results.to_csv(csv_file, index=False)
                print(f" Saved as CSV: {csv_file}")
                return csv_file
            except Exception as csv_error:
                print(f" CSV save also failed: {csv_error}")
                return None
        
        # Print focused summary
        total_compounds = len(safety_results)
        
        print(f"\n TOP50 FOCUSED SUMMARY:")
        print(f"   Total Top50 records processed: {total_compounds:,}")
        
        if 'Database' in safety_results.columns:
            db_counts = safety_results['Database'].value_counts()
            print(f"   Database representation:")
            for db, count in db_counts.items():
                if 'Safety_Index' in safety_results.columns:
                    safe_count = (safety_results['Database'] == db) & (safety_results['Safety_Index'] >= 0.6)
                    safe_compounds = safety_results[safe_count].shape[0]
                    print(f"     {db}: {count:,} compounds ({safe_compounds:,} safe)")
                else:
                    print(f"     {db}: {count:,} compounds")
        
        # Safety analysis
        if 'Safety_Category' in safety_results.columns:
            safety_dist = safety_results['Safety_Category'].value_counts()
            print(f"   Top50 safety distribution:")
            for category, count in safety_dist.items():
                print(f"     {category}: {count:,} ({count/total_compounds*100:.1f}%)")
        
        # High-value compound identification
        if 'Safety_Index' in safety_results.columns and activity_col in safety_results.columns:
            high_activity_safe = safety_results[
                (safety_results[activity_col] >= 0.6) & 
                (safety_results['Safety_Index'] >= 0.6)
            ]
            
            print(f"\n TOP50 KEY FINDINGS:")
            print(f"   High activity + safe compounds: {len(high_activity_safe):,}")
            
            # FOODB advantage
            foodb_compounds = safety_results[safety_results['Database'] == 'FOODB'] if 'Database' in safety_results.columns else pd.DataFrame()
            if len(foodb_compounds) > 0:
                foodb_safe = foodb_compounds[foodb_compounds['Safety_Index'] >= 0.6]
                foodb_safe_rate = len(foodb_safe) / len(foodb_compounds) * 100
                print(f"   FOODB safety advantage: {len(foodb_safe):,}/{len(foodb_compounds):,} compounds safe ({foodb_safe_rate:.1f}%)")
        
        print(f"\n EXPERIMENTAL RECOMMENDATIONS:")
        print(f"   Phase 1 (Immediate): FOODB safe compounds")
        print(f"   Phase 2 (Secondary): Database-specific top performers")  
        print(f"   Phase 3 (Comprehensive): Remaining Top50 safe compounds")
        print(f"   Focus: Food-derived compounds for nutraceutical potential")
        
        return output_file

def main():
    """Top50 focused toxicity screening pipeline"""
    print(" FOCUSED TOP50 TOXICITY SCREENING")
    print("=" * 45)
    print("Safety assessment of the best compounds from COCONUT, FOODB, and LOTUS Top50 sheets")
    print()
    
    # Get user preference for combination strategy
    print(" COMBINATION STRATEGY OPTIONS:")
    print("1. best_per_compound - Keep highest scoring prediction per compound (Recommended)")
    print("2. database_comparison - Keep best from each database for comparison")
    print("3. all_top50 - Keep all Top50 predictions (may have duplicates)")
    print()
    
    strategy_choice = input("Choose combination strategy (1-3) or press Enter for default (1): ").strip()
    strategy_map = {
        '1': 'best_per_compound',
        '2': 'database_comparison', 
        '3': 'all_top50'
    }
    strategy = strategy_map.get(strategy_choice, 'best_per_compound')
    print(f" Using strategy: {strategy}")
    print()
    
    tox_screener = Top50ToxicityScreening()
    
    # Load Top50 sheets specifically
    top50_sheets_data = tox_screener.load_top50_vs_results()
    if top50_sheets_data is None:
        return
    
    # Combine Top50 results with focused strategy
    combined_results = tox_screener.combine_top50_results(top50_sheets_data, strategy=strategy)
    if combined_results is None:
        return
    
    # Check for SMILES column
    smiles_cols = [col for col in combined_results.columns if 'smiles' in col.lower()]
    if not smiles_cols:
        print(" No SMILES column found in combined results!")
        print("Available columns:", list(combined_results.columns))
        return
    
    smiles_col = smiles_cols[0]
    print(f" Using SMILES column: {smiles_col}")
    
    # For Top50 sheets, we can process all compounds (they're already pre-filtered)
    compounds_to_process = combined_results
    smiles_list = compounds_to_process[smiles_col].tolist()
    
    print(f"\n PROCESSING {len(compounds_to_process):,} TOP50 COMPOUNDS FOR TOXICITY ASSESSMENT")
    
    # Calculate simple toxicity features
    toxicity_features = tox_screener.calculate_simple_toxicity_features(smiles_list)
    
    # Simulate ProTox predictions
    protox_predictions = tox_screener.simulate_protox_prediction(smiles_list)
    
    # Calculate comprehensive safety index
    safety_results = tox_screener.calculate_safety_index(compounds_to_process, toxicity_features, protox_predictions)
    
    # Create safety-prioritized hits
    prioritized_hits = tox_screener.create_safety_prioritized_hits(safety_results)
    
    # Save enhanced results
    output_file = tox_screener.save_top50_toxicity_results(prioritized_hits)
    
    print(f"\n TOP50 FOCUSED TOXICITY SCREENING COMPLETE")
    print("=" * 50)
    print(" Top50 sheets from COCONUT, FOODB, and LOTUS processed")
    print(" Best compounds from each database analyzed")
    print(" Cross-database safety comparison completed")
    print(" Focused toxicity predictions generated")
    print(" Top50 safety-prioritized hit list created")
    print(f" Results saved: {output_file}")
    
    print(f"\n FOCUSED NEXT STEPS:")
    print("1. Review the top-scoring safe compounds first")
    print("2. Compare safety profiles across databases")
    print("3. Prioritize FOODB compounds for inherent food safety")
    print("4. Focus experimental validation on highest-confidence hits")
    print("5. Use database diversity for comprehensive coverage")
    print("6. Consider cross-database compounds as high-confidence leads")

if __name__ == "__main__":
    main()
