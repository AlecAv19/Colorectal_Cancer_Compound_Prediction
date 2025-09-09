# -*- coding: utf-8 -*-
"""
Created on Tue Aug 19 09:24:19 2025

@author: ASUS TUF F15
"""

"""
GDSC Z-SCORE MULTI-CELL LINE ANALYZER

This script analyzes GDSC data with Z-score based activity classification:
- Z-score < -2.0: Sensitive (Active = 1)
- Z-score > 2.0: Resistant (Active = 0) 
- -2.0 ≤ Z-score ≤ 2.0: Inconclusive (excluded)

Expected columns: ID, Drug Name, Targets, Z Score, Count

Usage: python gdsc_zscore_analyzer.py
"""

import pandas as pd
import numpy as np
import os
import glob
from collections import defaultdict, Counter
import warnings
warnings.filterwarnings('ignore')

class GDSCZScoreAnalyzer:
    """
    Analyze GDSC Z-score data across multiple cell lines
    """
    
    def __init__(self):
        self.cell_lines = ['SW620', 'RKO', 'LoVo', 'HCT-15', 'Colo-205']
        self.raw_data = {}
        self.processed_data = {}
        self.consistent_compounds = {}
        self.unified_dataset = None
        
        # Z-score thresholds
        self.sensitive_threshold = -2.0  # Z < -2.0 = sensitive (active)
        self.resistant_threshold = 2.0   # Z > 2.0 = resistant (inactive)
        
    def load_gdsc_zscore_files(self, data_directory='.'):
        """
        Load GDSC Z-score data files for each cell line
        """
        print("=== STEP 1: LOADING GDSC Z-SCORE DATA ===")
        print(f"Looking for data files in: {data_directory}")
        print(f"Target cell lines: {', '.join(self.cell_lines)}")
        print(f"Expected columns: ID, Drug Name, Targets, Z Score, Count")
        
        # File patterns for GDSC data
        file_patterns = [
            '*{cell_line}*.xlsx', '*{cell_line}*.csv',
            '*{cell_line}*.xls', 'GDSC*{cell_line}*',
            '{cell_line}*.xlsx', '{cell_line}*.csv',
            '*{cell_line}_*.xlsx', '*{cell_line}_*.csv'
        ]
        
        found_files = {}
        
        for cell_line in self.cell_lines:
            print(f"\n🔍 Searching for {cell_line} data...")
            
            cell_files = []
            for pattern in file_patterns:
                # Try exact case
                search_pattern = os.path.join(data_directory, pattern.format(cell_line=cell_line))
                matches = glob.glob(search_pattern)
                cell_files.extend(matches)
                
                # Try lowercase
                search_pattern_lower = os.path.join(data_directory, pattern.format(cell_line=cell_line.lower()))
                matches_lower = glob.glob(search_pattern_lower)
                cell_files.extend(matches_lower)
                
                # Try variations (e.g., HCT-15 vs HCT15)
                if '-' in cell_line:
                    cell_line_alt = cell_line.replace('-', '')
                    search_pattern_alt = os.path.join(data_directory, pattern.format(cell_line=cell_line_alt))
                    matches_alt = glob.glob(search_pattern_alt)
                    cell_files.extend(matches_alt)
            
            # Remove duplicates
            cell_files = list(set(cell_files))
            
            if cell_files:
                print(f"  Found {len(cell_files)} potential files:")
                for f in cell_files:
                    print(f"    📄 {os.path.basename(f)}")
                
                # Select the first file
                selected_file = cell_files[0]
                found_files[cell_line] = selected_file
                print(f"  ✅ Selected: {os.path.basename(selected_file)}")
            else:
                print(f"  ❌ No files found for {cell_line}")
        
        if not found_files:
            print("\n❌ No GDSC data files found!")
            print("Please ensure your GDSC files are in the current directory")
            print("Expected naming examples:")
            print("  • SW620.xlsx, RKO.csv, LoVo_data.xlsx")
            print("  • GDSC_SW620.xlsx, HCT-15_results.csv")
            print("  • Colo-205.xlsx")
            return False
        
        # Load and validate the files
        print(f"\n📂 Loading and validating {len(found_files)} files...")
        
        for cell_line, file_path in found_files.items():
            try:
                print(f"\n🔄 Processing {cell_line}: {os.path.basename(file_path)}")
                
                # Load file
                if file_path.endswith(('.xlsx', '.xls')):
                    # Handle Excel files
                    excel_file = pd.ExcelFile(file_path)
                    sheet_names = excel_file.sheet_names
                    print(f"  📋 Available sheets: {sheet_names}")
                    
                    # Try to find the right sheet
                    data_sheet = None
                    for sheet in sheet_names:
                        df_test = pd.read_excel(file_path, sheet_name=sheet, nrows=5)
                        if self._validate_gdsc_columns(df_test):
                            data_sheet = sheet
                            break
                    
                    if data_sheet:
                        df = pd.read_excel(file_path, sheet_name=data_sheet)
                        print(f"  ✅ Using sheet: {data_sheet}")
                    else:
                        df = pd.read_excel(file_path)
                        print(f"  ⚠️  Using default sheet")
                        
                elif file_path.endswith('.csv'):
                    df = pd.read_csv(file_path)
                    
                else:
                    print(f"  ❌ Unsupported file format")
                    continue
                
                print(f"  📊 Loaded: {len(df)} rows, {len(df.columns)} columns")
                print(f"  📝 Columns: {list(df.columns)}")
                
                # Validate and clean column names
                df_cleaned = self._clean_and_validate_data(df, cell_line)
                
                if df_cleaned is not None:
                    self.raw_data[cell_line] = {
                        'data': df_cleaned,
                        'file_path': file_path,
                        'original_rows': len(df),
                        'processed_rows': len(df_cleaned)
                    }
                    print(f"  ✅ Successfully processed: {len(df_cleaned)} valid rows")
                else:
                    print(f"  ❌ Failed validation for {cell_line}")
                    
            except Exception as e:
                print(f"  ❌ Error loading {cell_line}: {e}")
                continue
        
        if self.raw_data:
            print(f"\n✅ Successfully loaded {len(self.raw_data)} cell line datasets:")
            total_compounds = 0
            for cell_line, info in self.raw_data.items():
                compounds = len(info['data'])
                total_compounds += compounds
                print(f"  {cell_line}: {compounds} compounds")
            print(f"  Total: {total_compounds} compound-cell line combinations")
            return True
        else:
            print("\n❌ No valid data files loaded")
            return False
    
    def _validate_gdsc_columns(self, df):
        """
        Check if dataframe has expected GDSC columns
        """
        expected_patterns = {
            'id': ['id', 'drug_id', 'compound_id'],
            'name': ['drug name', 'name', 'compound_name', 'drug_name'],
            'targets': ['targets', 'target', 'pathway'],
            'zscore': ['z score', 'zscore', 'z_score', 'score'],
            'count': ['count', 'n', 'num', 'frequency']
        }
        
        columns_lower = [col.lower().strip() for col in df.columns]
        
        found_patterns = 0
        for pattern_type, patterns in expected_patterns.items():
            for pattern in patterns:
                if any(pattern in col for col in columns_lower):
                    found_patterns += 1
                    break
        
        return found_patterns >= 3  # Need at least 3 key columns
    
    def _clean_and_validate_data(self, df, cell_line):
        """
        Clean and standardize GDSC data columns
        """
        print(f"    🧹 Cleaning data for {cell_line}...")
        
        # Create column mapping
        column_mapping = {}
        columns_lower = {col.lower().strip(): col for col in df.columns}
        
        # Map ID column
        id_patterns = ['id', 'drug_id', 'compound_id', 'drugid']
        for pattern in id_patterns:
            if pattern in columns_lower:
                column_mapping['ID'] = columns_lower[pattern]
                break
        
        # Map Drug Name column
        name_patterns = ['drug name', 'name', 'compound_name', 'drug_name', 'drugname']
        for pattern in name_patterns:
            if pattern in columns_lower:
                column_mapping['Drug_Name'] = columns_lower[pattern]
                break
        
        # Map Targets column
        target_patterns = ['targets', 'target', 'pathway', 'mechanism']
        for pattern in target_patterns:
            if pattern in columns_lower:
                column_mapping['Targets'] = columns_lower[pattern]
                break
        
        # Map Z Score column
        zscore_patterns = ['z score', 'zscore', 'z_score', 'score', 'z-score']
        for pattern in zscore_patterns:
            if pattern in columns_lower:
                column_mapping['Z_Score'] = columns_lower[pattern]
                break
        
        # Map Count column
        count_patterns = ['count', 'n', 'num', 'frequency', 'observations']
        for pattern in count_patterns:
            if pattern in columns_lower:
                column_mapping['Count'] = columns_lower[pattern]
                break
        
        print(f"    📋 Column mapping: {column_mapping}")
        
        # Check if we have the essential columns
        if 'Z_Score' not in column_mapping:
            print(f"    ❌ No Z-Score column found")
            return None
        
        if 'ID' not in column_mapping and 'Drug_Name' not in column_mapping:
            print(f"    ❌ No compound identifier found")
            return None
        
        # Create cleaned dataframe
        cleaned_df = pd.DataFrame()
        cleaned_df['cell_line'] = cell_line
        
        # Map columns
        for standard_name, original_name in column_mapping.items():
            if original_name in df.columns:
                cleaned_df[standard_name] = df[original_name]
        
        # Add missing columns with defaults
        for col in ['ID', 'Drug_Name', 'Targets', 'Z_Score', 'Count']:
            if col not in cleaned_df.columns:
                cleaned_df[col] = None
        
        # Clean and validate Z-scores
        if cleaned_df['Z_Score'].notna().sum() > 0:
            # Convert to numeric, handling any string formats
            z_scores = pd.to_numeric(cleaned_df['Z_Score'], errors='coerce')
            cleaned_df['Z_Score'] = z_scores
            
            # Remove rows with invalid Z-scores
            valid_zscores = cleaned_df['Z_Score'].notna()
            cleaned_df = cleaned_df[valid_zscores].copy()
            
            print(f"    📈 Valid Z-scores: {len(cleaned_df)}")
            print(f"    📊 Z-score range: {cleaned_df['Z_Score'].min():.2f} to {cleaned_df['Z_Score'].max():.2f}")
        else:
            print(f"    ❌ No valid Z-scores found")
            return None
        
        # Create compound identifier (prefer ID, fallback to Drug_Name)
        if cleaned_df['ID'].notna().sum() > 0:
            cleaned_df['compound_identifier'] = cleaned_df['ID'].astype(str)
            print(f"    🆔 Using ID as identifier: {cleaned_df['ID'].notna().sum()} compounds")
        elif cleaned_df['Drug_Name'].notna().sum() > 0:
            cleaned_df['compound_identifier'] = cleaned_df['Drug_Name'].astype(str)
            print(f"    🏷️  Using Drug_Name as identifier: {cleaned_df['Drug_Name'].notna().sum()} compounds")
        else:
            print(f"    ❌ No valid compound identifiers")
            return None
        
        # Filter out compounds without identifiers
        valid_compounds = cleaned_df['compound_identifier'].notna() & (cleaned_df['compound_identifier'] != 'nan')
        cleaned_df = cleaned_df[valid_compounds].copy()
        
        print(f"    ✅ Final dataset: {len(cleaned_df)} compounds with valid data")
        
        return cleaned_df if len(cleaned_df) > 0 else None
    
    def classify_activity_by_zscore(self):
        """
        Classify compounds as active/inactive based on Z-scores
        """
        print(f"\n=== STEP 2: CLASSIFYING ACTIVITY BY Z-SCORES ===")
        print(f"Activity classification rules:")
        print(f"  • Sensitive (Active): Z-score < {self.sensitive_threshold}")
        print(f"  • Resistant (Inactive): Z-score > {self.resistant_threshold}")
        print(f"  • Inconclusive: {self.sensitive_threshold} ≤ Z-score ≤ {self.resistant_threshold} (excluded)")
        
        for cell_line, data_info in self.raw_data.items():
            print(f"\n🔄 Processing {cell_line}...")
            df = data_info['data'].copy()
            
            # Classify based on Z-score
            conditions = [
                df['Z_Score'] < self.sensitive_threshold,  # Sensitive
                df['Z_Score'] > self.resistant_threshold   # Resistant
            ]
            choices = [1, 0]  # 1 = active, 0 = inactive
            
            df['activity'] = np.select(conditions, choices, default=np.nan)
            
            # Count classifications
            sensitive_count = (df['activity'] == 1).sum()
            resistant_count = (df['activity'] == 0).sum()
            inconclusive_count = df['activity'].isna().sum()
            
            print(f"  📊 Classification results:")
            print(f"    Sensitive (active): {sensitive_count}")
            print(f"    Resistant (inactive): {resistant_count}")
            print(f"    Inconclusive (excluded): {inconclusive_count}")
            
            # Keep only conclusive classifications
            conclusive_data = df[df['activity'].notna()].copy()
            
            print(f"  ✅ Kept {len(conclusive_data)} compounds with conclusive activity")
            
            self.processed_data[cell_line] = conclusive_data
        
        return True
    
    def find_multi_cellline_compounds(self, min_cell_lines=3):
        """
        Find compounds tested in multiple cell lines with consistent activity
        """
        print(f"\n=== STEP 3: FINDING MULTI-CELL LINE COMPOUNDS (≥{min_cell_lines}) ===")
        
        # Collect all compound data
        compound_data = defaultdict(list)
        
        for cell_line, df in self.processed_data.items():
            print(f"\n📊 Analyzing {cell_line}: {len(df)} compounds")
            
            for idx, row in df.iterrows():
                compound_id = row['compound_identifier']
                compound_data[compound_id].append({
                    'cell_line': cell_line,
                    'z_score': row['Z_Score'],
                    'activity': row['activity'],
                    'drug_name': row.get('Drug_Name', ''),
                    'targets': row.get('Targets', ''),
                    'count': row.get('Count', 0)
                })
        
        print(f"\n🔍 Compound consistency analysis:")
        print(f"  Total unique compounds: {len(compound_data)}")
        
        # Filter by minimum cell line requirement
        multi_cellline_compounds = {}
        
        for compound_id, measurements in compound_data.items():
            cell_lines_tested = set(m['cell_line'] for m in measurements)
            
            if len(cell_lines_tested) >= min_cell_lines:
                # Analyze activity consistency
                activities = [m['activity'] for m in measurements]
                z_scores = [m['z_score'] for m in measurements]
                
                # Calculate consistency metrics
                activity_counts = Counter(activities)
                most_common_activity = activity_counts.most_common(1)[0]
                consistency_ratio = most_common_activity[1] / len(activities)
                
                multi_cellline_compounds[compound_id] = {
                    'measurements': measurements,
                    'cell_lines': sorted(list(cell_lines_tested)),
                    'num_cell_lines': len(cell_lines_tested),
                    'activities': activities,
                    'z_scores': z_scores,
                    'consensus_activity': most_common_activity[0],
                    'consistency_ratio': consistency_ratio,
                    'mean_z_score': np.mean(z_scores),
                    'std_z_score': np.std(z_scores),
                    'drug_name': measurements[0]['drug_name'],
                    'targets': measurements[0]['targets']
                }
        
        print(f"  Compounds in ≥{min_cell_lines} cell lines: {len(multi_cellline_compounds)}")
        
        # Show distribution by cell line count
        cellline_distribution = defaultdict(int)
        for compound_info in multi_cellline_compounds.values():
            cellline_distribution[compound_info['num_cell_lines']] += 1
        
        print(f"\n📈 Distribution by number of cell lines:")
        for num_lines in sorted(cellline_distribution.keys(), reverse=True):
            count = cellline_distribution[num_lines]
            print(f"  {num_lines} cell lines: {count} compounds")
        
        # Show consistency analysis
        high_consistency = sum(1 for info in multi_cellline_compounds.values() 
                              if info['consistency_ratio'] >= 0.8)
        medium_consistency = sum(1 for info in multi_cellline_compounds.values() 
                                if 0.6 <= info['consistency_ratio'] < 0.8)
        low_consistency = sum(1 for info in multi_cellline_compounds.values() 
                             if info['consistency_ratio'] < 0.6)
        
        print(f"\n🎯 Activity consistency analysis:")
        print(f"  High consistency (≥80%): {high_consistency} compounds")
        print(f"  Medium consistency (60-80%): {medium_consistency} compounds")
        print(f"  Low consistency (<60%): {low_consistency} compounds")
        
        self.consistent_compounds = multi_cellline_compounds
        return len(multi_cellline_compounds) > 0
    
    def create_unified_zscore_dataset(self, min_consistency=0.6):
        """
        Create unified dataset with quality filtering
        """
        print(f"\n=== STEP 4: CREATING UNIFIED DATASET ===")
        print(f"Minimum consistency threshold: {min_consistency}")
        
        if not self.consistent_compounds:
            print("❌ No consistent compounds available")
            return False
        
        unified_records = []
        
        for compound_id, compound_info in self.consistent_compounds.items():
            # Apply quality filters
            if compound_info['consistency_ratio'] < min_consistency:
                continue  # Skip inconsistent compounds
            
            # Create unified record
            record = {
                'compound_id': compound_id,
                'drug_name': compound_info['drug_name'],
                'targets': compound_info['targets'],
                'activity': compound_info['consensus_activity'],
                'consistency_ratio': compound_info['consistency_ratio'],
                'num_cell_lines': compound_info['num_cell_lines'],
                'cell_lines_tested': ','.join(compound_info['cell_lines']),
                'mean_z_score': compound_info['mean_z_score'],
                'std_z_score': compound_info['std_z_score'],
                'min_z_score': min(compound_info['z_scores']),
                'max_z_score': max(compound_info['z_scores']),
                'z_score_range': max(compound_info['z_scores']) - min(compound_info['z_scores']),
                'num_measurements': len(compound_info['measurements'])
            }
            
            # Add individual cell line Z-scores
            for cell_line in self.cell_lines:
                cell_measurements = [m for m in compound_info['measurements'] 
                                   if m['cell_line'] == cell_line]
                if cell_measurements:
                    record[f'{cell_line}_Z_Score'] = cell_measurements[0]['z_score']
                    record[f'{cell_line}_Activity'] = cell_measurements[0]['activity']
                else:
                    record[f'{cell_line}_Z_Score'] = None
                    record[f'{cell_line}_Activity'] = None
            
            unified_records.append(record)
        
        if unified_records:
            self.unified_dataset = pd.DataFrame(unified_records)
            
            # Calculate final statistics
            total_compounds = len(self.unified_dataset)
            active_compounds = (self.unified_dataset['activity'] == 1).sum()
            inactive_compounds = (self.unified_dataset['activity'] == 0).sum()
            
            print(f"✅ Unified dataset created:")
            print(f"  Total compounds: {total_compounds}")
            print(f"  Active (sensitive): {active_compounds}")
            print(f"  Inactive (resistant): {inactive_compounds}")
            print(f"  Activity ratio (A:I): 1:{inactive_compounds/active_compounds:.1f}" if active_compounds > 0 else "")
            
            # Quality metrics
            high_quality = (self.unified_dataset['consistency_ratio'] >= 0.8).sum()
            tested_in_4plus = (self.unified_dataset['num_cell_lines'] >= 4).sum()
            tested_in_all = (self.unified_dataset['num_cell_lines'] == 5).sum()
            
            print(f"\n📊 Quality metrics:")
            print(f"  High consistency (≥80%): {high_quality}")
            print(f"  Tested in ≥4 cell lines: {tested_in_4plus}")
            print(f"  Tested in all 5 cell lines: {tested_in_all}")
            
            return True
        else:
            print("❌ No records passed quality filters")
            return False
    
    def save_zscore_results(self, output_file="GDSC_ZScore_Unified_Dataset.xlsx", min_cell_lines=3):
        """
        Save comprehensive Z-score analysis results
        """
        print(f"\n=== STEP 5: SAVING RESULTS ===")
        
        try:
            with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
                
                # Main unified dataset
                if self.unified_dataset is not None:
                    self.unified_dataset.to_excel(writer, sheet_name='Unified_Dataset', index=False)
                    
                    # Separate by activity
                    active_compounds = self.unified_dataset[self.unified_dataset['activity'] == 1]
                    inactive_compounds = self.unified_dataset[self.unified_dataset['activity'] == 0]
                    
                    if not active_compounds.empty:
                        active_compounds.to_excel(writer, sheet_name='Sensitive_Compounds', index=False)
                    if not inactive_compounds.empty:
                        inactive_compounds.to_excel(writer, sheet_name='Resistant_Compounds', index=False)
                    
                    # High quality subsets
                    high_consistency = self.unified_dataset[self.unified_dataset['consistency_ratio'] >= 0.8]
                    tested_all_lines = self.unified_dataset[self.unified_dataset['num_cell_lines'] == 5]
                    
                    if not high_consistency.empty:
                        high_consistency.to_excel(writer, sheet_name='High_Consistency', index=False)
                    if not tested_all_lines.empty:
                        tested_all_lines.to_excel(writer, sheet_name='All_CellLines', index=False)
                
                # Individual cell line processed data
                for cell_line, df in self.processed_data.items():
                    sheet_name = f'{cell_line}_Processed'
                    df.to_excel(writer, sheet_name=sheet_name, index=False)
                
                # Raw data for reference
                for cell_line, data_info in self.raw_data.items():
                    df = data_info['data']
                    sheet_name = f'{cell_line}_Raw'
                    df.to_excel(writer, sheet_name=sheet_name, index=False)
                
                # Analysis summary
                summary_data = []
                for cell_line in self.cell_lines:
                    if cell_line in self.raw_data:
                        raw_info = self.raw_data[cell_line]
                        processed_df = self.processed_data.get(cell_line, pd.DataFrame())
                        
                        sensitive_count = (processed_df['activity'] == 1).sum() if not processed_df.empty else 0
                        resistant_count = (processed_df['activity'] == 0).sum() if not processed_df.empty else 0
                        
                        summary_data.append({
                            'Cell_Line': cell_line,
                            'File_Path': os.path.basename(raw_info['file_path']),
                            'Raw_Compounds': raw_info['original_rows'],
                            'Valid_ZScores': raw_info['processed_rows'],
                            'Sensitive_Compounds': sensitive_count,
                            'Resistant_Compounds': resistant_count,
                            'Conclusive_Compounds': sensitive_count + resistant_count
                        })
                    else:
                        summary_data.append({
                            'Cell_Line': cell_line,
                            'File_Path': 'NOT FOUND',
                            'Raw_Compounds': 0,
                            'Valid_ZScores': 0,
                            'Sensitive_Compounds': 0,
                            'Resistant_Compounds': 0,
                            'Conclusive_Compounds': 0
                        })
                
                summary_df = pd.DataFrame(summary_data)
                summary_df.to_excel(writer, sheet_name='CellLine_Summary', index=False)
                
                # Processing parameters and results
                if self.unified_dataset is not None:
                    results_summary = pd.DataFrame({
                        'Parameter': [
                            'Sensitive_Threshold', 'Resistant_Threshold',
                            'Min_Cell_Lines_Required', 'Min_Consistency_Required',
                            'Total_Cell_Lines_Analyzed', 'Total_Unified_Compounds',
                            'Sensitive_Compounds', 'Resistant_Compounds',
                            'High_Consistency_Compounds', 'All_CellLines_Compounds',
                            'Processing_Timestamp'
                        ],
                        'Value': [
                            self.sensitive_threshold, self.resistant_threshold,
                            min_cell_lines, 0.6,  # min_consistency used in create_unified_zscore_dataset
                            len(self.raw_data), len(self.unified_dataset),
                            (self.unified_dataset['activity'] == 1).sum(),
                            (self.unified_dataset['activity'] == 0).sum(),
                            (self.unified_dataset['consistency_ratio'] >= 0.8).sum(),
                            (self.unified_dataset['num_cell_lines'] == 5).sum(),
                            pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')
                        ]
                    })
                    results_summary.to_excel(writer, sheet_name='Processing_Results', index=False)
            
            print(f"✅ Results saved to: {output_file}")
            
            print(f"\n🎉 GDSC Z-SCORE ANALYSIS COMPLETED!")
            print(f"📁 Output file: {output_file}")
            
            if self.unified_dataset is not None:
                print(f"\n🎯 FINAL SUMMARY:")
                print(f"  • Cell lines analyzed: {len(self.raw_data)}")
                print(f"  • Unified compounds: {len(self.unified_dataset)}")
                print(f"  • Sensitive (active): {(self.unified_dataset['activity'] == 1).sum()}")
                print(f"  • Resistant (inactive): {(self.unified_dataset['activity'] == 0).sum()}")
                print(f"  • High consistency: {(self.unified_dataset['consistency_ratio'] >= 0.8).sum()}")
                print(f"  • Tested in all 5 cell lines: {(self.unified_dataset['num_cell_lines'] == 5).sum()}")
            
            print(f"\n🚀 NEXT STEPS:")
            print(f"  1. Review unified dataset in '{output_file}'")
            print(f"  2. Use for molecular property calculation and SMILES lookup")
            print(f"  3. Generate property-matched decoys")
            print(f"  4. Train machine learning models with robust activity labels")
            
            return output_file
            
        except Exception as e:
            print(f"❌ Error saving results: {e}")
            import traceback
            traceback.print_exc()
            return None

def run_gdsc_zscore_analysis(min_cell_lines=2, data_directory='.'):
    """
    Run complete GDSC Z-score analysis workflow
    """
    print("🧬 GDSC Z-SCORE MULTI-CELL LINE ANALYZER")
    print("Activity Classification Based on Z-Scores")
    print("=" * 65)
    print(f"📊 Z-Score Thresholds:")
    print(f"   • Sensitive (Active): Z < -2.0")
    print(f"   • Resistant (Inactive): Z > 2.0")
    print(f"   • Inconclusive: -2.0 ≤ Z ≤ 2.0 (excluded)")
    print(f"⚙️  Min cell lines required: {min_cell_lines}")
    print(f"📂 Data directory: {data_directory}")
    print(f"🎯 Target cell lines: SW620, RKO, LoVo, HCT-15, Colo-205")
    print("=" * 65)
    
    analyzer = GDSCZScoreAnalyzer()
    
    # Step 1: Load GDSC Z-score files
    if not analyzer.load_gdsc_zscore_files(data_directory):
        print("\n❌ Failed to load GDSC files")
        print("Please ensure files contain columns: ID, Drug Name, Targets, Z Score, Count")
        return None
    
    # Step 2: Classify activity by Z-scores
    if not analyzer.classify_activity_by_zscore():
        print("\n❌ Failed to classify activities")
        return None
    
    # Step 3: Find multi-cell line compounds
    if not analyzer.find_multi_cellline_compounds(min_cell_lines):
        print(f"\n❌ No compounds found in ≥{min_cell_lines} cell lines")
        print("Try reducing min_cell_lines parameter")
        return None
    
    # Step 4: Create unified dataset
    if not analyzer.create_unified_zscore_dataset():
        print("\n❌ Failed to create unified dataset")
        return None
    
    # Step 5: Save results
    output_file = analyzer.save_zscore_results(
        output_file="GDSC_ZScore_Unified_Dataset.xlsx",
        min_cell_lines=min_cell_lines
    )
    
    return output_file

if __name__ == "__main__":
    # Configuration
    MIN_CELL_LINES = 2  # Require compounds tested in at least 3 cell lines
    DATA_DIRECTORY = '.'  # Current directory
    
    print("🚀 Starting GDSC Z-Score Analysis...")
    
    try:
        output_file = run_gdsc_zscore_analysis(
            min_cell_lines=MIN_CELL_LINES,
            data_directory=DATA_DIRECTORY
        )
        
        if output_file:
            print(f"\n🎉 SUCCESS! GDSC Z-score analysis completed!")
            print(f"📁 Results: {output_file}")
            print(f"\n💡 Your data is now ready for:")
            print(f"   • SMILES lookup (ChEMBL API or database)")
            print(f"   • Molecular property calculation")
            print(f"   • Decoy generation")
            print(f"   • Machine learning model training")
            print(f"\n🔬 Data Quality Achieved:")
            print(f"   • Consensus activity labels from multiple cell lines")
            print(f"   • Consistent Z-score based classification")
            print(f"   • Quality metrics for each compound")
            print(f"   • Ready for robust ML training")
        else:
            print(f"\n❌ Analysis failed - check error messages above")
            print(f"\n🔧 Troubleshooting tips:")
            print(f"   • Ensure files are named with cell line identifiers")
            print(f"   • Check column names match GDSC format")
            print(f"   • Verify Z-score values are numeric")
            print(f"   • Try reducing min_cell_lines to 2")
    
    except Exception as e:
        print(f"\n❌ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        print(f"\n📧 If the error persists, please check:")
        print(f"   • File formats (Excel/CSV)")
        print(f"   • Column names and data types") 
        print(f"   • File permissions and paths")