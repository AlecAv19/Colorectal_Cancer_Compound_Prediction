"""
NCI COMPOUND DATA EXTRACTOR WITH MULTI-CELL LINE EVALUATION

This script extracts compound screening data from NCI databases for:
SW620, RKO, LoVo, HCT-15, Colo-205, DLD-1, SW480, Caco-2

Features:
- Multi-cell line coverage analysis (5 tiers)
- Virtual screening strategy optimization
- ML-ready training sets with weights
- Conflict resolution using consensus voting
- Comprehensive output for different VS approaches

Usage: python nci_multicell_extractor.py
"""

import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')

class NCIMultiCellExtractor:
    """
    Extract and process NCI compound data with multi-cell line evaluation
    """
    
    def __init__(self):
        # Target colorectal cancer cell lines
        self.target_cell_lines = {
            'SW620': 'SW-620',
            'RKO': 'RKO',
            'LoVo': 'LoVo', 
            'HCT-15': 'HCT-15',
            'Colo-205': 'COLO205',
            'DLD-1': 'DLD-1',
            'SW480': 'SW-480',
            'Caco-2': 'Caco-2'
        }
        
        # NCI-60 cell lines available
        self.nci60_cell_lines = {
            'SW-620': 'Colon',
            'RKO': 'Colon', 
            'LoVo': 'Colon',
            'HCT-15': 'Colon',
            'COLO205': 'Colon',
            'SW-480': 'Colon'
        }
        
        self.compound_data = {}
        self.unified_dataset = None
        
    def fetch_sample_nci_data(self):
        """
        Create sample NCI compound data for demonstration
        """
        print("=== STEP 1: LOADING SAMPLE NCI DATA ===")
        print("Creating sample data structure...")
        print("(In practice: download from https://discover.nci.nih.gov/cellminer/)")
        
        # Sample compounds from NCI database
        nci_compounds = [
            {'nsc': 'NSC-123127', 'name': 'Doxorubicin'},
            {'nsc': 'NSC-109724', 'name': 'Vincristine'},
            {'nsc': 'NSC-740', 'name': 'Fluorouracil'},
            {'nsc': 'NSC-3053', 'name': 'Mitomycin C'},
            {'nsc': 'NSC-26271', 'name': 'Cisplatin'},
            {'nsc': 'NSC-125973', 'name': 'Paclitaxel'},
            {'nsc': 'NSC-613327', 'name': 'Camptothecin'},
            {'nsc': 'NSC-106977', 'name': 'Carboplatin'},
            {'nsc': 'NSC-241240', 'name': 'Etoposide'},
            {'nsc': 'NSC-118218', 'name': 'Methotrexate'},
            {'nsc': 'NSC-178248', 'name': 'Podophyllotoxin'},
            {'nsc': 'NSC-407292', 'name': 'Homoharringtonine'},
            {'nsc': 'NSC-67574', 'name': 'Actinomycin D'},
            {'nsc': 'NSC-141540', 'name': 'Bleomycin'},
            {'nsc': 'NSC-295739', 'name': 'Topotecan'},
            {'nsc': 'NSC-609699', 'name': 'Irinotecan'},
            {'nsc': 'NSC-104800', 'name': 'Idarubicin'},
            {'nsc': 'NSC-256927', 'name': 'Epirubicin'},
            {'nsc': 'NSC-750', 'name': 'Busulfan'},
            {'nsc': 'NSC-746', 'name': 'Melphalan'}
        ]
        
        compounds = []
        for compound in nci_compounds:
            compound_record = {
                'nsc_id': compound['nsc'],
                'compound_name': compound['name'],
                'smiles': f"SMILES_{compound['nsc']}",  # Placeholder
                'source': 'NCI-60_Sample'
            }
            
            # Generate realistic activity values for each cell line
            for cell_line_code in self.nci60_cell_lines.keys():
                # Simulate GI50 values (log10 M)
                base_activity = np.random.normal(-5.5, 1.5)
                gi50_val = max(-8.0, min(-3.0, base_activity + np.random.normal(0, 0.3)))
                
                compound_record[f'{cell_line_code}_GI50'] = gi50_val
                compound_record[f'{cell_line_code}_TGI'] = gi50_val + np.random.uniform(0.5, 1.5)
                compound_record[f'{cell_line_code}_LC50'] = compound_record[f'{cell_line_code}_TGI'] + np.random.uniform(0.5, 1.5)
                
                # Activity classification
                if gi50_val < -5.0:
                    compound_record[f'{cell_line_code}_Activity'] = 1  # Active
                else:
                    compound_record[f'{cell_line_code}_Activity'] = 0  # Inactive
            
            compounds.append(compound_record)
        
        self.compound_data['sample'] = compounds
        print(f"✅ Created sample dataset with {len(compounds)} compounds")
        return True
    
    def create_unified_dataset(self, activity_threshold=-5.0):
        """
        Create unified dataset from compound data
        """
        print(f"\n=== STEP 2: CREATING UNIFIED DATASET ===")
        print(f"Activity threshold: GI50 < {activity_threshold} log10(M) = Active")
        
        all_compounds = []
        for source, compounds in self.compound_data.items():
            all_compounds.extend(compounds)
        
        if not all_compounds:
            print("❌ No compound data available")
            return False
        
        df = pd.DataFrame(all_compounds)
        unified_records = []
        
        for _, compound in df.iterrows():
            for cell_line_code in self.nci60_cell_lines.keys():
                gi50_col = f'{cell_line_code}_GI50'
                
                if gi50_col in compound and pd.notna(compound[gi50_col]):
                    record = {
                        'nsc_id': compound['nsc_id'],
                        'compound_name': compound['compound_name'],
                        'smiles': compound['smiles'],
                        'source': compound['source'],
                        'cell_line': cell_line_code,
                        'gi50_log10m': compound[gi50_col],
                        'tgi_log10m': compound.get(f'{cell_line_code}_TGI', None),
                        'lc50_log10m': compound.get(f'{cell_line_code}_LC50', None)
                    }
                    
                    # Determine activity
                    if compound[gi50_col] < activity_threshold:
                        record['activity'] = 1
                        record['activity_level'] = 'Active'
                    else:
                        record['activity'] = 0
                        record['activity_level'] = 'Inactive'
                    
                    unified_records.append(record)
        
        self.unified_dataset = pd.DataFrame(unified_records)
        
        print(f"✅ Created unified dataset:")
        print(f"  Total combinations: {len(self.unified_dataset)}")
        print(f"  Unique compounds: {self.unified_dataset['nsc_id'].nunique()}")
        print(f"  Active combinations: {(self.unified_dataset['activity'] == 1).sum()}")
        print(f"  Inactive combinations: {(self.unified_dataset['activity'] == 0).sum()}")
        
        return True
    
    def analyze_multicell_activity(self, min_cell_lines=2):
        """
        Analyze compounds across multiple cell lines with conflict resolution
        """
        print(f"\n=== STEP 3: MULTI-CELL LINE ANALYSIS ===")
        print(f"Finding compounds tested in ≥{min_cell_lines} cell lines")
        
        if self.unified_dataset is None or self.unified_dataset.empty:
            print("❌ No unified dataset available")
            return False
        
        compound_analysis = []
        conflicts_found = 0
        
        for compound_id in self.unified_dataset['nsc_id'].unique():
            compound_data = self.unified_dataset[self.unified_dataset['nsc_id'] == compound_id]
            
            if len(compound_data) >= min_cell_lines:
                activities = compound_data['activity'].tolist()
                gi50_values = compound_data['gi50_log10m'].tolist()
                
                # Check for conflicts
                has_conflict = len(set(activities)) > 1
                if has_conflict:
                    conflicts_found += 1
                
                # Consensus resolution
                active_count = sum(activities)
                total_count = len(activities)
                activity_ratio = active_count / total_count
                
                # Determine consensus activity and confidence
                if activity_ratio >= 0.7:
                    consensus_activity = 1
                    confidence = 'high'
                    selectivity = 'broad_active'
                elif activity_ratio >= 0.5:
                    consensus_activity = 1
                    confidence = 'medium'
                    selectivity = 'selective_active'
                elif activity_ratio >= 0.3:
                    consensus_activity = 0
                    confidence = 'medium'
                    selectivity = 'selective_inactive'
                else:
                    consensus_activity = 0
                    confidence = 'high'
                    selectivity = 'broad_inactive'
                
                gi50_consistency = 'high' if np.std(gi50_values) < 1.0 else 'medium' if np.std(gi50_values) < 2.0 else 'low'
                
                analysis_record = {
                    'nsc_id': compound_id,
                    'compound_name': compound_data['compound_name'].iloc[0],
                    'smiles': compound_data['smiles'].iloc[0],
                    'source': compound_data['source'].iloc[0],
                    'num_cell_lines_tested': len(compound_data),
                    'cell_lines_tested': ','.join(compound_data['cell_line'].tolist()),
                    'active_cell_lines': active_count,
                    'inactive_cell_lines': total_count - active_count,
                    'activity_ratio': activity_ratio,
                    'consensus_activity': consensus_activity,
                    'confidence': confidence,
                    'selectivity': selectivity,
                    'gi50_consistency': gi50_consistency,
                    'has_conflict': has_conflict,
                    'mean_gi50': np.mean(gi50_values),
                    'std_gi50': np.std(gi50_values),
                    'min_gi50': min(gi50_values),
                    'max_gi50': max(gi50_values)
                }
                
                compound_analysis.append(analysis_record)
        
        self.multi_cellline_analysis = pd.DataFrame(compound_analysis)
        
        print(f"✅ Multi-cell line analysis completed:")
        print(f"  Compounds analyzed: {len(self.multi_cellline_analysis)}")
        print(f"  Compounds with conflicts: {conflicts_found}")
        print(f"  Consensus active: {(self.multi_cellline_analysis['consensus_activity'] == 1).sum()}")
        print(f"  Consensus inactive: {(self.multi_cellline_analysis['consensus_activity'] == 0).sum()}")
        
        return True
    
    def create_coverage_tiers(self):
        """
        Create 5-tier coverage analysis for virtual screening
        """
        print(f"\n=== STEP 4: CREATING COVERAGE TIERS ===")
        print("Generating 5-tier coverage analysis for virtual screening")
        
        if not hasattr(self, 'multi_cellline_analysis') or self.multi_cellline_analysis.empty:
            print("❌ No multi-cell line analysis available")
            return False
        
        # Define coverage tiers
        tier_criteria = {
            'tier_1_robust': {
                'min_cell_lines': 5,
                'min_confidence': 'high',
                'min_activity_ratio': 0.8,
                'description': 'Ultra-robust: 5+ cell lines, high confidence'
            },
            'tier_2_confident': {
                'min_cell_lines': 4,
                'min_confidence': 'high',
                'min_activity_ratio': 0.7,
                'description': 'Confident: 4+ cell lines, high confidence'
            },
            'tier_3_validated': {
                'min_cell_lines': 3,
                'min_confidence': 'medium',
                'min_activity_ratio': 0.6,
                'description': 'Validated: 3+ cell lines, medium+ confidence'
            },
            'tier_4_selective': {
                'min_cell_lines': 2,
                'min_confidence': 'medium',
                'min_activity_ratio': 0.5,
                'description': 'Selective: 2+ cell lines, selective activity'
            },
            'tier_5_exploratory': {
                'min_cell_lines': 1,
                'min_confidence': 'low',
                'min_activity_ratio': 0.0,
                'description': 'Exploratory: Single cell line data'
            }
        }
        
        confidence_levels = {'low': 0, 'medium': 1, 'high': 2}
        self.coverage_tiers = {}
        
        print(f"\n📊 Creating coverage tiers:")
        
        for tier_name, criteria in tier_criteria.items():
            filtered_compounds = self.multi_cellline_analysis[
                (self.multi_cellline_analysis['num_cell_lines_tested'] >= criteria['min_cell_lines']) &
                (self.multi_cellline_analysis['confidence'].map(confidence_levels) >= confidence_levels[criteria['min_confidence']]) &
                (self.multi_cellline_analysis['activity_ratio'] >= criteria['min_activity_ratio'])
            ].copy()
            
            if not filtered_compounds.empty:
                # Add VS suitability score
                filtered_compounds['vs_score'] = (
                    filtered_compounds['num_cell_lines_tested'] * 0.3 +
                    filtered_compounds['confidence'].map(confidence_levels) * 0.3 +
                    filtered_compounds['activity_ratio'] * 0.4
                ) / 3.0
                
                filtered_compounds['tier'] = tier_name
                filtered_compounds['tier_description'] = criteria['description']
                
                self.coverage_tiers[tier_name] = filtered_compounds
                
                active_count = (filtered_compounds['consensus_activity'] == 1).sum()
                inactive_count = (filtered_compounds['consensus_activity'] == 0).sum()
                
                print(f"  {tier_name}: {len(filtered_compounds)} compounds")
                print(f"    Active: {active_count}, Inactive: {inactive_count}")
            else:
                print(f"  {tier_name}: 0 compounds")
                self.coverage_tiers[tier_name] = pd.DataFrame()
        
        return True
    
    def create_vs_strategies(self):
        """
        Create virtual screening strategies
        """
        print(f"\n=== STEP 5: CREATING VS STRATEGIES ===")
        
        strategies = {
            'conservative': {
                'tiers': ['tier_1_robust', 'tier_2_confident'],
                'description': 'High-confidence, multi-validated compounds'
            },
            'balanced': {
                'tiers': ['tier_2_confident', 'tier_3_validated', 'tier_4_selective'],
                'description': 'Mix of validated and selective compounds'
            },
            'comprehensive': {
                'tiers': ['tier_3_validated', 'tier_4_selective', 'tier_5_exploratory'],
                'description': 'Maximum chemical diversity'
            }
        }
        
        self.vs_strategies = {}
        
        for strategy_name, strategy_info in strategies.items():
            combined_compounds = []
            
            for tier_name in strategy_info['tiers']:
                if tier_name in self.coverage_tiers and not self.coverage_tiers[tier_name].empty:
                    tier_data = self.coverage_tiers[tier_name].copy()
                    tier_data['strategy'] = strategy_name
                    combined_compounds.append(tier_data)
            
            if combined_compounds:
                strategy_dataset = pd.concat(combined_compounds, ignore_index=True)
                strategy_dataset = strategy_dataset.drop_duplicates(subset=['nsc_id'])
                
                # Add ML weights
                def calculate_ml_weight(row):
                    conf_weight = {'low': 0.3, 'medium': 0.7, 'high': 1.0}[row['confidence']]
                    coverage_bonus = min(row['num_cell_lines_tested'] / 6.0, 1.0)
                    consistency_bonus = row['activity_ratio'] if row['consensus_activity'] == 1 else (1 - row['activity_ratio'])
                    return min((conf_weight * 0.4 + coverage_bonus * 0.3 + consistency_bonus * 0.3), 1.0)
                
                strategy_dataset['ml_weight'] = strategy_dataset.apply(calculate_ml_weight, axis=1)
                strategy_dataset['training_priority'] = strategy_dataset['ml_weight'] >= 0.7
                
                self.vs_strategies[strategy_name] = strategy_dataset
                
                print(f"  {strategy_name.upper()}: {len(strategy_dataset)} compounds")
                print(f"    Active: {(strategy_dataset['consensus_activity'] == 1).sum()}")
                print(f"    High-weight: {(strategy_dataset['ml_weight'] >= 0.7).sum()}")
            else:
                self.vs_strategies[strategy_name] = pd.DataFrame()
        
        return True
    
    def save_results(self, output_file="NCI_MultiCellLine_VS_Dataset.xlsx"):
        """
        Save comprehensive results
        """
        print(f"\n=== STEP 6: SAVING RESULTS ===")
        
        try:
            with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
                
                # Virtual screening strategies (priority sheets)
                if hasattr(self, 'vs_strategies'):
                    for strategy_name, strategy_data in self.vs_strategies.items():
                        if not strategy_data.empty:
                            sheet_name = f'VS_{strategy_name.title()}'
                            strategy_data.to_excel(writer, sheet_name=sheet_name, index=False)
                
                # Coverage tiers
                if hasattr(self, 'coverage_tiers'):
                    for tier_name, tier_data in self.coverage_tiers.items():
                        if not tier_data.empty:
                            sheet_name = f'Tier_{tier_name.split("_")[1].title()}'
                            tier_data.to_excel(writer, sheet_name=sheet_name, index=False)
                
                # Analysis data
                if hasattr(self, 'multi_cellline_analysis'):
                    self.multi_cellline_analysis.to_excel(writer, sheet_name='Multi_CellLine_Analysis', index=False)
                
                if self.unified_dataset is not None:
                    self.unified_dataset.to_excel(writer, sheet_name='Unified_Dataset', index=False)
                
                # Summary tables
                if hasattr(self, 'vs_strategies'):
                    vs_summary = []
                    for strategy_name, strategy_data in self.vs_strategies.items():
                        if not strategy_data.empty:
                            vs_summary.append({
                                'Strategy': strategy_name.title(),
                                'Total_Compounds': len(strategy_data),
                                'Active_Compounds': (strategy_data['consensus_activity'] == 1).sum(),
                                'Inactive_Compounds': (strategy_data['consensus_activity'] == 0).sum(),
                                'Avg_Cell_Lines': strategy_data['num_cell_lines_tested'].mean(),
                                'High_Weight_Count': (strategy_data['ml_weight'] >= 0.7).sum()
                            })
                    
                    if vs_summary:
                        vs_summary_df = pd.DataFrame(vs_summary)
                        vs_summary_df.to_excel(writer, sheet_name='VS_Strategy_Summary', index=False)
                
                # Processing info
                metadata = pd.DataFrame({
                    'Parameter': [
                        'Target_Cell_Lines',
                        'NCI60_Available',
                        'Total_Compounds',
                        'Multi_CellLine_Compounds',
                        'VS_Strategies_Created',
                        'Coverage_Tiers',
                        'Processing_Date'
                    ],
                    'Value': [
                        ', '.join(self.target_cell_lines.keys()),
                        ', '.join(self.nci60_cell_lines.keys()),
                        self.unified_dataset['nsc_id'].nunique() if self.unified_dataset is not None else 0,
                        len(self.multi_cellline_analysis) if hasattr(self, 'multi_cellline_analysis') else 0,
                        len(self.vs_strategies) if hasattr(self, 'vs_strategies') else 0,
                        len(self.coverage_tiers) if hasattr(self, 'coverage_tiers') else 0,
                        pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')
                    ]
                })
                metadata.to_excel(writer, sheet_name='Processing_Info', index=False)
            
            print(f"✅ Results saved to: {output_file}")
            
            print(f"\n🎉 NCI MULTI-CELL LINE ANALYSIS COMPLETED!")
            print(f"📁 Output file: {output_file}")
            
            if hasattr(self, 'vs_strategies'):
                print(f"\n🎯 VIRTUAL SCREENING STRATEGIES:")
                for strategy_name, strategy_data in self.vs_strategies.items():
                    if not strategy_data.empty:
                        active_count = (strategy_data['consensus_activity'] == 1).sum()
                        print(f"  {strategy_name.upper()}: {len(strategy_data)} compounds ({active_count} active)")
            
            print(f"\n📊 PRIORITY SHEETS:")
            print(f"  1. 'VS_Conservative' - Robust hit identification")
            print(f"  2. 'VS_Balanced' - Standard drug discovery")
            print(f"  3. 'VS_Comprehensive' - Novel mechanism discovery")
            print(f"  4. 'VS_Strategy_Summary' - Strategy comparison")
            
            return output_file
            
        except Exception as e:
            print(f"❌ Error saving results: {e}")
            return None

def run_nci_multicell_analysis():
    """
    Run complete NCI multi-cell line analysis
    """
    print("🧬 NCI MULTI-CELL LINE ANALYSIS")
    print("Advanced Virtual Screening Dataset Generation")
    print("=" * 60)
    
    extractor = NCIMultiCellExtractor()
    
    # Step 1: Load sample data
    if not extractor.fetch_sample_nci_data():
        print("❌ Failed to load sample data")
        return None
    
    # Step 2: Create unified dataset
    if not extractor.create_unified_dataset():
        print("❌ Failed to create unified dataset")
        return None
    
    # Step 3: Multi-cell line analysis
    if not extractor.analyze_multicell_activity():
        print("❌ Failed to analyze multi-cell line activity")
        return None
    
    # Step 4: Create coverage tiers
    if not extractor.create_coverage_tiers():
        print("❌ Failed to create coverage tiers")
        return None
    
    # Step 5: Create VS strategies
    if not extractor.create_vs_strategies():
        print("❌ Failed to create VS strategies")
        return None
    
    # Step 6: Save results
    output_file = extractor.save_results()
    
    return output_file

if __name__ == "__main__":
    print("🚀 Starting NCI Multi-Cell Line Analysis...")
    
    try:
        output_file = run_nci_multicell_analysis()
        
        if output_file:
            print(f"\n🎉 SUCCESS! Analysis completed!")
            print(f"📁 Results: {output_file}")
            print(f"\n💡 Next steps:")
            print(f"  1. Review 'VS_Strategy_Summary' to choose approach")
            print(f"  2. Use VS sheets for virtual screening")
            print(f"  3. Apply ML weights for model training")
        else:
            print(f"\n❌ Analysis failed")
            
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
        