# -*- coding: utf-8 -*-
"""
Created on Wed Aug 20 22:59:16 2025

@author: ASUS TUF F15
"""
"""
Simple Multi-Model Virtual Screening Comparison
No graphics dependencies - pure pandas analysis
"""

import pandas as pd
import numpy as np
import os
import glob
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

class SimpleModelComparison:
    """
    Compare virtual screening results from multiple models without graphics
    """
    
    def __init__(self):
        self.results_data = {}
        self.comparison_summary = None
        self.cell_line_info = {
            'RKO': {'type': 'Colon', 'description': 'Colon carcinoma'},
            'HCT-15': {'type': 'Colon', 'description': 'Colon carcinoma'},
            'SW620': {'type': 'Colon', 'description': 'Colon adenocarcinoma (metastatic)'},
            'LoVo': {'type': 'Colon', 'description': 'Colon adenocarcinoma'},
            'DLD-1': {'type': 'Colon', 'description': 'Colon adenocarcinoma'},
            'Colo-205': {'type': 'Colon', 'description': 'Colon adenocarcinoma'},
            'Caco-2': {'type': 'Colon', 'description': 'Colon adenocarcinoma'},
            'SW480': {'type': 'Colon', 'description': 'Colon adenocarcinoma'}
        }
    
    def find_and_load_results(self):
        """Find and load all VS results files"""
        print("🔍 MULTI-MODEL VS RESULTS COMPARISON")
        print("=" * 45)
        
        # Look for VS results files
        result_files = glob.glob("VS_Results_*.xlsx")
        result_files.sort()
        
        print(f"📂 Found {len(result_files)} results files:")
        
        loaded_results = {}
        
        for i, file in enumerate(result_files, 1):
            file_size = os.path.getsize(file) / 1024  # KB
            mod_time = datetime.fromtimestamp(os.path.getmtime(file))
            
            # Extract model identifier
            model_id = "Unknown"
            for cell_line in self.cell_line_info.keys():
                if cell_line in file:
                    model_id = cell_line
                    break
            
            print(f"  {i:2d}. {file}")
            print(f"      Model: {model_id} | Size: {file_size:.0f} KB | Modified: {mod_time.strftime('%H:%M')}")
            
            # Load the file
            try:
                # Try different sheet names
                sheet_options = ['All_Results', 'Sheet1', None]
                df = None
                
                for sheet in sheet_options:
                    try:
                        if sheet:
                            df = pd.read_excel(file, sheet_name=sheet)
                        else:
                            df = pd.read_excel(file)
                        break
                    except:
                        continue
                
                if df is None:
                    print(f"      ❌ Could not load")
                    continue
                
                # Validate required columns
                if 'Probability' not in df.columns or 'Compound_ID' not in df.columns:
                    print(f"      ❌ Missing required columns")
                    continue
                
                # Store results
                loaded_results[model_id] = {
                    'filename': file,
                    'data': df,
                    'model_id': model_id,
                    'cell_line_type': self.cell_line_info.get(model_id, {}).get('type', 'Unknown'),
                    'description': self.cell_line_info.get(model_id, {}).get('description', 'Unknown')
                }
                
                print(f"      ✅ Loaded: {df.shape[0]:,} compounds")
                
            except Exception as e:
                print(f"      ❌ Error: {e}")
        
        self.results_data = loaded_results
        print(f"\n✅ Successfully loaded {len(loaded_results)} model results")
        return loaded_results
    
    def analyze_all_models(self):
        """Comprehensive analysis of all models"""
        print(f"\n📊 COMPREHENSIVE MODEL ANALYSIS")
        print("=" * 35)
        
        analysis_results = []
        
        for model_id, result_info in self.results_data.items():
            print(f"\n🧬 ANALYZING: {model_id} - {result_info['description']}")
            print("-" * 60)
            
            df = result_info['data']
            
            # Basic statistics
            total_compounds = len(df)
            active_predictions = df['Active'].sum() if 'Active' in df.columns else (df['Probability'] >= 0.5).sum()
            active_percentage = (active_predictions / total_compounds) * 100
            
            mean_prob = df['Probability'].mean()
            median_prob = df['Probability'].median()
            std_prob = df['Probability'].std()
            min_prob = df['Probability'].min()
            max_prob = df['Probability'].max()
            prob_range = max_prob - min_prob
            
            print(f"📊 Basic Statistics:")
            print(f"   Total compounds: {total_compounds:,}")
            print(f"   Active predictions: {active_predictions:,} ({active_percentage:.1f}%)")
            print(f"   Score range: {min_prob:.3f} - {max_prob:.3f} (range: {prob_range:.3f})")
            print(f"   Mean/Median: {mean_prob:.3f} / {median_prob:.3f}")
            
            # Percentile thresholds
            p95 = df['Probability'].quantile(0.95)
            p90 = df['Probability'].quantile(0.90)
            p75 = df['Probability'].quantile(0.75)
            
            print(f"📈 Percentile Thresholds:")
            print(f"   95th percentile: {p95:.3f}")
            print(f"   90th percentile: {p90:.3f}")
            print(f"   75th percentile: {p75:.3f}")
            
            # Known compound analysis
            known_patterns = [
                'Paclitaxel', 'Podophyllotoxin', 'Doxorubicin', 'Cisplatin',
                'Actinomycin', 'Etoposide', 'Fluorouracil', 'Carboplatin', 'CHEMBL'
            ]
            
            known_compounds = df[df['Compound_ID'].str.contains('|'.join(known_patterns), case=False, na=False)]
            n_known = len(known_compounds)
            
            print(f"💊 Known Compounds Analysis:")
            if n_known > 0:
                known_active = known_compounds['Active'].sum() if 'Active' in known_compounds.columns else (known_compounds['Probability'] >= 0.5).sum()
                known_mean_score = known_compounds['Probability'].mean()
                known_max_score = known_compounds['Probability'].max()
                
                print(f"   Found: {n_known} known compounds")
                print(f"   Active: {known_active}/{n_known} ({known_active/n_known*100:.1f}%)")
                print(f"   Average score: {known_mean_score:.3f}")
                print(f"   Best score: {known_max_score:.3f}")
                
                # Show top known compounds
                top_known = known_compounds.nlargest(5, 'Probability')
                print(f"   Top known compounds:")
                for j, (_, row) in enumerate(top_known.iterrows(), 1):
                    status = "ACTIVE" if (row['Active'] if 'Active' in row else row['Probability'] >= 0.5) else "inactive"
                    print(f"     {j}. {row['Compound_ID']}: {row['Probability']:.3f} ({status})")
            else:
                known_mean_score = 0
                print(f"   No known compounds found")
            
            # DECOY analysis
            decoy_compounds = df[df['Compound_ID'].str.contains('DECOY', case=False, na=False)]
            n_decoy = len(decoy_compounds)
            
            print(f"🎯 DECOY Analysis:")
            if n_decoy > 0:
                decoy_active = decoy_compounds['Active'].sum() if 'Active' in decoy_compounds.columns else (decoy_compounds['Probability'] >= 0.5).sum()
                decoy_active_pct = (decoy_active / n_decoy * 100)
                decoy_mean_score = decoy_compounds['Probability'].mean()
                
                print(f"   DECOY compounds: {n_decoy:,}")
                print(f"   DECOY active: {decoy_active:,} ({decoy_active_pct:.1f}%)")
                print(f"   DECOY mean score: {decoy_mean_score:.3f}")
            else:
                decoy_active_pct = 0
                print(f"   No DECOY compounds found")
            
            # Top hits analysis
            print(f"🏆 Top 10 Hits:")
            top_hits = df.nlargest(10, 'Probability')
            for j, (_, row) in enumerate(top_hits.iterrows(), 1):
                compound = row['Compound_ID']
                score = row['Probability']
                status = "ACTIVE" if (row['Active'] if 'Active' in row else score >= 0.5) else "inactive"
                
                # Classify compound
                if 'DECOY' in compound:
                    comp_type = "DECOY"
                elif any(drug in compound for drug in known_patterns):
                    comp_type = "KNOWN"
                else:
                    comp_type = "OTHER"
                
                print(f"     {j:2d}. [{comp_type}] {compound}: {score:.3f} ({status})")
            
            # Store analysis results
            analysis_results.append({
                'Model_ID': model_id,
                'Cell_Line_Type': result_info['cell_line_type'],
                'Description': result_info['description'],
                'Total_Compounds': total_compounds,
                'Active_Predictions': active_predictions,
                'Active_Percentage': active_percentage,
                'Mean_Probability': mean_prob,
                'Median_Probability': median_prob,
                'Std_Probability': std_prob,
                'Min_Probability': min_prob,
                'Max_Probability': max_prob,
                'Probability_Range': prob_range,
                'P95_Threshold': p95,
                'P90_Threshold': p90,
                'P75_Threshold': p75,
                'Known_Compounds': n_known,
                'Known_Active': known_active if n_known > 0 else 0,
                'Known_Mean_Score': known_mean_score,
                'DECOY_Compounds': n_decoy,
                'DECOY_Active': decoy_active if n_decoy > 0 else 0,
                'DECOY_Active_Pct': decoy_active_pct,
                'Score_Diversity': std_prob / mean_prob if mean_prob > 0 else 0,
                'Filename': result_info['filename']
            })
        
        self.comparison_summary = pd.DataFrame(analysis_results)
        return self.comparison_summary
    
    def rank_models(self):
        """Rank models by quality metrics"""
        print(f"\n🏆 MODEL QUALITY RANKING")
        print("=" * 25)
        
        if self.comparison_summary is None:
            print("❌ No analysis data available")
            return None
        
        df = self.comparison_summary.copy()
        
        # Calculate quality scores
        print("🔧 Calculating quality scores...")
        print("   Criteria: Range + Known Performance + DECOY Discrimination + Hit Rate + Diversity")
        
        quality_data = []
        
        for _, row in df.iterrows():
            score = 0
            details = []
            
            # 1. Probability Range (25 points max)
            range_score = min(row['Probability_Range'] * 50, 25)
            score += range_score
            if row['Probability_Range'] > 0.3:
                details.append(f"Excellent range ({row['Probability_Range']:.3f})")
            elif row['Probability_Range'] > 0.2:
                details.append(f"Good range ({row['Probability_Range']:.3f})")
            elif row['Probability_Range'] > 0.1:
                details.append(f"Fair range ({row['Probability_Range']:.3f})")
            
            # 2. Known Compound Performance (30 points max)
            if row['Known_Compounds'] > 0:
                known_score = row['Known_Mean_Score'] * 30
                score += known_score
                if row['Known_Mean_Score'] > 0.7:
                    details.append(f"Excellent known drug scores ({row['Known_Mean_Score']:.3f})")
                elif row['Known_Mean_Score'] > 0.5:
                    details.append(f"Good known drug scores ({row['Known_Mean_Score']:.3f})")
                elif row['Known_Mean_Score'] > 0.4:
                    details.append(f"Fair known drug scores ({row['Known_Mean_Score']:.3f})")
            
            # 3. DECOY Discrimination (20 points max)
            decoy_score = max(0, 20 - row['DECOY_Active_Pct'])
            score += decoy_score
            if row['DECOY_Active_Pct'] < 2:
                details.append(f"Excellent DECOY discrimination ({row['DECOY_Active_Pct']:.1f}%)")
            elif row['DECOY_Active_Pct'] < 5:
                details.append(f"Good DECOY discrimination ({row['DECOY_Active_Pct']:.1f}%)")
            elif row['DECOY_Active_Pct'] < 10:
                details.append(f"Fair DECOY discrimination ({row['DECOY_Active_Pct']:.1f}%)")
            
            # 4. Hit Rate (15 points max)
            hit_rate = row['Active_Percentage']
            if 1 <= hit_rate <= 10:
                hit_score = 15
                details.append(f"Ideal hit rate ({hit_rate:.1f}%)")
            elif 0.5 <= hit_rate <= 15:
                hit_score = 10
                details.append(f"Good hit rate ({hit_rate:.1f}%)")
            elif hit_rate <= 25:
                hit_score = 5
                details.append(f"High hit rate ({hit_rate:.1f}%)")
            else:
                hit_score = 0
                details.append(f"Very high hit rate ({hit_rate:.1f}%)")
            score += hit_score
            
            # 5. Score Diversity (10 points max)
            diversity_score = min(row['Score_Diversity'] * 20, 10)
            score += diversity_score
            if row['Score_Diversity'] > 0.4:
                details.append("Good score diversity")
            
            quality_data.append({
                'Model_ID': row['Model_ID'],
                'Description': row['Description'],
                'Quality_Score': score,
                'Details': details,
                'Range_Score': range_score,
                'Known_Score': known_score if 'known_score' in locals() else 0,
                'DECOY_Score': decoy_score,
                'Hit_Rate_Score': hit_score,
                'Diversity_Score': diversity_score
            })
        
        # Create and sort ranking
        ranking_df = pd.DataFrame(quality_data)
        ranking_df = ranking_df.sort_values('Quality_Score', ascending=False)
        
        print(f"\n📊 FINAL RANKINGS:")
        print("=" * 50)
        
        for i, (_, row) in enumerate(ranking_df.iterrows(), 1):
            model_data = df[df['Model_ID'] == row['Model_ID']].iloc[0]
            
            print(f"\n🥇 RANK {i}: {row['Model_ID']}")
            print(f"   Cell Line: {row['Description']}")
            print(f"   Quality Score: {row['Quality_Score']:.1f}/100")
            print(f"   Hit Rate: {model_data['Active_Percentage']:.1f}% | Range: {model_data['Probability_Range']:.3f}")
            print(f"   Known Drugs: {model_data['Known_Mean_Score']:.3f} avg | DECOY Rate: {model_data['DECOY_Active_Pct']:.1f}%")
            print(f"   Strengths: {', '.join(row['Details']) if row['Details'] else 'None specific'}")
        
        return ranking_df
    
    def find_consensus_hits(self):
        """Find compounds that multiple models agree on"""
        print(f"\n🤝 CONSENSUS ANALYSIS")
        print("=" * 20)
        
        if len(self.results_data) < 2:
            print("❌ Need at least 2 models for consensus")
            return None
        
        # Find common compounds
        common_compounds = None
        for model_id, result_info in self.results_data.items():
            compound_set = set(result_info['data']['Compound_ID'])
            if common_compounds is None:
                common_compounds = compound_set
            else:
                common_compounds = common_compounds.intersection(compound_set)
        
        print(f"📊 Common compounds across all models: {len(common_compounds):,}")
        
        if len(common_compounds) < 100:
            print("❌ Too few common compounds for meaningful consensus")
            return None
        
        # Calculate consensus for each compound
        consensus_results = []
        
        for compound in common_compounds:
            scores = []
            predictions = []
            
            for model_id, result_info in self.results_data.items():
                df = result_info['data']
                compound_row = df[df['Compound_ID'] == compound]
                
                if len(compound_row) > 0:
                    score = compound_row['Probability'].iloc[0]
                    prediction = compound_row['Active'].iloc[0] if 'Active' in compound_row.columns else score >= 0.5
                    scores.append(score)
                    predictions.append(prediction)
            
            if len(scores) == len(self.results_data):
                consensus_results.append({
                    'Compound_ID': compound,
                    'Mean_Score': np.mean(scores),
                    'Std_Score': np.std(scores),
                    'Min_Score': np.min(scores),
                    'Max_Score': np.max(scores),
                    'Models_Active': sum(predictions),
                    'Total_Models': len(predictions),
                    'Agreement_Rate': sum(predictions) / len(predictions),
                    'Consensus_Active': sum(predictions) > len(predictions) / 2
                })
        
        consensus_df = pd.DataFrame(consensus_results)
        consensus_df = consensus_df.sort_values('Mean_Score', ascending=False)
        
        # Show consensus statistics
        high_consensus = consensus_df[consensus_df['Agreement_Rate'] >= 0.8]
        
        print(f"\n📊 CONSENSUS STATISTICS:")
        print(f"   Total analyzed: {len(consensus_df):,} compounds")
        print(f"   High consensus (≥80%): {len(high_consensus):,} compounds")
        print(f"   Consensus active: {consensus_df['Consensus_Active'].sum():,} compounds")
        
        # Show top consensus hits
        print(f"\n🎯 TOP 15 CONSENSUS HITS:")
        top_consensus = consensus_df.head(15)
        
        for i, (_, row) in enumerate(top_consensus.iterrows(), 1):
            compound = row['Compound_ID']
            mean_score = row['Mean_Score']
            agreement = row['Agreement_Rate']
            models_active = row['Models_Active']
            total_models = row['Total_Models']
            
            # Classify compound
            if 'DECOY' in compound:
                comp_type = "DECOY"
            elif any(drug in compound for drug in ['Paclitaxel', 'Podophyllotoxin', 'Doxorubicin', 'CHEMBL']):
                comp_type = "KNOWN"
            else:
                comp_type = "OTHER"
            
            status = "✅ CONSENSUS" if row['Consensus_Active'] else "❌ consensus"
            
            print(f"  {i:2d}. [{comp_type}] {compound}")
            print(f"      Score: {mean_score:.3f} | Agreement: {agreement:.0%} ({models_active}/{total_models}) | {status}")
        
        return consensus_df
    
    def create_final_recommendation(self):
        """Generate final recommendation"""
        print(f"\n🎯 FINAL RECOMMENDATION")
        print("=" * 25)
        
        ranking = self.rank_models()
        if ranking is None or len(ranking) == 0:
            print("❌ Cannot generate recommendation - no ranking data")
            return
        
        best_model = ranking.iloc[0]
        best_data = self.comparison_summary[self.comparison_summary['Model_ID'] == best_model['Model_ID']].iloc[0]
        
        print(f"🏆 RECOMMENDED MODEL: {best_model['Model_ID']}")
        print(f"   Cell Line: {best_model['Description']}")
        print(f"   Quality Score: {best_model['Quality_Score']:.1f}/100")
        print(f"   File: {best_data['Filename']}")
        
        print(f"\n📊 WHY THIS MODEL:")
        for detail in best_model['Details']:
            print(f"   ✅ {detail}")
        
        print(f"\n🎯 USAGE RECOMMENDATIONS:")
        if best_data['Active_Percentage'] < 5:
            print("   • Use 95th percentile as threshold for hit selection")
        elif best_data['Active_Percentage'] < 10:
            print("   • Current 0.5 threshold appears appropriate")
        else:
            print("   • Consider higher threshold (0.6-0.7) to reduce false positives")
        
        if best_data['Known_Mean_Score'] > 0.6:
            print("   • Model shows excellent performance on known drugs")
        elif best_data['Known_Mean_Score'] > 0.5:
            print("   • Model shows good performance on known drugs")
        else:
            print("   • ⚠ Model may need calibration - known drugs scoring low")
        
        if len(ranking) > 1:
            print(f"\n🥈 ALTERNATIVE OPTIONS:")
            for i in range(1, min(4, len(ranking))):
                alt = ranking.iloc[i]
                alt_data = self.comparison_summary[self.comparison_summary['Model_ID'] == alt['Model_ID']].iloc[0]
                print(f"   {i+1}. {alt['Model_ID']} (Score: {alt['Quality_Score']:.1f}) - {alt['Description']}")

def main():
    """Main comparison function"""
    print("🔬 SIMPLE MULTI-MODEL COMPARISON")
    print("=" * 40)
    print("Compare virtual screening results without graphics")
    print()
    
    comparator = SimpleModelComparison()
    
    # Load all results
    results = comparator.find_and_load_results()
    if not results:
        print("❌ No results files could be loaded!")
        return
    
    # Analyze all models
    analysis = comparator.analyze_all_models()
    
    # Rank models
    ranking = comparator.rank_models()
    
    # Find consensus
    consensus = comparator.find_consensus_hits()
    
    # Final recommendation
    comparator.create_final_recommendation()
    
    # Save summary
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    summary_file = f'Model_Comparison_Summary_{timestamp}.xlsx'
    
    print(f"\n💾 SAVING SUMMARY")
    print("=" * 15)
    
    with pd.ExcelWriter(summary_file, engine='openpyxl') as writer:
        if analysis is not None:
            analysis.to_excel(writer, sheet_name='Model_Analysis', index=False)
        if ranking is not None:
            ranking.to_excel(writer, sheet_name='Quality_Ranking', index=False)
        if consensus is not None:
            consensus.head(500).to_excel(writer, sheet_name='Consensus_Hits', index=False)
    
    print(f"✅ Summary saved: {summary_file}")
    
    print(f"\n🎉 COMPARISON COMPLETE!")
    print("=" * 25)
    print(f"✅ Analyzed {len(results)} models")
    print(f"✅ Generated quality rankings")
    print(f"✅ Identified best model for your research")

if __name__ == "__main__":
    main()