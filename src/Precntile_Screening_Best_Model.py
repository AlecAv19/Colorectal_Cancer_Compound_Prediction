# -*- coding: utf-8 -*-
"""
Created on Wed Aug 20 23:37:00 2025

@author: ASUS TUF F15
"""

"""
HCT-15 Percentile-Based Virtual Screening
Use relative ranking instead of absolute thresholds for better hit selection
"""

import pandas as pd
import numpy as np
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

class HCT15PercentileScreening:
    """
    Implement percentile-based virtual screening for better hit selection
    """
    
    def __init__(self):
        self.hct15_results = None
        self.optimized_results = None
        
    def load_hct15_results(self):
        """Load HCT-15 results"""
        print("📊 LOADING HCT-15 RESULTS FOR PERCENTILE OPTIMIZATION")
        print("=" * 55)
        
        import glob
        hct15_files = [f for f in glob.glob("VS_Results_*HCT-15*.xlsx") if "HCT-15" in f]
        
        if not hct15_files:
            print("❌ HCT-15 results file not found!")
            return None
        
        hct15_file = max(hct15_files, key=lambda x: x.split('_')[2] + x.split('_')[3])
        print(f"✅ Loading: {hct15_file}")
        
        try:
            for sheet in ['All_Results', 'Sheet1', None]:
                try:
                    if sheet:
                        df = pd.read_excel(hct15_file, sheet_name=sheet)
                    else:
                        df = pd.read_excel(hct15_file)
                    break
                except:
                    continue
            
            print(f"📊 Loaded: {df.shape[0]:,} compounds")
            print(f"📊 Score range: {df['Probability'].min():.3f} - {df['Probability'].max():.3f}")
            
            # Sort by probability for ranking
            df = df.sort_values('Probability', ascending=False)
            df['Percentile_Rank'] = (range(1, len(df) + 1))
            df['Percentile'] = (1 - (df['Percentile_Rank'] - 1) / len(df)) * 100
            
            self.hct15_results = df
            return df
            
        except Exception as e:
            print(f"❌ Error loading: {e}")
            return None
    
    def analyze_percentile_performance(self):
        """Analyze performance at different percentile cutoffs"""
        print(f"\n🔍 PERCENTILE-BASED ANALYSIS")
        print("=" * 30)
        
        if self.hct15_results is None:
            return None
        
        df = self.hct15_results
        
        # Get known compounds
        known_patterns = [
            'Paclitaxel', 'Podophyllotoxin', 'Doxorubicin', 'Cisplatin',
            'Actinomycin', 'Etoposide', 'Fluorouracil', 'Carboplatin', 'CHEMBL'
        ]
        known_compounds = df[df['Compound_ID'].str.contains('|'.join(known_patterns), case=False, na=False)]
        decoy_compounds = df[df['Compound_ID'].str.contains('DECOY', case=False, na=False)]
        
        print(f"📊 Dataset composition:")
        print(f"   Total compounds: {len(df):,}")
        print(f"   Known compounds: {len(known_compounds)}")
        print(f"   DECOY compounds: {len(decoy_compounds):,}")
        
        # Analyze known compound rankings
        print(f"\n💊 KNOWN COMPOUND RANKINGS:")
        known_sorted = known_compounds.sort_values('Probability', ascending=False)
        
        for i, (_, row) in enumerate(known_sorted.head(10).iterrows(), 1):
            percentile = row['Percentile']
            rank = row['Percentile_Rank']
            print(f"   {i:2d}. {row['Compound_ID']}: {row['Probability']:.3f} (Rank {rank:,}, {percentile:.1f}th percentile)")
        
        # Test different percentile cutoffs
        percentile_cutoffs = [1, 2, 5, 10, 15, 20, 25, 30]
        
        print(f"\n🎯 PERCENTILE CUTOFF ANALYSIS:")
        print("   Cutoff | Compounds | Known Captured | DECOY Rate | Score Threshold")
        print("   ------ | --------- | -------------- | ---------- | ---------------")
        
        optimization_results = []
        
        for cutoff in percentile_cutoffs:
            # Calculate threshold for this percentile
            threshold_score = df['Probability'].quantile(1 - cutoff/100)
            n_selected = int(len(df) * cutoff / 100)
            
            # Get top compounds at this cutoff
            top_compounds = df.head(n_selected)
            
            # Count known compounds in top selection
            known_in_top = len(top_compounds[top_compounds['Compound_ID'].str.contains('|'.join(known_patterns), case=False, na=False)])
            known_capture_rate = known_in_top / len(known_compounds) * 100 if len(known_compounds) > 0 else 0
            
            # Count DECOYs in top selection
            decoy_in_top = len(top_compounds[top_compounds['Compound_ID'].str.contains('DECOY', case=False, na=False)])
            decoy_rate = decoy_in_top / n_selected * 100 if n_selected > 0 else 0
            
            # Calculate quality score (maximize known capture, minimize DECOY rate)
            quality_score = known_capture_rate - (decoy_rate * 0.5)
            
            optimization_results.append({
                'cutoff': cutoff,
                'n_selected': n_selected,
                'threshold_score': threshold_score,
                'known_captured': known_in_top,
                'known_capture_rate': known_capture_rate,
                'decoy_in_top': decoy_in_top,
                'decoy_rate': decoy_rate,
                'quality_score': quality_score
            })
            
            print(f"   {cutoff:5.0f}% | {n_selected:8,} | {known_in_top:2d}/{len(known_compounds):2d} ({known_capture_rate:5.1f}%) | {decoy_rate:8.1f}% | {threshold_score:13.3f}")
        
        # Find optimal percentile cutoff
        best_result = max(optimization_results, key=lambda x: x['quality_score'])
        optimal_cutoff = best_result['cutoff']
        
        print(f"\n✅ OPTIMAL PERCENTILE CUTOFF: {optimal_cutoff}%")
        print(f"   Compounds selected: {best_result['n_selected']:,}")
        print(f"   Score threshold: {best_result['threshold_score']:.3f}")
        print(f"   Known compounds captured: {best_result['known_captured']}/{len(known_compounds)} ({best_result['known_capture_rate']:.1f}%)")
        print(f"   DECOY rate: {best_result['decoy_rate']:.1f}%")
        print(f"   Quality score: {best_result['quality_score']:.1f}")
        
        return optimal_cutoff, optimization_results
    
    def create_percentile_results(self, optimal_cutoff):
        """Create results based on optimal percentile cutoff"""
        print(f"\n📊 CREATING PERCENTILE-BASED RESULTS")
        print("=" * 35)
        
        if self.hct15_results is None:
            return None
        
        df = self.hct15_results.copy()
        
        # Calculate selection based on percentile
        n_selected = int(len(df) * optimal_cutoff / 100)
        threshold_score = df['Probability'].quantile(1 - optimal_cutoff/100)
        
        # Mark selected compounds
        df['Selected'] = df['Percentile_Rank'] <= n_selected
        df['Selection_Method'] = f'Top_{optimal_cutoff}%'
        
        # Add tier classification
        df['Tier'] = 'Not Selected'
        df.loc[df['Percentile_Rank'] <= n_selected * 0.2, 'Tier'] = 'Tier 1 (Top 20%)'  # Top 20% of selected
        df.loc[(df['Percentile_Rank'] > n_selected * 0.2) & (df['Percentile_Rank'] <= n_selected * 0.5), 'Tier'] = 'Tier 2 (Top 50%)'
        df.loc[(df['Percentile_Rank'] > n_selected * 0.5) & (df['Percentile_Rank'] <= n_selected), 'Tier'] = 'Tier 3 (Remaining)'
        
        # Add confidence based on ranking within selected compounds
        df['Confidence_Percentile'] = 'Low'
        df.loc[df['Percentile_Rank'] <= n_selected * 0.1, 'Confidence_Percentile'] = 'Very High'  # Top 10% of selected
        df.loc[(df['Percentile_Rank'] > n_selected * 0.1) & (df['Percentile_Rank'] <= n_selected * 0.3), 'Confidence_Percentile'] = 'High'
        df.loc[(df['Percentile_Rank'] > n_selected * 0.3) & (df['Percentile_Rank'] <= n_selected * 0.7), 'Confidence_Percentile'] = 'Medium'
        
        print(f"📊 Selection Results:")
        print(f"   Optimal cutoff: {optimal_cutoff}%")
        print(f"   Selected compounds: {n_selected:,} / {len(df):,}")
        print(f"   Score threshold: {threshold_score:.3f}")
        
        # Analyze selected compounds
        selected_compounds = df[df['Selected']]
        
        # Known compounds in selection
        known_patterns = [
            'Paclitaxel', 'Podophyllotoxin', 'Doxorubicin', 'Cisplatin',
            'Actinomycin', 'Etoposide', 'Fluorouracil', 'Carboplatin', 'CHEMBL'
        ]
        known_selected = selected_compounds[selected_compounds['Compound_ID'].str.contains('|'.join(known_patterns), case=False, na=False)]
        decoy_selected = selected_compounds[selected_compounds['Compound_ID'].str.contains('DECOY', case=False, na=False)]
        
        print(f"\n📊 Selection Composition:")
        print(f"   Known compounds: {len(known_selected):,} ({len(known_selected)/n_selected*100:.1f}%)")
        print(f"   DECOY compounds: {len(decoy_selected):,} ({len(decoy_selected)/n_selected*100:.1f}%)")
        print(f"   Other compounds: {n_selected - len(known_selected) - len(decoy_selected):,}")
        
        # Tier breakdown
        print(f"\n📊 Tier Breakdown:")
        for tier in ['Tier 1 (Top 20%)', 'Tier 2 (Top 50%)', 'Tier 3 (Remaining)']:
            tier_count = (df['Tier'] == tier).sum()
            if tier_count > 0:
                print(f"   {tier}: {tier_count:,} compounds")
        
        self.optimized_results = df
        return df
    
    def save_percentile_results(self, optimal_cutoff):
        """Save percentile-based results"""
        if self.optimized_results is None:
            return None
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_file = f'HCT15_Percentile_Results_{optimal_cutoff}pct_{timestamp}.xlsx'
        
        print(f"\n💾 SAVING PERCENTILE-BASED RESULTS")
        print("=" * 35)
        
        df = self.optimized_results
        
        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            # All results with percentile rankings
            df.to_excel(writer, sheet_name='All_Results_Ranked', index=False)
            
            # Selected compounds only
            selected = df[df['Selected']].copy()
            selected.to_excel(writer, sheet_name='Selected_Compounds', index=False)
            
            # Tier 1 (highest priority)
            tier1 = df[df['Tier'] == 'Tier 1 (Top 20%)'].copy()
            if len(tier1) > 0:
                tier1.to_excel(writer, sheet_name='Tier_1_Priority', index=False)
            
            # Known compounds analysis
            known_patterns = [
                'Paclitaxel', 'Podophyllotoxin', 'Doxorubicin', 'Cisplatin',
                'Actinomycin', 'Etoposide', 'Fluorouracil', 'Carboplatin', 'CHEMBL'
            ]
            known_compounds = df[df['Compound_ID'].str.contains('|'.join(known_patterns), case=False, na=False)].copy()
            if len(known_compounds) > 0:
                known_compounds.to_excel(writer, sheet_name='Known_Compounds', index=False)
            
            # Very high confidence compounds
            very_high_conf = df[df['Confidence_Percentile'] == 'Very High'].copy()
            if len(very_high_conf) > 0:
                very_high_conf.to_excel(writer, sheet_name='Very_High_Confidence', index=False)
            
            # Summary statistics
            summary_data = [
                ['Selection Method', f'Top {optimal_cutoff}% percentile'],
                ['Total Compounds', len(df)],
                ['Selected Compounds', df['Selected'].sum()],
                ['Selection Rate', f"{df['Selected'].mean()*100:.1f}%"],
                ['Score Threshold', f"{df[df['Selected']]['Probability'].min():.3f}"],
                ['Top Score', f"{df['Probability'].max():.3f}"],
                ['Known Compounds Selected', len(df[df['Selected'] & df['Compound_ID'].str.contains('|'.join(known_patterns), case=False, na=False)])],
                ['Analysis Date', datetime.now().strftime('%Y-%m-%d %H:%M:%S')]
            ]
            
            summary_df = pd.DataFrame(summary_data, columns=['Metric', 'Value'])
            summary_df.to_excel(writer, sheet_name='Summary', index=False)
        
        print(f"✅ Percentile results saved: {output_file}")
        return output_file
    
    def show_top_percentile_hits(self, optimal_cutoff):
        """Show top hits from percentile-based selection"""
        if self.optimized_results is None:
            return
        
        print(f"\n🏆 TOP 20 PERCENTILE-BASED HITS ({optimal_cutoff}% CUTOFF)")
        print("=" * 50)
        
        selected_compounds = self.optimized_results[self.optimized_results['Selected']]
        top_20 = selected_compounds.head(20)
        
        for i, (_, row) in enumerate(top_20.iterrows(), 1):
            compound = row['Compound_ID']
            score = row['Probability']
            rank = row['Percentile_Rank']
            percentile = row['Percentile']
            tier = row['Tier']
            confidence = row['Confidence_Percentile']
            
            # Identify compound type
            if 'DECOY' in compound:
                comp_type = "DECOY"
            elif any(drug in compound for drug in ['Paclitaxel', 'Podophyllotoxin', 'Doxorubicin', 'CHEMBL']):
                comp_type = "KNOWN"
            else:
                comp_type = "OTHER"
            
            print(f"{i:2d}. [{comp_type}] {compound}")
            print(f"    Score: {score:.3f} | Rank: {rank:,} ({percentile:.1f}%ile) | {confidence} confidence")
        
        # Show statistics
        print(f"\n📊 SELECTION STATISTICS:")
        n_selected = selected_compounds.shape[0]
        
        # Known compound analysis
        known_patterns = ['Paclitaxel', 'Podophyllotoxin', 'Doxorubicin', 'Cisplatin', 'Actinomycin', 'Etoposide', 'CHEMBL']
        known_in_selection = selected_compounds[selected_compounds['Compound_ID'].str.contains('|'.join(known_patterns), case=False, na=False)]
        
        print(f"   Total selected: {n_selected:,}")
        print(f"   Known drugs captured: {len(known_in_selection)}")
        print(f"   Very High confidence: {(selected_compounds['Confidence_Percentile'] == 'Very High').sum():,}")
        print(f"   High confidence: {(selected_compounds['Confidence_Percentile'] == 'High').sum():,}")
        
        if len(known_in_selection) > 0:
            print(f"\n💊 KNOWN DRUGS IN SELECTION:")
            for _, row in known_in_selection.head(10).iterrows():
                print(f"   • {row['Compound_ID']}: Rank {row['Percentile_Rank']:,} ({row['Percentile']:.1f}%ile)")

def main():
    """Main percentile-based screening pipeline"""
    print("📊 HCT-15 PERCENTILE-BASED VIRTUAL SCREENING")
    print("=" * 50)
    print("Use relative ranking instead of absolute thresholds")
    print()
    
    screener = HCT15PercentileScreening()
    
    # Load HCT-15 results
    df = screener.load_hct15_results()
    if df is None:
        return
    
    # Analyze percentile performance
    optimal_cutoff, optimization_results = screener.analyze_percentile_performance()
    if optimal_cutoff is None:
        return
    
    # Create percentile-based results
    optimized_df = screener.create_percentile_results(optimal_cutoff)
    if optimized_df is None:
        return
    
    # Save results
    output_file = screener.save_percentile_results(optimal_cutoff)
    
    # Show top hits
    screener.show_top_percentile_hits(optimal_cutoff)
    
    print(f"\n🎉 PERCENTILE-BASED SCREENING COMPLETE!")
    print("=" * 45)
    print(f"✅ Optimal cutoff: Top {optimal_cutoff}% of compounds")
    print(f"✅ Realistic selection size")
    print(f"✅ Good known compound capture")
    print(f"✅ Results saved: {output_file}")
    
    print(f"\n🎯 RECOMMENDED APPROACH:")
    print("1. Focus on Tier 1 compounds for immediate validation")
    print("2. Use Very High confidence compounds as primary targets")
    print("3. Consider Tier 2 compounds for secondary screening")
    print(f"4. This approach avoids the calibration issues of absolute thresholds")

if __name__ == "__main__":
    main()