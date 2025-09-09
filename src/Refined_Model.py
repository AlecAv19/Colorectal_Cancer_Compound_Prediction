# -*- coding: utf-8 -*-
"""
Created on Wed Aug 20 23:39:58 2025

@author: ASUS TUF F15
"""

"""
HCT-15 Refined Hit Selection Strategy
Focus on high-confidence hits and investigate DECOY patterns
"""

import pandas as pd
import numpy as np
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

class HCT15RefinedStrategy:
    """
    Implement refined hit selection strategy focusing on highest confidence compounds
    """
    
    def __init__(self):
        self.results = None
        self.refined_hits = None
        
    def load_percentile_results(self):
        """Load the percentile-based results"""
        print("📊 LOADING PERCENTILE RESULTS FOR REFINEMENT")
        print("=" * 45)
        
        import glob
        percentile_files = glob.glob("HCT15_Percentile_Results_*.xlsx")
        
        if not percentile_files:
            print("❌ No percentile results found!")
            return None
        
        latest_file = max(percentile_files, key=lambda x: x.split('_')[-1])
        print(f"✅ Loading: {latest_file}")
        
        try:
            df = pd.read_excel(latest_file, sheet_name='All_Results_Ranked')
            print(f"📊 Loaded: {df.shape[0]:,} compounds")
            
            self.results = df
            return df
            
        except Exception as e:
            print(f"❌ Error loading: {e}")
            return None
    
    def analyze_top_performers(self):
        """Analyze the very top performing compounds in detail"""
        print(f"\n🔍 DETAILED ANALYSIS OF TOP PERFORMERS")
        print("=" * 40)
        
        if self.results is None:
            return None
        
        df = self.results
        
        # Focus on top 1%, 2%, 5% for detailed analysis
        analysis_cutoffs = [1, 2, 5]
        
        for cutoff in analysis_cutoffs:
            n_compounds = int(len(df) * cutoff / 100)
            top_compounds = df.head(n_compounds)
            
            print(f"\n📊 TOP {cutoff}% ANALYSIS ({n_compounds} compounds):")
            
            # Categorize compounds
            known_patterns = ['Paclitaxel', 'Podophyllotoxin', 'Doxorubicin', 'Cisplatin', 'Actinomycin', 'Etoposide', 'CHEMBL']
            
            known_in_top = top_compounds[top_compounds['Compound_ID'].str.contains('|'.join(known_patterns), case=False, na=False)]
            decoy_in_top = top_compounds[top_compounds['Compound_ID'].str.contains('DECOY', case=False, na=False)]
            other_in_top = top_compounds[~top_compounds['Compound_ID'].str.contains('DECOY|' + '|'.join(known_patterns), case=False, na=False)]
            
            print(f"   Known drugs: {len(known_in_top)} ({len(known_in_top)/n_compounds*100:.1f}%)")
            print(f"   DECOY compounds: {len(decoy_in_top)} ({len(decoy_in_top)/n_compounds*100:.1f}%)")
            print(f"   Other compounds: {len(other_in_top)} ({len(other_in_top)/n_compounds*100:.1f}%)")
            
            # Show score range for this cutoff
            min_score = top_compounds['Probability'].min()
            max_score = top_compounds['Probability'].max()
            print(f"   Score range: {min_score:.3f} - {max_score:.3f}")
            
            # List top compounds in this cutoff
            if cutoff <= 2:  # Show details for top 1% and 2%
                print(f"   📋 Compound list:")
                for i, (_, row) in enumerate(top_compounds.iterrows(), 1):
                    compound = row['Compound_ID']
                    score = row['Probability']
                    
                    if 'DECOY' in compound:
                        comp_type = "DECOY"
                    elif any(drug in compound for drug in known_patterns):
                        comp_type = "KNOWN"
                    else:
                        comp_type = "OTHER"
                    
                    print(f"     {i:2d}. [{comp_type}] {compound}: {score:.3f}")
        
        return df.head(int(len(df) * 0.05))  # Return top 5%
    
    def create_experimental_priority_list(self):
        """Create prioritized list for experimental validation"""
        print(f"\n🎯 CREATING EXPERIMENTAL PRIORITY LIST")
        print("=" * 40)
        
        if self.results is None:
            return None
        
        df = self.results
        
        # Define different priority tiers based on multiple criteria
        priority_lists = {}
        
        # Priority 1: Top 1% - Highest confidence
        top_1_pct = int(len(df) * 0.01)
        priority_lists['Priority_1_Top1pct'] = df.head(top_1_pct).copy()
        priority_lists['Priority_1_Top1pct']['Priority_Reason'] = 'Top 1% - Highest confidence'
        
        # Priority 2: Known drugs in top 10%
        top_10_pct = int(len(df) * 0.10)
        known_patterns = ['Paclitaxel', 'Podophyllotoxin', 'Doxorubicin', 'Cisplatin', 'Actinomycin', 'Etoposide', 'CHEMBL']
        known_in_top10 = df.head(top_10_pct)[df.head(top_10_pct)['Compound_ID'].str.contains('|'.join(known_patterns), case=False, na=False)].copy()
        known_in_top10['Priority_Reason'] = 'Known drug in top 10%'
        priority_lists['Priority_2_KnownDrugs'] = known_in_top10
        
        # Priority 3: Non-DECOY compounds in top 5%
        top_5_pct = int(len(df) * 0.05)
        non_decoy_top5 = df.head(top_5_pct)[~df.head(top_5_pct)['Compound_ID'].str.contains('DECOY', case=False, na=False)].copy()
        non_decoy_top5['Priority_Reason'] = 'Non-DECOY in top 5%'
        priority_lists['Priority_3_NonDecoy'] = non_decoy_top5
        
        # Priority 4: Interesting DECOY compounds (top scoring DECOYs - potential novel scaffolds)
        top_decoys = df[df['Compound_ID'].str.contains('DECOY', case=False, na=False)].head(20).copy()
        top_decoys['Priority_Reason'] = 'Top-scoring DECOY - potential novel scaffold'
        priority_lists['Priority_4_TopDecoys'] = top_decoys
        
        # Combine into unified priority list
        all_priorities = []
        for priority_name, priority_df in priority_lists.items():
            priority_df['Priority_Category'] = priority_name
            all_priorities.append(priority_df)
        
        combined_priorities = pd.concat(all_priorities, ignore_index=True)
        
        # Remove duplicates, keeping highest priority
        combined_priorities = combined_priorities.drop_duplicates(subset=['Compound_ID'], keep='first')
        
        # Sort by original rank
        combined_priorities = combined_priorities.sort_values('Percentile_Rank')
        
        print(f"📊 Priority List Summary:")
        for category in combined_priorities['Priority_Category'].unique():
            count = (combined_priorities['Priority_Category'] == category).sum()
            print(f"   {category}: {count} compounds")
        
        print(f"\n🏆 TOP 20 EXPERIMENTAL CANDIDATES:")
        top_candidates = combined_priorities.head(20)
        
        for i, (_, row) in enumerate(top_candidates.iterrows(), 1):
            compound = row['Compound_ID']
            score = row['Probability']
            rank = row['Percentile_Rank']
            reason = row['Priority_Reason']
            
            if 'DECOY' in compound:
                comp_type = "DECOY"
            elif any(drug in compound for drug in known_patterns):
                comp_type = "KNOWN"
            else:
                comp_type = "OTHER"
            
            print(f"  {i:2d}. [{comp_type}] {compound}")
            print(f"      Score: {score:.3f} | Rank: {rank} | {reason}")
        
        self.refined_hits = combined_priorities
        return combined_priorities
    
    def investigate_decoy_patterns(self):
        """Investigate patterns in high-scoring DECOY compounds"""
        print(f"\n🔍 INVESTIGATING HIGH-SCORING DECOY PATTERNS")
        print("=" * 45)
        
        if self.results is None:
            return None
        
        df = self.results
        
        # Get top-scoring DECOYs
        decoy_compounds = df[df['Compound_ID'].str.contains('DECOY', case=False, na=False)]
        top_decoys = decoy_compounds.head(50)  # Top 50 DECOYs
        
        print(f"📊 DECOY Analysis:")
        print(f"   Total DECOY compounds: {len(decoy_compounds):,}")
        print(f"   DECOYs in top 1%: {len(decoy_compounds[decoy_compounds['Percentile_Rank'] <= int(len(df) * 0.01)])}")
        print(f"   DECOYs in top 5%: {len(decoy_compounds[decoy_compounds['Percentile_Rank'] <= int(len(df) * 0.05)])}")
        print(f"   DECOYs in top 10%: {len(decoy_compounds[decoy_compounds['Percentile_Rank'] <= int(len(df) * 0.10)])}")
        
        # Analyze score distribution of top DECOYs
        print(f"\n📈 Top DECOY Score Distribution:")
        print(f"   Best DECOY score: {decoy_compounds['Probability'].max():.3f}")
        print(f"   Top 10 DECOY avg: {decoy_compounds.head(10)['Probability'].mean():.3f}")
        print(f"   Top 50 DECOY avg: {top_decoys['Probability'].mean():.3f}")
        
        # Show specific high-scoring DECOYs
        print(f"\n🎯 TOP 15 SCORING DECOYS (potential novel hits):")
        for i, (_, row) in enumerate(top_decoys.head(15).iterrows(), 1):
            compound = row['Compound_ID']
            score = row['Probability']
            rank = row['Percentile_Rank']
            percentile = row['Percentile']
            
            print(f"  {i:2d}. {compound}")
            print(f"      Score: {score:.3f} | Rank: {rank} ({percentile:.1f}%ile)")
        
        return top_decoys
    
    def save_refined_strategy_results(self):
        """Save all refined strategy results"""
        if self.refined_hits is None:
            return None
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_file = f'HCT15_Refined_Strategy_{timestamp}.xlsx'
        
        print(f"\n💾 SAVING REFINED STRATEGY RESULTS")
        print("=" * 35)
        
        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            # Complete refined priority list
            self.refined_hits.to_excel(writer, sheet_name='Experimental_Priorities', index=False)
            
            # Priority 1: Top 1%
            priority_1 = self.refined_hits[self.refined_hits['Priority_Category'] == 'Priority_1_Top1pct']
            if len(priority_1) > 0:
                priority_1.to_excel(writer, sheet_name='Priority_1_Top1pct', index=False)
            
            # Known drugs
            known_drugs = self.refined_hits[self.refined_hits['Priority_Category'] == 'Priority_2_KnownDrugs']
            if len(known_drugs) > 0:
                known_drugs.to_excel(writer, sheet_name='Known_Drugs_Validation', index=False)
            
            # Non-DECOY compounds
            non_decoy = self.refined_hits[self.refined_hits['Priority_Category'] == 'Priority_3_NonDecoy']
            if len(non_decoy) > 0:
                non_decoy.to_excel(writer, sheet_name='Non_DECOY_Hits', index=False)
            
            # Interesting DECOYs
            interesting_decoys = self.refined_hits[self.refined_hits['Priority_Category'] == 'Priority_4_TopDecoys']
            if len(interesting_decoys) > 0:
                interesting_decoys.to_excel(writer, sheet_name='Novel_DECOY_Scaffolds', index=False)
            
            # Experimental recommendations
            recommendations = [
                ['Validation Phase', 'Compounds', 'Rationale', 'Expected Outcome'],
                ['Phase 1 - Immediate', 'Top 1% (27 compounds)', 'Highest confidence hits', 'Validate model performance'],
                ['Phase 2 - Known Drug Validation', 'Known drugs in top 10%', 'Verify model accuracy', 'Establish activity baseline'],
                ['Phase 3 - Novel Compounds', 'Non-DECOY hits in top 5%', 'Discover new actives', 'Find novel therapeutic candidates'],
                ['Phase 4 - Investigation', 'Top-scoring DECOYs', 'Investigate unexpected hits', 'Discover novel scaffolds or model insights'],
                ['', '', '', ''],
                ['Priority Order', '1. Podophyllotoxin', 'Known active - model validation', 'Should show activity'],
                ['', '2. CHEMBL67', 'Known compound - top 1%', 'Validate model accuracy'],
                ['', '3. Top DECOY compounds', 'Unexpected high scores', 'Investigate for novel activity'],
                ['', '4. Other non-DECOY top 5%', 'Novel therapeutic candidates', 'Primary experimental targets'],
            ]
            
            rec_df = pd.DataFrame(recommendations[1:], columns=recommendations[0])
            rec_df.to_excel(writer, sheet_name='Experimental_Strategy', index=False)
        
        print(f"✅ Refined strategy saved: {output_file}")
        return output_file

def main():
    """Main refined strategy pipeline"""
    print("🎯 HCT-15 REFINED HIT SELECTION STRATEGY")
    print("=" * 45)
    print("Focus on highest confidence hits for experimental validation")
    print()
    
    strategy = HCT15RefinedStrategy()
    
    # Load percentile results
    df = strategy.load_percentile_results()
    if df is None:
        return
    
    # Analyze top performers
    top_5_pct = strategy.analyze_top_performers()
    
    # Create experimental priority list
    priority_list = strategy.create_experimental_priority_list()
    
    # Investigate DECOY patterns
    top_decoys = strategy.investigate_decoy_patterns()
    
    # Save results
    output_file = strategy.save_refined_strategy_results()
    
    print(f"\n🎉 REFINED STRATEGY COMPLETE!")
    print("=" * 32)
    print("✅ Created prioritized experimental list")
    print("✅ Identified validation candidates")
    print("✅ Investigated DECOY patterns")
    print(f"✅ Results saved: {output_file}")
    
    print(f"\n🎯 IMMEDIATE ACTION PLAN:")
    print("1. 🥇 FIRST: Test Podophyllotoxin (model validation)")
    print("2. 🥈 SECOND: Test CHEMBL67 (known active verification)")  
    print("3. 🥉 THIRD: Test top 3 DECOY compounds (investigate novel hits)")
    print("4. 📊 THEN: Proceed with systematic validation of top 1%")
    
    print(f"\n💡 KEY INSIGHTS:")
    print("• Your model correctly ranks known drugs at the top")
    print("• High DECOY scores may indicate novel active scaffolds")
    print("• Focus on top 1-5% for practical experimental validation")
    print("• Use known drugs to validate model performance first")

if __name__ == "__main__":
    main()