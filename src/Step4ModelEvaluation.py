# -*- coding: utf-8 -*-
"""
Created on Mon Aug 18 12:17:04 2025

@author: ASUS TUF F15
"""

import pandas as pd
import numpy as np
import pickle
import glob
import os

def evaluate_all_models():
    """
    Evaluate all trained models and identify the best one for deployment
    """
    print("=== STEP 3A MODEL EVALUATION ===")
    print(f"Working directory: {os.getcwd()}")
    
    # Find all result files
    result_files = glob.glob("Step3A_*_BCR.txt")
    
    if not result_files:
        print("No result files found! Make sure you're in the correct directory.")
        return None
    
    print(f"\nFound {len(result_files)} result files:")
    
    all_results = {}
    
    # Evaluate individual models
    fitness_functions = ['RF', 'DTREE', 'KNN', 'BNB', 'GNB']
    
    for ff in fitness_functions:
        filename = f"Step3A_GA_Results_{ff}_BCR.txt"
        if os.path.exists(filename):
            try:
                df = pd.read_csv(filename, sep='\t')
                best_idx = df['EXT_BCR'].idxmax()
                best_performance = df.iloc[best_idx]
                
                all_results[f"Individual_{ff}"] = {
                    'file': filename,
                    'model_type': 'Individual',
                    'algorithm': ff,
                    'ext_bcr': best_performance['EXT_BCR'],
                    'ext_acc': best_performance['EXT_ACC'],
                    'ext_f1': best_performance['EXT_F1'],
                    'ext_sensitivity': best_performance['EXT_SE'],
                    'ext_specificity': best_performance['EXT_SP'],
                    'generalization_gap': best_performance['GA_BCR'] - best_performance['EXT_BCR'],
                    'best_model_index': best_idx
                }
                
                print(f"  ✓ {ff}: External BCR = {best_performance['EXT_BCR']:.3f}")
                
            except Exception as e:
                print(f"  ✗ {ff}: Error reading file - {e}")
    
    # Evaluate ensemble model
    ensemble_file = "Step3A_Ensembles_GA_BCR.txt"
    if os.path.exists(ensemble_file):
        try:
            df = pd.read_csv(ensemble_file, sep='\t')
            best_idx = df['EXT_BCR'].idxmax()
            best_performance = df.iloc[best_idx]
            
            all_results['Ensemble'] = {
                'file': ensemble_file,
                'model_type': 'Ensemble',
                'algorithm': 'Multiple',
                'ext_bcr': best_performance['EXT_BCR'],
                'ext_acc': best_performance['EXT_ACC'],
                'ext_f1': best_performance['EXT_F1'],
                'ext_sensitivity': best_performance['EXT_SE'],
                'ext_specificity': best_performance['EXT_SP'],
                'generalization_gap': best_performance['GA_BCR'] - best_performance['EXT_BCR'],
                'best_model_index': best_idx
            }
            
            print(f"  ✓ Ensemble: External BCR = {best_performance['EXT_BCR']:.3f}")
            
        except Exception as e:
            print(f"  ✗ Ensemble: Error reading file - {e}")
    
    if not all_results:
        print("No valid results found!")
        return None
    
    # Create summary DataFrame
    summary_df = pd.DataFrame(all_results).T
    summary_df = summary_df.sort_values('ext_bcr', ascending=False)
    
    print(f"\n=== MODEL PERFORMANCE RANKING ===")
    print(f"{'Rank':<4} {'Model':<15} {'Algorithm':<8} {'Ext_BCR':<8} {'Ext_ACC':<8} {'Ext_F1':<8} {'Gen_Gap':<8}")
    print("-" * 70)
    
    for i, (model_name, row) in enumerate(summary_df.iterrows(), 1):
        print(f"{i:<4} {model_name:<15} {row['algorithm']:<8} {row['ext_bcr']:<8.3f} "
              f"{row['ext_acc']:<8.3f} {row['ext_f1']:<8.3f} {row['generalization_gap']:<8.3f}")
    
    # Recommend best model
    best_model = summary_df.index[0]
    best_row = summary_df.iloc[0]
    
    print(f"\n=== RECOMMENDED MODEL FOR DEPLOYMENT ===")
    print(f"Model: {best_model}")
    print(f"Type: {best_row['model_type']}")
    print(f"Algorithm: {best_row['algorithm']}")
    print(f"External BCR: {best_row['ext_bcr']:.3f}")
    print(f"External Accuracy: {best_row['ext_acc']:.3f}")
    print(f"External F1: {best_row['ext_f1']:.3f}")
    print(f"Sensitivity: {best_row['ext_sensitivity']:.3f}")
    print(f"Specificity: {best_row['ext_specificity']:.3f}")
    print(f"Generalization Gap: {best_row['generalization_gap']:.3f}")
    
    # Performance assessment
    print(f"\n=== PERFORMANCE ASSESSMENT ===")
    bcr = best_row['ext_bcr']
    gen_gap = abs(best_row['generalization_gap'])
    
    if bcr >= 0.8:
        perf_level = "Excellent"
    elif bcr >= 0.7:
        perf_level = "Good"
    elif bcr >= 0.6:
        perf_level = "Moderate"
    else:
        perf_level = "Poor"
    
    if gen_gap <= 0.05:
        gen_level = "Excellent generalization"
    elif gen_gap <= 0.1:
        gen_level = "Good generalization"
    elif gen_gap <= 0.15:
        gen_level = "Moderate generalization"
    else:
        gen_level = "Poor generalization"
    
    print(f"Performance Level: {perf_level} (BCR = {bcr:.3f})")
    print(f"Generalization: {gen_level} (Gap = {gen_gap:.3f})")
    
    # Deployment recommendation
    if bcr >= 0.7 and gen_gap <= 0.15:
        print("✅ RECOMMENDED FOR DEPLOYMENT")
        print("This model is suitable for virtual screening.")
    elif bcr >= 0.6:
        print("⚠️  CONDITIONAL DEPLOYMENT")
        print("Model may be useful but consider additional validation.")
    else:
        print("❌ NOT RECOMMENDED FOR DEPLOYMENT")
        print("Model performance is too low for reliable virtual screening.")
    
    return summary_df, best_model, best_row

def check_model_files():
    """
    Check which model files are available for deployment
    """
    print(f"\n=== AVAILABLE MODEL FILES FOR DEPLOYMENT ===")
    
    # Check for pickle files containing trained models
    pkl_files = glob.glob("Step3A_*.pkl")
    
    print(f"Found {len(pkl_files)} model files:")
    for pkl_file in pkl_files:
        file_size = os.path.getsize(pkl_file) / (1024*1024)  # Size in MB
        print(f"  ✓ {pkl_file} ({file_size:.1f} MB)")
    
    # Check for application domain file
    ad_file = "Step3A_ApplicationDomain_Conditions.pkl"
    if os.path.exists(ad_file):
        print(f"  ✓ Application Domain file: {ad_file}")
    else:
        print(f"  ✗ Application Domain file missing: {ad_file}")
    
    # Check for normalization file
    norm_files = glob.glob("*Step2_Scaler.pkl")
    if norm_files:
        print(f"  ✓ Normalization file: {norm_files[0]}")
    else:
        print(f"  ✗ Normalization file missing")
    
    return pkl_files

if __name__ == "__main__":
    # Run evaluation
    try:
        summary_df, best_model, best_row = evaluate_all_models()
        
        if summary_df is not None:
            # Check available model files
            pkl_files = check_model_files()
            
            # Save summary
            summary_df.to_csv("Model_Performance_Summary.csv")
            print(f"\n📊 Performance summary saved to: Model_Performance_Summary.csv")
            
        else:
            print("Could not complete evaluation.")
            
    except Exception as e:
        print(f"Error during evaluation: {e}")
        print("Please make sure you're in the correct directory with Step 3A results.")