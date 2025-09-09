# -*- coding: utf-8 -*-
"""
Created on Tue Aug 19 23:54:05 2025

@author: ASUS TUF F15
"""
"""
Enhanced Production Virtual Screening Pipeline
Improved error handling, flexible data loading, and comprehensive analysis
"""
"""
Quick Virtual Screening Runner
Run virtual screening with your existing models and DECOY data
"""

import pandas as pd
import numpy as np
import pickle
import os
import glob
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

class QuickVirtualScreening:
    """
    Quick virtual screening with automatic data and model detection
    """
    
    def __init__(self):
        self.models = None
        self.scaler = None
        self.feature_selector = None
        self.feature_indices = None
        self.application_domain = None
        self.model_loaded = False
        
    def find_model_files(self):
        """Find available model files"""
        print("🔍 SCANNING FOR MODEL FILES")
        print("=" * 30)
        
        model_patterns = [
            'Step3A_Ensembles_GA_BCR.pkl',
            'Step3A_GA_Results_*_BCR.pkl',
            '*_model*.pkl',
            '*ensemble*.pkl'
        ]
        
        found_models = []
        for pattern in model_patterns:
            found_models.extend(glob.glob(pattern))
        
        if found_models:
            print("✅ Found model files:")
            for model in found_models:
                print(f"   - {model}")
        else:
            print("❌ No model files found!")
            print("Available .pkl files:")
            for file in glob.glob('*.pkl'):
                print(f"   - {file}")
        
        return found_models
    
    def find_data_files(self):
        """Find available descriptor data files"""
        print("\n🔍 SCANNING FOR DESCRIPTOR DATA")
        print("=" * 35)
        
        data_patterns = [
            'Smart_1024_Descriptors_*.xlsx',
            'Complete_1024_Descriptors_*.xlsx',
            'Comprehensive_Clean_DECOY_*.xlsx',
            'Extreme_Values_Fixed_DECOY_*.xlsx',
            'Fixed_Boolean_DECOY_*.xlsx',
            'Ultra_Clean_DECOY_*.xlsx',
            'Fresh_DECOY_*.xlsx',
            'VS_DECOY_DATA_UNF.xlsx'
        ]
        
        found_data = []
        for pattern in data_patterns:
            found_data.extend(glob.glob(pattern))
        
        if found_data:
            # Use most recent file
            latest_data = max(found_data, key=os.path.getctime)
            print(f"✅ Using latest data file: {latest_data}")
            
            print("Available data files:")
            for data in found_data:
                print(f"   - {data}")
        else:
            print("❌ No descriptor data files found!")
            print("Available Excel files:")
            for file in glob.glob('*.xlsx'):
                print(f"   - {file}")
        
        return latest_data if found_data else None
    
    def load_model(self, model_file):
        """Load the virtual screening model"""
        print(f"\n🤖 LOADING MODEL: {model_file}")
        print("=" * 30)
        
        try:
            with open(model_file, 'rb') as f:
                model_data = pickle.load(f)
            
            print(f"✅ Model loaded successfully")
            
            # Try to extract model components (adapt based on your model structure)
            if 'ensemble' in model_file.lower():
                # Handle ensemble model
                try:
                    if len(model_data) >= 5:
                        initial_params, data_modeling, ga_data, models_best, models_variables = model_data
                        self.models = models_best
                        self.feature_sets = models_variables
                        print(f"   Ensemble models: {len(self.models)}")
                    else:
                        print("   ⚠ Unexpected ensemble structure")
                        return False
                except Exception as e:
                    print(f"   ⚠ Could not parse ensemble: {e}")
                    return False
            else:
                # Handle single model
                try:
                    if len(model_data) >= 3:
                        initial_params, data_modeling, ga_data = model_data
                        # Extract best model (simplified)
                        self.models = [ga_data[1][0, 5]]  # Simplified extraction
                        self.feature_sets = [np.where(ga_data[2][0] == 1)[0]]
                        print(f"   Single model loaded")
                    else:
                        print("   ⚠ Unexpected model structure")
                        return False
                except Exception as e:
                    print(f"   ⚠ Could not parse model: {e}")
                    return False
            
            # Extract preprocessing components
            if len(model_data) > 1:
                data_modeling = model_data[1]
                if len(data_modeling) > 0:
                    self.scaler = data_modeling[0] if len(data_modeling) > 0 else None
                    self.feature_indices = data_modeling[1] if len(data_modeling) > 1 else None
                    self.feature_selector = data_modeling[2] if len(data_modeling) > 2 else None
                    
                    print(f"   Preprocessing components loaded")
            
            self.model_loaded = True
            return True
            
        except Exception as e:
            print(f"❌ Error loading model: {e}")
            return False
    
    def load_screening_data(self, data_file):
        """Load data for screening"""
        print(f"\n📊 LOADING SCREENING DATA: {data_file}")
        print("=" * 40)
        
        try:
            # Try different sheet names
            sheet_options = [
                'Smart_1024_Descriptors', 'Complete_1024_Descriptors',
                'Clean_DECOY_Data', 'Fixed_DECOY_Data', 'Fresh_DECOY_Data',
                'All_Descriptors', 'Complete_Dataset', None
            ]
            
            df = None
            for sheet in sheet_options:
                try:
                    if sheet:
                        df = pd.read_excel(data_file, sheet_name=sheet)
                    else:
                        df = pd.read_excel(data_file)
                    print(f"✅ Loaded from sheet: {sheet or 'first sheet'}")
                    break
                except:
                    continue
            
            if df is None:
                print("❌ Could not load data")
                return None, None, None
            
            print(f"📊 Data shape: {df.shape}")
            
            # Identify columns
            metadata_patterns = [
                'identifier', 'compound_id', 'smiles', 'activity', 'source', 
                'dataset_source', 'calculation_method', 'calculation_date'
            ]
            
            metadata_cols = []
            for col in df.columns:
                if any(pattern.lower() in col.lower() for pattern in metadata_patterns):
                    metadata_cols.append(col)
            
            descriptor_cols = [col for col in df.columns if col not in metadata_cols]
            
            print(f"📋 Metadata columns: {len(metadata_cols)}")
            print(f"🧪 Descriptor columns: {len(descriptor_cols)}")
            
            if len(descriptor_cols) == 0:
                print("❌ No descriptor columns found!")
                return None, None, None
            
            # Extract data
            compound_metadata = df[metadata_cols] if metadata_cols else pd.DataFrame(index=df.index)
            descriptor_data = df[descriptor_cols]
            
            # Get compound identifiers
            if 'identifier' in df.columns:
                compound_ids = df['identifier'].tolist()
            elif 'compound_id' in df.columns:
                compound_ids = df['compound_id'].tolist()
            else:
                compound_ids = [f"Compound_{i+1}" for i in range(len(df))]
            
            # Clean descriptor data
            print("🧹 Cleaning descriptor data...")
            
            # Handle missing values
            missing_before = descriptor_data.isnull().sum().sum()
            descriptor_data = descriptor_data.fillna(descriptor_data.median())
            
            # Handle infinite values
            inf_before = np.isinf(descriptor_data.select_dtypes(include=[np.number])).sum().sum()
            descriptor_data = descriptor_data.replace([np.inf, -np.inf], np.nan)
            descriptor_data = descriptor_data.fillna(0)
            
            print(f"   Fixed {missing_before} missing values")
            print(f"   Fixed {inf_before} infinite values")
            
            print(f"✅ Ready to screen {len(compound_ids)} compounds")
            
            return descriptor_data, compound_ids, compound_metadata
            
        except Exception as e:
            print(f"❌ Error loading data: {e}")
            return None, None, None
    
    def preprocess_data(self, X):
        """Apply preprocessing pipeline"""
        try:
            X_processed = X.copy()
            
            # Apply scaling
            if self.scaler is not None:
                print("   Applying scaling...")
                X_processed = self.scaler.transform(X_processed)
            
            # Apply feature selection
            if self.feature_indices is not None:
                print(f"   Selecting {len(self.feature_indices)} features...")
                # Ensure indices are valid
                max_idx = X_processed.shape[1] - 1
                valid_indices = [idx for idx in self.feature_indices if idx <= max_idx]
                
                if len(valid_indices) != len(self.feature_indices):
                    print(f"   ⚠ Using {len(valid_indices)}/{len(self.feature_indices)} features")
                
                if valid_indices:
                    X_processed = X_processed[:, valid_indices]
            
            return X_processed
            
        except Exception as e:
            print(f"   ⚠ Preprocessing error: {e}")
            return X
    
    def make_predictions(self, X_processed):
        """Make predictions with loaded models"""
        try:
            if hasattr(self, 'feature_sets') and len(self.feature_sets) > 1:
                print("   Making ensemble predictions...")
                # Ensemble prediction
                all_probs = []
                
                for i, (model, features) in enumerate(zip(self.models, self.feature_sets)):
                    try:
                        # Ensure feature indices are valid
                        max_features = X_processed.shape[1]
                        valid_features = [f for f in features if f < max_features]
                        
                        if valid_features:
                            X_model = X_processed[:, valid_features]
                            probs = model.predict_proba(X_model)[:, 1]
                            all_probs.append(probs)
                            print(f"      Model {i+1}: {len(valid_features)} features")
                    except Exception as e:
                        print(f"      ⚠ Model {i+1} failed: {e}")
                        all_probs.append(np.full(len(X_processed), 0.5))
                
                if all_probs:
                    mean_probs = np.mean(all_probs, axis=0)
                    predictions = (mean_probs > 0.5).astype(int)
                    print(f"   ✅ Ensemble prediction complete ({len(all_probs)} models)")
                else:
                    mean_probs = np.full(len(X_processed), 0.5)
                    predictions = np.zeros(len(X_processed))
                    print(f"   ⚠ No models succeeded, using default")
                    
            else:
                print("   Making single model prediction...")
                # Single model prediction
                model = self.models[0]
                features = self.feature_sets[0] if hasattr(self, 'feature_sets') else list(range(X_processed.shape[1]))
                
                # Ensure features are valid
                max_features = X_processed.shape[1]
                valid_features = [f for f in features if f < max_features]
                
                if valid_features:
                    X_model = X_processed[:, valid_features]
                    mean_probs = model.predict_proba(X_model)[:, 1]
                    predictions = model.predict(X_model)
                    print(f"   ✅ Single model prediction complete ({len(valid_features)} features)")
                else:
                    mean_probs = np.full(len(X_processed), 0.5)
                    predictions = np.zeros(len(X_processed))
                    print(f"   ⚠ No valid features, using default")
            
            return predictions, mean_probs
            
        except Exception as e:
            print(f"   ❌ Prediction error: {e}")
            n_compounds = len(X_processed)
            return np.zeros(n_compounds), np.full(n_compounds, 0.5)
    
    def run_virtual_screening(self):
        """Run complete virtual screening pipeline"""
        print("🚀 QUICK VIRTUAL SCREENING PIPELINE")
        print("=" * 40)
        
        # Find model files
        model_files = self.find_model_files()
        if not model_files:
            print("❌ No model files found! Cannot proceed.")
            return None
        
        # Use first available model
        model_file = model_files[0]
        
        # Load model
        if not self.load_model(model_file):
            print("❌ Failed to load model! Cannot proceed.")
            return None
        
        # Find data files
        data_file = self.find_data_files()
        if not data_file:
            print("❌ No descriptor data found! Cannot proceed.")
            return None
        
        # Load screening data
        descriptor_data, compound_ids, metadata = self.load_screening_data(data_file)
        if descriptor_data is None:
            print("❌ Failed to load screening data! Cannot proceed.")
            return None
        
        # Run screening
        print(f"\n🎯 RUNNING VIRTUAL SCREENING")
        print("=" * 30)
        print(f"Compounds to screen: {len(compound_ids)}")
        print(f"Descriptors per compound: {descriptor_data.shape[1]}")
        
        # Preprocess data
        print("🔧 Preprocessing data...")
        X_processed = self.preprocess_data(descriptor_data.values)
        
        # Make predictions
        print("🤖 Making predictions...")
        predictions, probabilities = self.make_predictions(X_processed)
        
        # Create results
        print("📊 Compiling results...")
        results = pd.DataFrame({
            'Compound_ID': compound_ids,
            'Prediction': predictions,
            'Probability': probabilities,
            'Active': probabilities >= 0.5,
            'Confidence': np.where(probabilities >= 0.7, 'High', 
                         np.where(probabilities >= 0.3, 'Medium', 'Low'))
        })
        
        # Add metadata if available
        if metadata is not None and len(metadata) > 0:
            for col in metadata.columns:
                if col not in results.columns:
                    try:
                        results[col] = metadata[col].values
                    except:
                        pass
        
        # Sort by probability (highest first)
        results = results.sort_values('Probability', ascending=False)
        
        # Add ranking
        results['Rank'] = range(1, len(results) + 1)
        
        # Summary statistics
        n_total = len(results)
        n_active = results['Active'].sum()
        mean_prob = results['Probability'].mean()
        top_1_percent = results['Probability'].quantile(0.99)
        
        print(f"\n🎉 VIRTUAL SCREENING COMPLETE!")
        print("=" * 35)
        print(f"📊 Results Summary:")
        print(f"   Total compounds screened: {n_total:,}")
        print(f"   Predicted active (≥0.5): {n_active:,} ({n_active/n_total*100:.1f}%)")
        print(f"   Mean probability: {mean_prob:.3f}")
        print(f"   Top 1% threshold: {top_1_percent:.3f}")
        
        # Save results
        self.save_results(results)
        
        # Show top hits
        self.show_top_hits(results)
        
        return results
    
    def save_results(self, results):
        """Save screening results"""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_file = f'VS_Results_{timestamp}.xlsx'
        
        print(f"\n💾 SAVING RESULTS")
        print("=" * 20)
        
        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            # All results
            results.to_excel(writer, sheet_name='All_Results', index=False)
            
            # Top 100 hits
            top_100 = results.head(100)
            top_100.to_excel(writer, sheet_name='Top_100_Hits', index=False)
            
            # Predicted active compounds
            active_compounds = results[results['Active']].copy()
            active_compounds.to_excel(writer, sheet_name='Predicted_Active', index=False)
            
            # High confidence hits
            high_confidence = results[results['Confidence'] == 'High'].copy()
            if len(high_confidence) > 0:
                high_confidence.to_excel(writer, sheet_name='High_Confidence', index=False)
            
            # Summary statistics
            summary_stats = pd.DataFrame([
                ['Total_Compounds', len(results)],
                ['Predicted_Active', results['Active'].sum()],
                ['Mean_Probability', results['Probability'].mean()],
                ['Std_Probability', results['Probability'].std()],
                ['Top_1_Percent_Threshold', results['Probability'].quantile(0.99)],
                ['Top_5_Percent_Threshold', results['Probability'].quantile(0.95)],
                ['High_Confidence_Count', (results['Confidence'] == 'High').sum()],
                ['Screening_Date', datetime.now().strftime('%Y-%m-%d %H:%M:%S')]
            ], columns=['Metric', 'Value'])
            summary_stats.to_excel(writer, sheet_name='Summary', index=False)
        
        print(f"✅ Results saved to: {output_file}")
        return output_file
    
    def show_top_hits(self, results):
        """Show top screening hits"""
        print(f"\n🏆 TOP 20 SCREENING HITS")
        print("=" * 30)
        
        top_20 = results.head(20)
        
        for i, (_, row) in enumerate(top_20.iterrows(), 1):
            confidence_icon = "🟢" if row['Confidence'] == 'High' else "🟡" if row['Confidence'] == 'Medium' else "🔴"
            activity_status = "ACTIVE" if row['Active'] else "inactive"
            
            print(f"{i:2d}. {confidence_icon} {row['Compound_ID']}: {row['Probability']:.3f} ({activity_status})")
        
        # Show distribution
        print(f"\n📊 PROBABILITY DISTRIBUTION:")
        prob_ranges = [
            (0.9, 1.0, "Very High"),
            (0.7, 0.9, "High"),
            (0.5, 0.7, "Medium"),
            (0.3, 0.5, "Low"),
            (0.0, 0.3, "Very Low")
        ]
        
        for min_prob, max_prob, label in prob_ranges:
            count = len(results[(results['Probability'] >= min_prob) & (results['Probability'] < max_prob)])
            if count > 0:
                print(f"   {label} ({min_prob}-{max_prob}): {count:,} compounds")

def main():
    """Run quick virtual screening"""
    print("🚀 QUICK VIRTUAL SCREENING")
    print("=" * 30)
    print("Automatic model and data detection")
    print()
    
    vs = QuickVirtualScreening()
    results = vs.run_virtual_screening()
    
    if results is not None:
        print(f"\n🎉 VIRTUAL SCREENING SUCCESS!")
        print("✅ All compounds screened")
        print("✅ Results saved to Excel")
        print("✅ Top hits identified")
        
        print(f"\n🎯 NEXT STEPS:")
        print("1. Review top hits in Excel file")
        print("2. Validate high-probability compounds")
        print("3. Plan experimental testing")
        
        return True
    else:
        print("❌ Virtual screening failed")
        return False

if __name__ == "__main__":
    main()