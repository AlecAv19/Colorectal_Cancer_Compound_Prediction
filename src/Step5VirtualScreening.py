# -*- coding: utf-8 -*-
"""
Created on Mon Aug 18 16:37:48 2025

@author: ASUS TUF F15
"""
# -*- coding: utf-8 -*-
"""
FIXED Virtual Screening Pipeline with Proper PCA Handling
Fixes dimension mismatch issues in applicability domain checking
"""

import pandas as pd
import numpy as np
import pickle
import os
import warnings
warnings.filterwarnings('ignore')

class FixedVirtualScreeningPipeline:
    """
    Virtual Screening Pipeline with proper preprocessing chain and fixed PCA handling
    """
    
    def __init__(self, model_type='ensemble'):
        self.model_type = model_type
        self.models = None
        self.scaler = None
        self.feature_selector = None
        self.feature_indices = None
        self.application_domain = None
        self.is_loaded = False
        self.preprocessing_pipeline = {}
        
    def load_best_model(self):
        """Load model and understand the complete preprocessing pipeline"""
        print("=== LOADING MODEL WITH COMPLETE PREPROCESSING PIPELINE ===")
        
        # Load model files
        if self.model_type == 'ensemble':
            model_file = "Step3A_Ensembles_GA_BCR.pkl"
            results_file = "Step3A_Ensembles_GA_BCR.txt"
        else:
            model_file = f"Step3A_GA_Results_{self.model_type}_BCR.pkl"
            results_file = f"Step3A_GA_Results_{self.model_type}_BCR.txt"
        
        # Load model
        with open(model_file, 'rb') as f:
            model_data = pickle.load(f)
        
        # Extract model components
        if self.model_type == 'ensemble':
            initial_params, data_modeling, ga_data, models_best, models_variables = model_data
            results_df = pd.read_csv(results_file, sep='\t')
            best_idx = results_df['EXT_BCR'].idxmax()
            population = ga_data[2]
            best_ensemble = population[best_idx]
            selected_indices = np.where(best_ensemble == 1)[0]
            self.models = [models_best[i] for i in selected_indices]
            self.feature_sets = [models_variables[i] for i in selected_indices]
            print(f"✓ Loaded ensemble with {len(self.models)} models")
        else:
            initial_params, data_modeling, ga_data = model_data
            results_df = pd.read_csv(results_file, sep='\t')
            best_idx = results_df['EXT_BCR'].idxmax()
            population = ga_data[2]
            others = ga_data[1]
            best_features = population[best_idx]
            self.models = [others[best_idx, 5]]
            self.feature_sets = [np.where(best_features == 1)[0]]
            print(f"✓ Loaded {self.model_type} model")
        
        # Analyze preprocessing pipeline
        print(f"\n=== ANALYZING PREPROCESSING PIPELINE ===")
        
        if len(data_modeling) > 6:
            self.scaler = data_modeling[0]
            self.feature_indices = data_modeling[1] 
            self.feature_selector = data_modeling[2]
            
            # Store pipeline info
            self.preprocessing_pipeline['step1_input'] = self.scaler.n_features_in_
            print(f"Step 1 - Input: {self.preprocessing_pipeline['step1_input']} features")
            
            # Determine features after scaling (same number)
            self.preprocessing_pipeline['step2_after_scaling'] = self.preprocessing_pipeline['step1_input']
            print(f"Step 2 - After scaling: {self.preprocessing_pipeline['step2_after_scaling']} features")
            
            # Determine features after selection
            if hasattr(self.feature_selector, 'transform'):
                # Test the feature selector to see output size
                dummy_input = np.zeros((1, self.preprocessing_pipeline['step2_after_scaling']))
                try:
                    dummy_output = self.feature_selector.transform(dummy_input)
                    self.preprocessing_pipeline['step3_after_selection'] = dummy_output.shape[1]
                except:
                    # Fallback to feature indices if available
                    if isinstance(self.feature_indices, (list, np.ndarray)):
                        self.preprocessing_pipeline['step3_after_selection'] = len(self.feature_indices)
                    else:
                        self.preprocessing_pipeline['step3_after_selection'] = None
            elif isinstance(self.feature_indices, (list, np.ndarray)):
                self.preprocessing_pipeline['step3_after_selection'] = len(self.feature_indices)
            else:
                self.preprocessing_pipeline['step3_after_selection'] = None
            
            if self.preprocessing_pipeline['step3_after_selection']:
                print(f"Step 3 - After feature selection: {self.preprocessing_pipeline['step3_after_selection']} features")
            else:
                print(f"Step 3 - Feature selection: Could not determine output size")
        
        # Load and analyze application domain with detailed debugging
        ad_file = "Step3A_ApplicationDomain_Conditions.pkl"
        if os.path.exists(ad_file):
            with open(ad_file, 'rb') as f:
                ad_data = pickle.load(f)
            
            print(f"\n=== APPLICATION DOMAIN ANALYSIS ===")
            print(f"Application domain data structure: {len(ad_data)} components")
            
            self.application_domain = {
                'pca': ad_data[0],
                'max_dist': ad_data[1],
                'n_components': ad_data[2],
                'centroid': ad_data[3]
            }
            
            # Detailed PCA analysis
            pca = self.application_domain['pca']
            centroid = self.application_domain['centroid']
            n_components = self.application_domain['n_components']
            
            print(f"PCA input features: {pca.n_features_in_}")
            print(f"PCA components: {pca.n_components_}")
            print(f"Stored n_components: {n_components}")
            print(f"Centroid shape: {centroid.shape}")
            print(f"Max distance: {self.application_domain['max_dist']:.4f}")
            
            # Store PCA requirements
            self.preprocessing_pipeline['step4_pca_input'] = pca.n_features_in_
            self.preprocessing_pipeline['step4_pca_output'] = pca.n_components_
            self.preprocessing_pipeline['step4_centroid_dims'] = centroid.shape[1] if centroid.ndim > 1 else len(centroid)
            
            print(f"Step 4 - PCA expects: {self.preprocessing_pipeline['step4_pca_input']} features")
            print(f"Step 4 - PCA outputs: {self.preprocessing_pipeline['step4_pca_output']} components")
            print(f"Step 4 - Centroid has: {self.preprocessing_pipeline['step4_centroid_dims']} dimensions")
            
            # Check for dimension consistency issues
            if self.preprocessing_pipeline['step4_pca_output'] != self.preprocessing_pipeline['step4_centroid_dims']:
                print(f"⚠ DIMENSION MISMATCH DETECTED!")
                print(f"   PCA will output: {self.preprocessing_pipeline['step4_pca_output']} components")
                print(f"   But centroid has: {self.preprocessing_pipeline['step4_centroid_dims']} dimensions")
                print(f"   This will cause broadcasting errors!")
                
                # Try to fix by using only the available centroid dimensions
                if self.preprocessing_pipeline['step4_centroid_dims'] <= self.preprocessing_pipeline['step4_pca_output']:
                    print(f"💡 Will use only first {self.preprocessing_pipeline['step4_centroid_dims']} PCA components")
                else:
                    print(f"❌ Cannot fix: centroid has more dimensions than PCA output")
            
            # Check if pipeline is consistent
            if (self.preprocessing_pipeline.get('step3_after_selection') and 
                self.preprocessing_pipeline['step3_after_selection'] == self.preprocessing_pipeline['step4_pca_input']):
                print(f"✅ Pipeline is consistent: Feature selection → PCA dimensions match!")
            else:
                print(f"⚠ Pipeline inconsistency detected:")
                print(f"   Feature selection output: {self.preprocessing_pipeline.get('step3_after_selection')}")
                print(f"   PCA input requirement: {self.preprocessing_pipeline['step4_pca_input']}")
                
        else:
            print("⚠ Application domain not found")
        
        self.is_loaded = True
        
        # Show model performance
        best_performance = results_df.iloc[best_idx]
        print(f"\n=== MODEL PERFORMANCE ===")
        print(f"External BCR: {best_performance['EXT_BCR']:.3f}")
        print(f"External Accuracy: {best_performance['EXT_ACC']:.3f}")
        
        return True
    
    def preprocess_data_with_pipeline(self, X):
        """Apply the complete preprocessing pipeline step by step"""
        print(f"\n🔄 APPLYING COMPLETE PREPROCESSING PIPELINE")
        print(f"Input: {X.shape}")
        
        X_processed = X.copy()
        
        # Step 1: Check input dimensions
        expected_input = self.preprocessing_pipeline.get('step1_input')
        if expected_input and X_processed.shape[1] != expected_input:
            print(f"❌ Input dimension mismatch!")
            print(f"   Your data: {X_processed.shape[1]} features")
            print(f"   Expected: {expected_input} features") 
            raise ValueError(f"Input must have {expected_input} features")
        
        # Step 2: Apply scaling/normalization
        if self.scaler is not None:
            X_processed = self.scaler.transform(X_processed)
            print(f"Step 2 - After scaling: {X_processed.shape}")
        
        # Step 3: Apply feature selection
        if self.feature_selector is not None:
            try:
                X_processed = self.feature_selector.transform(X_processed)
                print(f"Step 3 - After feature selection: {X_processed.shape}")
            except Exception as e:
                print(f"⚠ Feature selector failed: {e}")
                # Try feature indices as backup
                if isinstance(self.feature_indices, (list, np.ndarray)):
                    X_processed = X_processed[:, self.feature_indices]
                    print(f"Step 3 - Using feature indices: {X_processed.shape}")
                else:
                    raise e
        elif isinstance(self.feature_indices, (list, np.ndarray)):
            X_processed = X_processed[:, self.feature_indices]
            print(f"Step 3 - Applied feature indices: {X_processed.shape}")
        
        # Step 4: Verify dimensions for PCA
        if self.application_domain:
            expected_pca_input = self.preprocessing_pipeline.get('step4_pca_input')
            if expected_pca_input and X_processed.shape[1] != expected_pca_input:
                print(f"❌ PCA dimension mismatch after preprocessing!")
                print(f"   After feature selection: {X_processed.shape[1]} features")
                print(f"   PCA expects: {expected_pca_input} features")
                raise ValueError(f"Feature selection should produce {expected_pca_input} features for PCA")
            else:
                print(f"✅ Dimensions ready for PCA: {X_processed.shape}")
        
        return X_processed
    
    def check_applicability_domain(self, X):
        """Check applicability domain with proper dimension handling"""
        if self.application_domain is None:
            return np.ones(len(X), dtype=bool), np.zeros(len(X))
        
        print(f"🔄 Checking applicability domain...")
        
        try:
            pca = self.application_domain['pca']
            max_dist = self.application_domain['max_dist']
            centroid = self.application_domain['centroid']
            
            # X should already be properly preprocessed at this point
            print(f"   PCA input: {X.shape[1]} features (expected: {pca.n_features_in_})")
            
            # Transform to PCA space
            X_pca = pca.transform(X)
            print(f"   PCA output: {X_pca.shape}")
            print(f"   Centroid shape: {centroid.shape}")
            
            # Handle dimension mismatch between PCA output and centroid
            pca_components = X_pca.shape[1]
            centroid_dims = centroid.shape[1] if centroid.ndim > 1 else len(centroid)
            
            if pca_components != centroid_dims:
                print(f"⚠ Fixing dimension mismatch:")
                print(f"   PCA components: {pca_components}")
                print(f"   Centroid dimensions: {centroid_dims}")
                
                if centroid_dims < pca_components:
                    # Use only the first N components to match centroid
                    X_pca = X_pca[:, :centroid_dims]
                    print(f"   Using first {centroid_dims} PCA components")
                elif centroid_dims > pca_components:
                    # Pad centroid or truncate it to match PCA
                    if centroid.ndim > 1:
                        centroid = centroid[:, :pca_components]
                    else:
                        centroid = centroid[:pca_components]
                    print(f"   Truncated centroid to {pca_components} dimensions")
            
            # Ensure centroid is 1D for broadcasting
            if centroid.ndim > 1:
                centroid = centroid.flatten()
            
            print(f"   Final PCA shape: {X_pca.shape}")
            print(f"   Final centroid shape: {centroid.shape}")
            
            # Calculate distances to centroid
            distances = np.sqrt(np.sum((X_pca - centroid) ** 2, axis=1))
            
            # Check which compounds are within domain
            within_domain = distances <= max_dist
            
            print(f"✅ Applicability domain check complete")
            print(f"   Compounds within domain: {within_domain.sum()}/{len(within_domain)} ({within_domain.mean()*100:.1f}%)")
            print(f"   Distance range: {distances.min():.3f} - {distances.max():.3f}")
            print(f"   Max allowed distance: {max_dist:.3f}")
            
            return within_domain, distances
            
        except Exception as e:
            print(f"❌ Applicability domain check failed: {e}")
            import traceback
            traceback.print_exc()
            return np.zeros(len(X), dtype=bool), np.full(len(X), float('inf'))
    
    def predict_compound(self, X):
        """Make predictions with complete pipeline"""
        X = np.array(X)
        if X.ndim == 1:
            X = X.reshape(1, -1)
        
        print(f"🤖 Making predictions for {X.shape[0]} compounds with {X.shape[1]} features")
        
        try:
            # Apply complete preprocessing pipeline
            X_processed = self.preprocess_data_with_pipeline(X)
            
            # Check applicability domain (after preprocessing)
            within_domain, distances = self.check_applicability_domain(X_processed)
            
            # Make predictions using the processed features
            if self.model_type == 'ensemble':
                all_probs = []
                print(f"🔮 Making ensemble predictions...")
                
                for i, (model, features) in enumerate(zip(self.models, self.feature_sets)):
                    try:
                        X_model = X_processed[:, features]
                        probs = model.predict_proba(X_model)[:, 1]
                        all_probs.append(probs)
                        print(f"   Model {i+1}: {probs.mean():.3f} avg probability")
                    except Exception as e:
                        print(f"   Model {i+1} failed: {e}")
                        all_probs.append(np.full(len(X_processed), 0.5))
                
                if all_probs:
                    mean_probs = np.mean(all_probs, axis=0)
                    predictions = (mean_probs > 0.5).astype(int)
                else:
                    mean_probs = np.full(len(X_processed), 0.5)
                    predictions = np.zeros(len(X_processed))
            else:
                model = self.models[0]
                features = self.feature_sets[0]
                X_model = X_processed[:, features]
                mean_probs = model.predict_proba(X_model)[:, 1]
                predictions = model.predict(X_model)
            
            print(f"✅ Predictions completed successfully")
            print(f"   Average probability: {mean_probs.mean():.3f}")
            print(f"   Predicted active: {predictions.sum()}/{len(predictions)} ({predictions.mean()*100:.1f}%)")
            
            return predictions, mean_probs, within_domain, distances
            
        except Exception as e:
            print(f"❌ Prediction pipeline failed: {e}")
            import traceback
            traceback.print_exc()
            # Return safe defaults
            n_compounds = len(X)
            return (np.zeros(n_compounds), 
                   np.full(n_compounds, 0.5), 
                   np.zeros(n_compounds, dtype=bool), 
                   np.full(n_compounds, float('inf')))
    
    def screen_compounds(self, compound_data, compound_ids=None, threshold=0.5):
        """Screen compounds with complete pipeline"""
        print("=== VIRTUAL SCREENING WITH COMPLETE PIPELINE ===")
        
        if isinstance(compound_data, pd.DataFrame):
            if compound_ids is None:
                compound_ids = compound_data.index.tolist()
            X = compound_data.values
        else:
            X = np.array(compound_data)
            if compound_ids is None:
                compound_ids = [f"Compound_{i+1}" for i in range(len(X))]
        
        print(f"Screening {len(X)} compounds with threshold {threshold}...")
        
        # Make predictions
        predictions, probabilities, within_domain, distances = self.predict_compound(X)
        
        # Create results
        results = pd.DataFrame({
            'Compound_ID': compound_ids,
            'Prediction': predictions,
            'Probability': probabilities,
            'Within_Domain': within_domain,
            'Domain_Distance': distances,
            'Active': probabilities >= threshold,
            'Reliable_Hit': (probabilities >= threshold) & within_domain
        })
        
        # Summary
        n_total = len(results)
        n_active = results['Active'].sum()
        n_within_domain = results['Within_Domain'].sum()
        n_reliable_hits = results['Reliable_Hit'].sum()
        
        print(f"\n=== SCREENING RESULTS ===")
        print(f"Total compounds: {n_total}")
        print(f"Predicted active: {n_active} ({n_active/n_total*100:.1f}%)")
        print(f"Within domain: {n_within_domain} ({n_within_domain/n_total*100:.1f}%)")
        print(f"Reliable hits: {n_reliable_hits} ({n_reliable_hits/n_total*100:.1f}%)")
        
        return results.sort_values('Probability', ascending=False)
    
    def save_results(self, results, filename="fixed_vs_results.xlsx"):
        """Save results with domain information"""
        with pd.ExcelWriter(filename, engine='openpyxl') as writer:
            # All results
            results.to_excel(writer, sheet_name='All_Results', index=False)
            
            # Predicted active
            active_results = results[results['Active']]
            active_results.to_excel(writer, sheet_name='Predicted_Active', index=False)
            
            # Reliable hits (active + within domain)
            reliable_hits = results[results['Reliable_Hit']]
            reliable_hits.to_excel(writer, sheet_name='Reliable_Hits', index=False)
            
            # Within domain only
            within_domain = results[results['Within_Domain']]
            within_domain.to_excel(writer, sheet_name='Within_Domain', index=False)
            
            # Summary
            summary = pd.DataFrame({
                'Metric': [
                    'Total_Compounds', 'Predicted_Active', 'Within_Domain', 
                    'Reliable_Hits', 'Hit_Rate_%', 'Domain_Coverage_%', 'Reliable_Hit_Rate_%'
                ],
                'Value': [
                    len(results), results['Active'].sum(), results['Within_Domain'].sum(),
                    results['Reliable_Hit'].sum(), results['Active'].mean() * 100,
                    results['Within_Domain'].mean() * 100, results['Reliable_Hit'].mean() * 100
                ]
            })
            summary.to_excel(writer, sheet_name='Summary', index=False)
        
        print(f"📊 Results saved to: {filename}")

def run_fixed_screening():
    """Run screening with fixed preprocessing pipeline"""
    print("🔧 FIXED VIRTUAL SCREENING WITH COMPLETE PIPELINE")
    print("=" * 55)
    
    # Initialize
    vs = FixedVirtualScreeningPipeline(model_type='ensemble')
    
    try:
        # Load model and analyze pipeline
        vs.load_best_model()
        
        # Load your 1024-feature data
        data_file = 'VSPipeline_1024Features_Input.xlsx'
        
        if os.path.exists(data_file):
            print(f"\n=== LOADING 1024-FEATURE DATASET ===")
            compounds_df = pd.read_excel(data_file, sheet_name='Complete_Dataset')
            
            # Separate metadata and features
            metadata_cols = ['compound_id', 'smiles', 'activity', 'source', 'confidence', 
                           'data_sources', 'vs_tier', 'original_source']
            descriptor_cols = [col for col in compounds_df.columns if col not in metadata_cols]
            
            compound_metadata = compounds_df[metadata_cols]
            compound_descriptors = compounds_df[descriptor_cols]
            
            # Clean data
            compound_descriptors = compound_descriptors.fillna(compound_descriptors.median())
            compound_descriptors = compound_descriptors.replace([np.inf, -np.inf], np.nan)
            compound_descriptors = compound_descriptors.fillna(0)
            
            compound_ids = compound_metadata['compound_id'].tolist()
            
            print(f"✅ Loaded {len(compounds_df)} compounds with {len(descriptor_cols)} features")
            
            # Verify dimensions match pipeline requirements
            expected_input = vs.preprocessing_pipeline.get('step1_input')
            if expected_input and len(descriptor_cols) == expected_input:
                print(f"✅ Perfect match! Data has {len(descriptor_cols)} features, pipeline expects {expected_input}")
            else:
                print(f"⚠ Dimension check: Data has {len(descriptor_cols)} features, pipeline expects {expected_input}")
                
        else:
            print(f"❌ File not found: {data_file}")
            print("Please run the 1024-feature calculator first!")
            return None
        
        # Run screening with a more reasonable threshold
        print(f"\n=== TESTING DIFFERENT THRESHOLDS ===")
        
        # First, get all results
        results_50 = vs.screen_compounds(compound_descriptors, compound_ids, threshold=0.5)
        
        # Show distribution of probabilities
        print(f"\n=== PROBABILITY DISTRIBUTION ===")
        print(f"Min probability: {results_50['Probability'].min():.3f}")
        print(f"Max probability: {results_50['Probability'].max():.3f}")
        print(f"Mean probability: {results_50['Probability'].mean():.3f}")
        print(f"Median probability: {results_50['Probability'].median():.3f}")
        
        # Show results at different thresholds
        for thresh in [0.3, 0.4, 0.5, 0.6, 0.7]:
            n_active = (results_50['Probability'] >= thresh).sum()
            n_reliable = ((results_50['Probability'] >= thresh) & results_50['Within_Domain']).sum()
            print(f"Threshold {thresh}: {n_active} active, {n_reliable} reliable hits")
        
        # Use a more reasonable threshold
        final_threshold = 0.5
        results = vs.screen_compounds(compound_descriptors, compound_ids, threshold=final_threshold)
        
        # Save results
        vs.save_results(results, "fixed_vs_results.xlsx")
        
        # Show top hits
        print(f"\n=== TOP 10 HITS BY PROBABILITY ===")
        top_hits = results.head(10)
        
        for idx, row in top_hits.iterrows():
            domain_status = "✓" if row['Within_Domain'] else "⚠"
            active_status = "🎯" if row['Active'] else "○"
            print(f"{active_status} {domain_status} {row['Compound_ID']}: {row['Probability']:.3f} (Distance: {row['Domain_Distance']:.2f})")
        
        return results
        
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return None

if __name__ == "__main__":
    results = run_fixed_screening()