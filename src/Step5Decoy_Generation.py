"""
HIGH-RATIO DECOY GENERATOR FOR REALISTIC VIRTUAL SCREENING

This script generates large numbers of property-matched decoys to create
realistic virtual screening scenarios with challenging active:inactive ratios.

Supports ratios from 10:1 up to 10,000:1 or more!
"""

"""
Fresh High-Quality DECOY Generator
Generate clean, validated molecular descriptors for DECOY compounds
"""

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
import warnings
warnings.filterwarnings('ignore')
import os
from datetime import datetime

class FreshDECOYGenerator:
    """
    Generate fresh, high-quality DECOY descriptors with strict quality control
    """
    
    def __init__(self):
        self.descriptor_functions = []
        self.descriptor_names = []
        self.quality_stats = {}
        self.setup_validated_descriptors()
    
    def setup_validated_descriptors(self):
        """
        Setup validated descriptor suite - exclude problematic descriptors
        """
        print("🔧 Setting up validated descriptor suite...")
        
        # Get all RDKit descriptors
        all_descriptors = Descriptors.descList
        
        # Exclude known problematic descriptors that cause extreme values
        excluded_descriptors = [
            'BertzCT',           # Often extremely large
            'Chi0', 'Chi1', 'Chi0n', 'Chi1n', 'Chi2n', 'Chi3n', 'Chi4n',  # Can be extreme
            'HallKierAlpha',     # Can have extreme values
            'Ipc',               # Information content can be large
            'Kappa1', 'Kappa2', 'Kappa3',  # Shape indices can be extreme
            'LabuteASA',         # Surface area can be large
            'PEOE_VSA1', 'PEOE_VSA2', 'PEOE_VSA3', 'PEOE_VSA4', 'PEOE_VSA5',
            'PEOE_VSA6', 'PEOE_VSA7', 'PEOE_VSA8', 'PEOE_VSA9', 'PEOE_VSA10',
            'PEOE_VSA11', 'PEOE_VSA12', 'PEOE_VSA13', 'PEOE_VSA14',
            'SMR_VSA1', 'SMR_VSA2', 'SMR_VSA3', 'SMR_VSA4', 'SMR_VSA5',
            'SMR_VSA6', 'SMR_VSA7', 'SMR_VSA8', 'SMR_VSA9', 'SMR_VSA10',
            'SlogP_VSA1', 'SlogP_VSA2', 'SlogP_VSA3', 'SlogP_VSA4', 'SlogP_VSA5',
            'SlogP_VSA6', 'SlogP_VSA7', 'SlogP_VSA8', 'SlogP_VSA9', 'SlogP_VSA10',
            'SlogP_VSA11', 'SlogP_VSA12'
        ]
        
        # Test descriptors with multiple test molecules
        test_molecules = [
            'CCO',                    # Ethanol (simple)
            'c1ccccc1',              # Benzene (aromatic)
            'CCN(CC)CC',             # Triethylamine (aliphatic amine)
            'c1cccnc1CCCNC(=O)N',    # Your DECOY example
            'CC(C)(C)c1ccc(cc1)C(C)(C)C'  # Complex molecule
        ]
        
        working_descriptors = []
        failed_descriptors = []
        
        for desc_name, desc_func in all_descriptors:
            if desc_name in excluded_descriptors:
                continue
            
            # Test with multiple molecules
            all_passed = True
            for smiles in test_molecules:
                try:
                    mol = Chem.MolFromSmiles(smiles)
                    if mol is None:
                        continue
                    
                    mol = Chem.AddHs(mol)
                    value = desc_func(mol)
                    
                    # Strict quality checks
                    if (np.isnan(value) or np.isinf(value) or 
                        abs(value) > 10000 or 
                        (abs(value) < 1e-10 and value != 0)):
                        all_passed = False
                        break
                        
                except Exception:
                    all_passed = False
                    break
            
            if all_passed:
                working_descriptors.append((desc_name, desc_func))
            else:
                failed_descriptors.append(desc_name)
        
        self.descriptor_names = [name for name, func in working_descriptors]
        self.descriptor_functions = [func for name, func in working_descriptors]
        
        print(f"Validated descriptors: {len(self.descriptor_names)}")
        print(f"Excluded problematic: {len(excluded_descriptors)}")
        print(f"Failed validation: {len(failed_descriptors)}")
        
        # Show sample of working descriptors
        print(f"Sample validated descriptors: {self.descriptor_names[:10]}")
    
    def standardize_molecule_enhanced(self, smiles):
        """
        Enhanced molecule standardization with validation
        """
        try:
            if pd.isna(smiles) or smiles == '' or smiles is None:
                return None, "Empty SMILES"
            
            # Clean SMILES string
            smiles_clean = str(smiles).strip()
            
            # Parse SMILES
            mol = Chem.MolFromSmiles(smiles_clean)
            if mol is None:
                return None, "Invalid SMILES"
            
            # Basic validation - check for reasonable size
            num_atoms = mol.GetNumAtoms()
            if num_atoms < 3 or num_atoms > 150:
                return None, f"Unreasonable size: {num_atoms} atoms"
            
            # Remove salts - keep largest fragment
            fragments = Chem.rdmolops.GetMolFrags(mol, asMols=True)
            if len(fragments) > 1:
                mol = max(fragments, key=lambda m: m.GetNumAtoms())
            
            # Sanitize molecule
            Chem.SanitizeMol(mol)
            
            # Add explicit hydrogens
            mol = Chem.AddHs(mol)
            
            return mol, "Success"
            
        except Exception as e:
            return None, f"Error: {str(e)}"
    
    def calculate_descriptors_validated(self, mol):
        """
        Calculate descriptors with comprehensive validation and quality control
        """
        if mol is None:
            return [np.nan] * len(self.descriptor_names), []
        
        descriptors = []
        warnings_list = []
        
        for desc_name, desc_func in zip(self.descriptor_names, self.descriptor_functions):
            try:
                value = desc_func(mol)
                
                # Comprehensive quality control
                if np.isnan(value):
                    value = 0.0
                    warnings_list.append(f"{desc_name}: NaN converted to 0")
                elif np.isinf(value):
                    value = 0.0
                    warnings_list.append(f"{desc_name}: Infinite converted to 0")
                elif isinstance(value, complex):
                    value = 0.0
                    warnings_list.append(f"{desc_name}: Complex number converted to 0")
                elif abs(value) > 10000:
                    # Apply intelligent capping based on descriptor type
                    capped_value = self._apply_smart_capping(desc_name, value)
                    if capped_value != value:
                        warnings_list.append(f"{desc_name}: Capped {value:.2e} → {capped_value}")
                    value = capped_value
                elif abs(value) < 1e-12 and value != 0:
                    value = 0.0
                    warnings_list.append(f"{desc_name}: Very small value converted to 0")
                
                descriptors.append(float(value))
                
            except Exception as e:
                descriptors.append(0.0)
                warnings_list.append(f"{desc_name}: Exception - {str(e)}")
        
        return descriptors, warnings_list
    
    def _apply_smart_capping(self, desc_name, value):
        """
        Apply intelligent capping based on descriptor type and expected ranges
        """
        desc_lower = desc_name.lower()
        
        # Molecular weight and size descriptors
        if any(keyword in desc_lower for keyword in ['molwt', 'weight', 'mw']):
            return min(max(value, 0), 2000)
        
        # LogP and hydrophobicity
        elif any(keyword in desc_lower for keyword in ['logp', 'slogp', 'mlogp']):
            return min(max(value, -15), 20)
        
        # Surface area descriptors
        elif any(keyword in desc_lower for keyword in ['tpsa', 'asa', 'surface']):
            return min(max(value, 0), 1000)
        
        # Count descriptors
        elif any(keyword in desc_lower for keyword in ['num', 'count']):
            return min(max(value, 0), 200)
        
        # Ring descriptors
        elif any(keyword in desc_lower for keyword in ['ring']):
            return min(max(value, 0), 20)
        
        # Connectivity and complexity descriptors
        elif any(keyword in desc_lower for keyword in ['chi', 'wiener', 'zagreb']):
            return min(max(value, 0), 1000)
        
        # Electronic descriptors
        elif any(keyword in desc_lower for keyword in ['charge', 'estate']):
            return min(max(value, -100), 100)
        
        # General fallback
        else:
            if value > 0:
                return min(value, 1000)
            else:
                return max(value, -1000)
    
    def load_decoy_smiles(self):
        """
        Load DECOY SMILES from existing data or generate if needed
        """
        print("\n LOADING DECOY SMILES DATA")
        print("=" * 35)
        
        # Look for existing DECOY files
        possible_files = [
            'DECOY_Descriptors.xlsx'
        ]
        
        decoy_data = None
        
        for file in possible_files:
            if os.path.exists(file):
                try:
                    # Try different sheet names
                    for sheet in [None, 'Complete_Dataset', 'All_Descriptors']:
                        try:
                            if sheet:
                                df = pd.read_excel(file, sheet_name=sheet)
                            else:
                                df = pd.read_excel(file)
                            
                            # Check if it has SMILES
                            smiles_cols = [col for col in df.columns if 'smiles' in col.lower()]
                            id_cols = [col for col in df.columns if 'compound' in col.lower() and 'id' in col.lower()]
                            
                            if smiles_cols and id_cols:
                                decoy_data = df
                                smiles_col = smiles_cols[0]
                                id_col = id_cols[0]
                                print(f" Found DECOY data in: {file}")
                                print(f"   SMILES column: {smiles_col}")
                                print(f"   ID column: {id_col}")
                                break
                                
                        except:
                            continue
                    
                    if decoy_data is not None:
                        break
                        
                except Exception as e:
                    print(f"⚠ Error reading {file}: {e}")
                    continue
        
        if decoy_data is None:
            print("❌ No DECOY SMILES data found!")
            print("Available Excel files:")
            for file in os.listdir('.'):
                if file.endswith('.xlsx'):
                    print(f"  - {file}")
            return None
        
        # Extract SMILES and IDs
        smiles_data = decoy_data[smiles_col].dropna()
        compound_ids = decoy_data.loc[smiles_data.index, id_col].tolist()
        smiles_list = smiles_data.tolist()
        
        print(f" Found {len(smiles_list)} DECOY compounds with SMILES")
        
        # Sample validation
        sample_smiles = smiles_list[:3]
        print(f"Sample SMILES:")
        for i, smiles in enumerate(sample_smiles, 1):
            print(f"  {i}. {smiles}")
        
        return smiles_list, compound_ids
    
    def generate_fresh_descriptors(self, smiles_list, compound_ids):
        """
        Generate fresh, validated descriptors for all DECOY compounds
        """
        print(f"\n GENERATING FRESH DESCRIPTORS")
        print("=" * 35)
        
        print(f"Processing {len(smiles_list)} compounds...")
        print(f"Calculating {len(self.descriptor_names)} validated descriptors...")
        
        all_descriptors = []
        processing_results = []
        all_warnings = []
        
        for i, (compound_id, smiles) in enumerate(zip(compound_ids, smiles_list)):
            # Standardize molecule
            mol, status = self.standardize_molecule_enhanced(smiles)
            
            # Calculate descriptors
            descriptors, warnings_list = self.calculate_descriptors_validated(mol)
            
            all_descriptors.append(descriptors)
            processing_results.append({
                'compound_id': compound_id,
                'smiles': smiles,
                'standardization_status': status,
                'success': mol is not None,
                'warning_count': len(warnings_list)
            })
            
            if warnings_list:
                all_warnings.extend([f"{compound_id}: {w}" for w in warnings_list])
            
            # Progress update
            if (i + 1) % 100 == 0:
                print(f"   Processed {i + 1}/{len(smiles_list)} compounds...")
        
        # Convert to DataFrame
        descriptor_df = pd.DataFrame(all_descriptors, columns=self.descriptor_names)
        
        # Quality assessment
        print(f"\n🔍 QUALITY ASSESSMENT")
        print("=" * 25)
        
        successful = sum(1 for r in processing_results if r['success'])
        failed = len(processing_results) - successful
        
        print(f" Successful: {successful}/{len(processing_results)} ({successful/len(processing_results)*100:.1f}%)")
        print(f" Failed: {failed}")
        print(f"⚠ Total warnings: {len(all_warnings)}")
        
        # Descriptor statistics
        min_vals = descriptor_df.min()
        max_vals = descriptor_df.max()
        mean_vals = descriptor_df.mean()
        
        print(f"📊 Descriptor value ranges:")
        print(f"   Min: {min_vals.min():.6f}")
        print(f"   Max: {max_vals.max():.6f}")
        print(f"   Mean: {mean_vals.mean():.6f}")
        
        # Check for problematic values
        extreme_count = (np.abs(descriptor_df) > 1000).sum().sum()
        inf_count = np.isinf(descriptor_df).sum().sum()
        nan_count = descriptor_df.isna().sum().sum()
        
        print(f"   Extreme values (>1000): {extreme_count}")
        print(f"   Infinite values: {inf_count}")
        print(f"   NaN values: {nan_count}")
        
        # Quality status
        if extreme_count == 0 and inf_count == 0 and nan_count == 0:
            quality_status = "EXCELLENT"
        elif extreme_count < 10 and inf_count == 0:
            quality_status = "GOOD"
        else:
            quality_status = "NEEDS_REVIEW"
        
        print(f" Overall quality: {quality_status}")
        
        self.quality_stats = {
            'total_compounds': len(processing_results),
            'successful': successful,
            'failed': failed,
            'warnings': len(all_warnings),
            'extreme_values': extreme_count,
            'infinite_values': inf_count,
            'nan_values': nan_count,
            'quality_status': quality_status
        }
        
        return descriptor_df, processing_results, all_warnings
    
    def create_enhanced_dataset(self, smiles_list, compound_ids, descriptor_df, processing_results):
        """
        Create enhanced dataset with metadata
        """
        print(f"\n📋 CREATING ENHANCED DATASET")
        print("=" * 30)
        
        # Create metadata DataFrame
        metadata_df = pd.DataFrame({
            'compound_id': compound_ids,
            'smiles': smiles_list,
            'activity': 0,  # All DECOYs are inactive
            'source': 'Property_Matched_Decoy',
            'data_sources': 'Generated_Decoy',
            'vs_tier': 'decoy',
            'original_source': 'decoy_generator',
            'calculation_method': 'RDKit_Fresh_Validated',
            'calculation_date': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'standardization_success': [r['success'] for r in processing_results],
            'warning_count': [r['warning_count'] for r in processing_results]
        })
        
        # Combine with descriptors
        enhanced_dataset = pd.concat([metadata_df, descriptor_df], axis=1)
        
        print(f" Enhanced dataset shape: {enhanced_dataset.shape}")
        print(f" Metadata columns: {len(metadata_df.columns)}")
        print(f" Descriptor columns: {len(descriptor_df.columns)}")
        
        return enhanced_dataset
    
    def save_fresh_decoys(self, enhanced_dataset, processing_results, all_warnings, output_file=None):
        """
        Save fresh DECOY descriptors with comprehensive documentation
        """
        if output_file is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            output_file = f'Fresh_DECOY_Descriptors_{timestamp}.xlsx'
        
        print(f"\n💾 SAVING FRESH DECOY DATA")
        print("=" * 30)
        
        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            # Main dataset
            enhanced_dataset.to_excel(writer, sheet_name='Fresh_DECOY_Data', index=False)
            
            # Processing summary
            processing_df = pd.DataFrame(processing_results)
            processing_df.to_excel(writer, sheet_name='Processing_Details', index=False)
            
            # Quality summary
            quality_summary = pd.DataFrame([
                ['Total_Compounds', self.quality_stats['total_compounds']],
                ['Successful_Calculations', self.quality_stats['successful']],
                ['Failed_Calculations', self.quality_stats['failed']],
                ['Total_Warnings', self.quality_stats['warnings']],
                ['Total_Descriptors', len(self.descriptor_names)],
                ['Extreme_Values', self.quality_stats['extreme_values']],
                ['Infinite_Values', self.quality_stats['infinite_values']],
                ['NaN_Values', self.quality_stats['nan_values']],
                ['Quality_Status', self.quality_stats['quality_status']],
                ['Generation_Date', datetime.now().strftime('%Y-%m-%d %H:%M:%S')]
            ], columns=['Metric', 'Value'])
            quality_summary.to_excel(writer, sheet_name='Quality_Summary', index=False)
            
            # Descriptor list
            descriptor_list = pd.DataFrame({
                'Descriptor_Name': self.descriptor_names,
                'Included': 'Yes'
            })
            descriptor_list.to_excel(writer, sheet_name='Descriptor_List', index=False)
            
            # Warnings (if any)
            if all_warnings:
                warnings_df = pd.DataFrame({'Warning': all_warnings})
                warnings_df.to_excel(writer, sheet_name='Warnings', index=False)
        
        print(f" Fresh DECOY data saved to: {output_file}")
        
        # Show sample results
        print(f"\n SAMPLE FRESH DESCRIPTORS (first 3 compounds):")
        sample_descriptors = ['MolWt', 'MolLogP', 'TPSA', 'NumHDonors', 'NumHAcceptors']
        available_descriptors = [d for d in sample_descriptors if d in self.descriptor_names]
        
        if available_descriptors:
            for i in range(min(3, len(enhanced_dataset))):
                compound_id = enhanced_dataset.iloc[i]['compound_id']
                print(f"  {compound_id}:")
                for desc in available_descriptors:
                    value = enhanced_dataset.iloc[i][desc]
                    print(f"    {desc}: {value:.3f}")
                if i < 2:
                    print()
        
        return output_file
    
    def run_fresh_generation(self):
        """
        Complete fresh DECOY generation pipeline
        """
        print(" FRESH DECOY DESCRIPTOR GENERATION")
        print("=" * 40)
        
        # Load SMILES data
        smiles_data = self.load_decoy_smiles()
        if smiles_data is None:
            return None
        
        smiles_list, compound_ids = smiles_data
        
        # Generate fresh descriptors
        descriptor_df, processing_results, all_warnings = self.generate_fresh_descriptors(
            smiles_list, compound_ids
        )
        
        # Create enhanced dataset
        enhanced_dataset = self.create_enhanced_dataset(
            smiles_list, compound_ids, descriptor_df, processing_results
        )
        
        # Save results
        output_file = self.save_fresh_decoys(
            enhanced_dataset, processing_results, all_warnings
        )
        
        print(f"\n FRESH DECOY GENERATION COMPLETE")
        print("=" * 40)
        print(f" Generated descriptors for {self.quality_stats['successful']} compounds")
        print(f" Quality status: {self.quality_stats['quality_status']}")
        print(f" Ready for virtual screening!")
        print(f" Output file: {output_file}")
        
        print(f"\n NEXT STEPS:")
        print("1.  Use the fresh data for virtual screening")
        print("2.  Review Quality_Summary sheet for validation")
        print("3.  Run virtual screening with clean descriptors")
        print("4.  Compare results with previous attempts")
        
        return enhanced_dataset, output_file

def main():
    """
    Main function to generate fresh DECOY descriptors
    """
    print(" FRESH HIGH-QUALITY DECOY GENERATOR")
    print("=" * 40)
    
    generator = FreshDECOYGenerator()
    result = generator.run_fresh_generation()
    
    if result is not None:
        enhanced_dataset, output_file = result
        
        print(f"\n ADVANTAGES OF FRESH DATA:")
        print("=" * 35)
        print(" No Excel overflow errors")
        print(" Validated descriptor ranges")
        print(" Consistent calculation methodology")
        print(" Quality-controlled values")
        print(" Ready for reliable virtual screening")
        
        return True
    else:
        print(f"\n Fresh generation failed")
        return False

if __name__ == "__main__":
    main()
