#!/usr/bin/env python3

"""
CENTAUR Machine Learning Pipeline (Simplified)
Implements Random Forest, Logistic Regression, and other ML models
for multi-omics biomarker discovery and cancer classification
"""

import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.metrics import roc_auc_score, precision_score, recall_score, f1_score, classification_report
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

class CENTAURMLPipeline:
    """CENTAUR Machine Learning Pipeline for multi-omics biomarker discovery"""
    
    def __init__(self, random_state=42):
        self.random_state = random_state
        self.models = {}
        self.scaler = StandardScaler()
        self.feature_importance = {}
        
    def load_data(self, feature_matrix_path):
        """Load multi-omics feature matrix"""
        self.data = pd.read_csv(feature_matrix_path, sep='\t')
        print(f"Loaded data with shape: {self.data.shape}")
        return self.data
    
    def prepare_features(self, target_column='sample_type'):
        """Prepare features and target for ML training"""
        # Separate features and target
        feature_columns = [col for col in self.data.columns 
                          if col not in ['sample_id', 'cancer_type', 'stage', 'sample_type']]
        
        self.X = self.data[feature_columns]
        self.y = self.data[target_column]
        
        # Convert target to binary (Cancer vs Control)
        self.y_binary = (self.y == 'Cancer').astype(int)
        
        print(f"Features: {len(feature_columns)}")
        print(f"Samples: {len(self.X)}")
        print(f"Cancer samples: {sum(self.y_binary)}")
        print(f"Control samples: {len(self.y_binary) - sum(self.y_binary)}")
        
        return self.X, self.y_binary
    
    def train_models(self):
        """Train multiple ML models"""
        
        # Initialize models
        self.models = {
            'RandomForest': RandomForestClassifier(
                n_estimators=100,
                max_depth=10,
                random_state=self.random_state
            ),
            'LogisticRegression': LogisticRegression(
                random_state=self.random_state,
                max_iter=1000
            ),
            'NaiveBayes': GaussianNB(),
            'SVM': SVC(
                probability=True,
                random_state=self.random_state
            )
        }
        
        # Scale features for models that need it
        X_scaled = self.scaler.fit_transform(self.X)
        
        # Train each model
        for name, model in self.models.items():
            print(f"\nTraining {name}...")
            
            if name in ['LogisticRegression', 'SVM']:
                model.fit(X_scaled, self.y_binary)
            else:
                model.fit(self.X, self.y_binary)
            
            print(f"{name} training completed")
    
    def evaluate_models(self):
        """Evaluate all trained models"""
        results = {}
        
        # Cross-validation setup
        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=self.random_state)
        
        for name, model in self.models.items():
            print(f"\nEvaluating {name}...")
            
            # Prepare data for evaluation
            if name in ['LogisticRegression', 'SVM']:
                X_eval = self.scaler.transform(self.X)
            else:
                X_eval = self.X
            
            # Cross-validation scores
            cv_scores = cross_val_score(model, X_eval, self.y_binary, cv=cv, scoring='roc_auc')
            
            # Predictions
            y_pred = model.predict(X_eval)
            y_pred_proba = model.predict_proba(X_eval)[:, 1]
            
            # Calculate metrics
            auc = roc_auc_score(self.y_binary, y_pred_proba)
            precision = precision_score(self.y_binary, y_pred)
            recall = recall_score(self.y_binary, y_pred)
            f1 = f1_score(self.y_binary, y_pred)
            
            results[name] = {
                'auc': auc,
                'precision': precision,
                'recall': recall,
                'f1_score': f1,
                'cv_auc_mean': cv_scores.mean(),
                'cv_auc_std': cv_scores.std(),
                'cv_scores': cv_scores
            }
            
            print(f"{name} - AUC: {auc:.3f}, Precision: {precision:.3f}, Recall: {recall:.3f}, F1: {f1:.3f}")
        
        self.evaluation_results = results
        return results
    
    def calculate_feature_importance(self):
        """Calculate feature importance for each model"""
        feature_names = self.X.columns.tolist()
        
        for name, model in self.models.items():
            if hasattr(model, 'feature_importances_'):
                # Tree-based models
                importance = model.feature_importances_
            elif hasattr(model, 'coef_'):
                # Linear models
                importance = np.abs(model.coef_[0])
            else:
                # Other models - use permutation importance
                from sklearn.inspection import permutation_importance
                if name in ['LogisticRegression', 'SVM']:
                    X_eval = self.scaler.transform(self.X)
                else:
                    X_eval = self.X
                
                perm_importance = permutation_importance(model, X_eval, self.y_binary, 
                                                       random_state=self.random_state)
                importance = perm_importance.importances_mean
            
            # Normalize importance scores
            importance = importance / np.sum(importance)
            
            # Create feature importance dataframe
            feature_importance_df = pd.DataFrame({
                'feature': feature_names,
                'importance': importance,
                'model': name
            }).sort_values('importance', ascending=False)
            
            feature_importance_df['rank'] = range(1, len(feature_importance_df) + 1)
            
            self.feature_importance[name] = feature_importance_df
            
            print(f"\n{name} - Top 5 Features:")
            print(feature_importance_df.head())
    
    def generate_predictions(self, new_data=None):
        """Generate predictions for new data"""
        if new_data is None:
            new_data = self.X
        
        predictions = {}
        
        for name, model in self.models.items():
            if name in ['LogisticRegression', 'SVM']:
                X_pred = self.scaler.transform(new_data)
            else:
                X_pred = new_data
            
            y_pred = model.predict(X_pred)
            y_pred_proba = model.predict_proba(X_pred)[:, 1]
            
            predictions[name] = {
                'predictions': y_pred,
                'probabilities': y_pred_proba
            }
        
        return predictions
    
    def save_results(self, output_dir='ml_results'):
        """Save all results to files"""
        import os
        os.makedirs(output_dir, exist_ok=True)
        
        # Save evaluation results
        eval_df = pd.DataFrame(self.evaluation_results).T
        eval_df.to_csv(f'{output_dir}/model_performance.tsv', sep='\t')
        
        # Save feature importance
        for name, importance_df in self.feature_importance.items():
            importance_df.to_csv(f'{output_dir}/{name}_feature_importance.tsv', sep='\t', index=False)
        
        # Save predictions
        predictions = self.generate_predictions()
        for name, pred_data in predictions.items():
            pred_df = pd.DataFrame({
                'sample_id': self.data['sample_id'],
                'true_label': self.y_binary,
                'predicted_label': pred_data['predictions'],
                'probability_cancer': pred_data['probabilities']
            })
            pred_df.to_csv(f'{output_dir}/{name}_predictions.tsv', sep='\t', index=False)
        
        print(f"Results saved to {output_dir}/")
    
    def create_summary_report(self, output_dir='ml_results'):
        """Create a summary report of ML results"""
        import os
        os.makedirs(output_dir, exist_ok=True)
        
        report = []
        report.append("# CENTAUR Machine Learning Results Summary")
        report.append("=" * 50)
        report.append("")
        
        # Model performance summary
        report.append("## Model Performance")
        report.append("")
        report.append("| Model | AUC | Precision | Recall | F1-Score | CV AUC (Mean ± Std) |")
        report.append("|-------|-----|-----------|--------|----------|---------------------|")
        
        for name, results in self.evaluation_results.items():
            report.append(f"| {name} | {results['auc']:.3f} | {results['precision']:.3f} | "
                         f"{results['recall']:.3f} | {results['f1_score']:.3f} | "
                         f"{results['cv_auc_mean']:.3f} ± {results['cv_auc_std']:.3f} |")
        
        report.append("")
        
        # Feature importance summary
        report.append("## Top Biomarkers by Model")
        report.append("")
        
        for name, importance_df in self.feature_importance.items():
            report.append(f"### {name}")
            report.append("")
            report.append("| Rank | Feature | Importance Score |")
            report.append("|------|---------|-----------------|")
            
            top_features = importance_df.head(5)
            for _, row in top_features.iterrows():
                report.append(f"| {row['rank']} | {row['feature']} | {row['importance']:.4f} |")
            
            report.append("")
        
        # Best performing model
        best_model = max(self.evaluation_results.items(), key=lambda x: x[1]['auc'])
        report.append("## Best Performing Model")
        report.append("")
        report.append(f"**{best_model[0]}** with AUC = {best_model[1]['auc']:.3f}")
        report.append("")
        
        # Clinical interpretation
        report.append("## Clinical Interpretation")
        report.append("")
        report.append("- **AUC > 0.8**: Good discriminative power for cancer detection")
        report.append("- **AUC > 0.9**: Excellent discriminative power for clinical use")
        report.append("- **Precision**: Proportion of predicted cancer cases that are actually cancer")
        report.append("- **Recall**: Proportion of actual cancer cases correctly identified")
        report.append("- **F1-Score**: Harmonic mean of precision and recall")
        report.append("")
        
        # Save report
        with open(f'{output_dir}/ml_results_summary.md', 'w') as f:
            f.write('\n'.join(report))
        
        print(f"Summary report saved to {output_dir}/ml_results_summary.md")

def main():
    """Main function to run the CENTAUR ML pipeline"""
    print("CENTAUR Machine Learning Pipeline")
    print("=" * 50)
    
    # Initialize pipeline
    pipeline = CENTAURMLPipeline()
    
    # For testing, create a sample feature matrix
    print("Creating sample feature matrix for testing...")
    
    # Generate sample data
    np.random.seed(42)
    n_samples = 100
    n_features = 15
    
    # Feature names based on our biomarker discovery
    feature_names = [
        'methylation_enhancer_CRXIP', 'methylation_SNP_upstream_PRK',
        'methylation_SNP_near_MYC', 'methylation_hypo_TMB_RNF',
        'methylation_hyper_TMB_PDL1', 'methylation_global_hypo_CpG',
        'fragmentomics_short_TFBS', 'fragmentomics_entropy',
        'fragmentomics_WPS_STAT1', 'fragmentomics_end_motif_CTCF',
        'variant_KRAS_mutation', 'variant_TMB', 'variant_DNA_burden',
        'expression_immune_activation', 'expression_t_cell_exhaustion'
    ]
    
    # Generate sample data
    data = []
    for i in range(n_samples):
        sample_id = f'SAMPLE_{i+1:03d}'
        cancer_type = np.random.choice(['CRC', 'Lung', 'Breast'])
        stage = np.random.choice(['Early', 'Late'])
        sample_type = np.random.choice(['Cancer', 'Control'])
        
        # Generate features based on sample type
        if sample_type == 'Cancer':
            features = np.random.normal(0.7, 0.2, n_features)
            features = np.clip(features, 0, 1)
        else:
            features = np.random.normal(0.3, 0.2, n_features)
            features = np.clip(features, 0, 1)
        
        row = [sample_id, cancer_type, stage, sample_type] + features.tolist()
        data.append(row)
    
    # Create DataFrame
    columns = ['sample_id', 'cancer_type', 'stage', 'sample_type'] + feature_names
    df = pd.DataFrame(data, columns=columns)
    
    # Save sample data
    df.to_csv('sample_feature_matrix.tsv', sep='\t', index=False)
    print("Sample feature matrix saved as 'sample_feature_matrix.tsv'")
    
    # Load data
    pipeline.load_data('sample_feature_matrix.tsv')
    
    # Prepare features
    X, y = pipeline.prepare_features()
    
    # Train models
    pipeline.train_models()
    
    # Evaluate models
    results = pipeline.evaluate_models()
    
    # Calculate feature importance
    pipeline.calculate_feature_importance()
    
    # Save results
    pipeline.save_results()
    
    # Create summary report
    pipeline.create_summary_report()
    
    print("\nCENTAUR ML Pipeline completed successfully!")
    print("Results saved to 'ml_results/' directory")

if __name__ == "__main__":
    main()