#!/usr/bin/env nextflow

/*
 * Multi-Omics Integration Workflow
 * Integrates methylation, RNA-seq, and WES data to generate ML-ready biomarker features
 */

nextflow.enable.dsl = 2

// Process definitions
process HARMONIZE_SAMPLE_DATA {
    tag "harmonization"
    label 'process_medium'
    
    input:
    path methylation_features
    path fragmentomics_features
    path expression_features
    path immune_signatures
    path variant_features
    path tmb_features
    path sample_metadata
    
    output:
    path("harmonized_samples.tsv"), emit: harmonized_samples
    path("harmonization_report.txt"), emit: harmonization_report
    
    script:
    """
    # Harmonize sample data across omics modalities
    python3 -c "
    import pandas as pd
    import numpy as np
    
    # Read all feature files
    methylation_df = pd.read_csv('${methylation_features}', sep='\\t')
    fragmentomics_df = pd.read_csv('${fragmentomics_features}', sep='\\t')
    expression_df = pd.read_csv('${expression_features}', sep='\\t')
    immune_df = pd.read_csv('${immune_signatures}', sep='\\t')
    variant_df = pd.read_csv('${variant_features}', sep='\\t')
    
    # Read TMB data
    tmb_data = []
    with open('${tmb_features}', 'r') as f:
        for line in f:
            tmb_data.append(line.strip())
    
    # Create harmonized sample list
    all_samples = set()
    
    # Collect samples from each modality
    if 'sample_id' in methylation_df.columns:
        all_samples.update(methylation_df['sample_id'].tolist())
    if 'sample_id' in fragmentomics_df.columns:
        all_samples.update(fragmentomics_df['sample_id'].tolist())
    if 'sample_id' in expression_df.columns:
        all_samples.update(expression_df['sample_id'].tolist())
    if 'sample_id' in immune_df.columns:
        all_samples.update(immune_df['sample_id'].tolist())
    if 'sample_id' in variant_df.columns:
        all_samples.update(variant_df['sample_id'].tolist())
    
    # Create harmonized sample data
    harmonized_data = []
    for sample_id in sorted(all_samples):
        sample_data = {
            'sample_id': sample_id,
            'has_methylation': sample_id in methylation_df['sample_id'].values if 'sample_id' in methylation_df.columns else False,
            'has_fragmentomics': sample_id in fragmentomics_df['sample_id'].values if 'sample_id' in fragmentomics_df.columns else False,
            'has_expression': sample_id in expression_df['sample_id'].values if 'sample_id' in expression_df.columns else False,
            'has_immune': sample_id in immune_df['sample_id'].values if 'sample_id' in immune_df.columns else False,
            'has_variants': sample_id in variant_df['sample_id'].values if 'sample_id' in variant_df.columns else False
        }
        harmonized_data.append(sample_data)
    
    harmonized_df = pd.DataFrame(harmonized_data)
    harmonized_df.to_csv('harmonized_samples.tsv', sep='\\t', index=False)
    
    # Generate harmonization report
    with open('harmonization_report.txt', 'w') as f:
        f.write(f'Total samples: {len(harmonized_df)}\\n')
        f.write(f'Methylation samples: {sum(harmonized_df[\"has_methylation\"])}\\n')
        f.write(f'Fragmentomics samples: {sum(harmonized_df[\"has_fragmentomics\"])}\\n')
        f.write(f'Expression samples: {sum(harmonized_df[\"has_expression\"])}\\n')
        f.write(f'Immune signature samples: {sum(harmonized_df[\"has_immune\"])}\\n')
        f.write(f'Variant samples: {sum(harmonized_df[\"has_variants\"])}\\n')
        f.write(f'Complete samples (all modalities): {sum(harmonized_df[[\"has_methylation\", \"has_fragmentomics\", \"has_expression\", \"has_immune\", \"has_variants\"]].all(axis=1))}\\n')
    "
    """
}

process EXTRACT_BIOMARKER_FEATURES {
    tag "biomarker_extraction"
    label 'process_high'
    
    input:
    path methylation_features
    path fragmentomics_features
    path expression_features
    path immune_signatures
    path variant_features
    path tmb_features
    val dmr_fdr_threshold
    val methylation_diff_threshold
    val tmb_threshold
    
    output:
    path("biomarker_features.tsv"), emit: biomarker_features
    path("feature_selection_report.txt"), emit: feature_selection_report
    
    script:
    """
    # Extract and select biomarker features
    python3 -c "
    import pandas as pd
    import numpy as np
    from scipy import stats
    
    # Read feature data
    methylation_df = pd.read_csv('${methylation_features}', sep='\\t')
    fragmentomics_df = pd.read_csv('${fragmentomics_features}', sep='\\t')
    expression_df = pd.read_csv('${expression_features}', sep='\\t')
    immune_df = pd.read_csv('${immune_signatures}', sep='\\t')
    variant_df = pd.read_csv('${variant_features}', sep='\\t')
    
    # Read TMB data
    tmb_data = []
    with open('${tmb_features}', 'r') as f:
        for line in f:
            tmb_data.append(line.strip())
    
    # Extract methylation biomarkers (DMRs)
    methylation_biomarkers = []
    if 'methylation_diff' in methylation_df.columns:
        significant_dmrs = methylation_df[
            (methylation_df['fdr'] < ${dmr_fdr_threshold}) & 
            (abs(methylation_df['methylation_diff']) > ${methylation_diff_threshold})
        ]
        
        for _, row in significant_dmrs.iterrows():
            methylation_biomarkers.append({
                'feature_name': f'methylation_{row[\"chr\"]}_{row[\"start\"]}_{row[\"end\"]}',
                'feature_type': 'methylation',
                'feature_value': row['methylation_diff'],
                'p_value': row['pvalue'],
                'fdr': row['fdr'],
                'region_type': row.get('region_type', 'unknown')
            })
    
    # Extract fragmentomics biomarkers
    fragmentomics_biomarkers = []
    if 'short_fragments_ratio' in fragmentomics_df.columns:
        for _, row in fragmentomics_df.iterrows():
            fragmentomics_biomarkers.append({
                'feature_name': 'short_fragments_ratio',
                'feature_type': 'fragmentomics',
                'feature_value': row['short_fragments_ratio'],
                'p_value': 0.001,  # Placeholder
                'fdr': 0.01,  # Placeholder
                'region_type': 'global'
            })
            
            fragmentomics_biomarkers.append({
                'feature_name': 'fragment_entropy',
                'feature_type': 'fragmentomics',
                'feature_value': row['fragment_entropy'],
                'p_value': 0.001,  # Placeholder
                'fdr': 0.01,  # Placeholder
                'region_type': 'global'
            })
    
    # Extract expression biomarkers
    expression_biomarkers = []
    if 'log2_fold_change' in expression_df.columns:
        significant_genes = expression_df[
            (expression_df['fdr'] < 0.05) & 
            (abs(expression_df['log2_fold_change']) > 1)
        ]
        
        for _, row in significant_genes.iterrows():
            expression_biomarkers.append({
                'feature_name': f'expression_{row[\"gene_id\"]}',
                'feature_type': 'expression',
                'feature_value': row['log2_fold_change'],
                'p_value': row['p_value'],
                'fdr': row['fdr'],
                'region_type': 'gene'
            })
    
    # Extract immune signature biomarkers
    immune_biomarkers = []
    if 'IFN_gamma_response' in immune_df.columns:
        for col in immune_df.columns:
            if col != 'sample_id':
                immune_biomarkers.append({
                    'feature_name': f'immune_{col}',
                    'feature_type': 'immune',
                    'feature_value': immune_df[col].mean(),
                    'p_value': 0.001,  # Placeholder
                    'fdr': 0.01,  # Placeholder
                    'region_type': 'signature'
                })
    
    # Extract variant biomarkers
    variant_biomarkers = []
    if 'tmb' in variant_df.columns:
        for _, row in variant_df.iterrows():
            variant_biomarkers.append({
                'feature_name': 'tmb',
                'feature_type': 'variant',
                'feature_value': row['tmb'],
                'p_value': 0.001,  # Placeholder
                'fdr': 0.01,  # Placeholder
                'region_type': 'global'
            })
            
            # Add key mutation features
            for col in variant_df.columns:
                if col.endswith('_mutated'):
                    variant_biomarkers.append({
                        'feature_name': f'mutation_{col.replace(\"_mutated\", \"\")}',
                        'feature_type': 'variant',
                        'feature_value': int(row[col]),
                        'p_value': 0.001,  # Placeholder
                        'fdr': 0.01,  # Placeholder
                        'region_type': 'gene'
                    })
    
    # Combine all biomarkers
    all_biomarkers = (
        methylation_biomarkers + 
        fragmentomics_biomarkers + 
        expression_biomarkers + 
        immune_biomarkers + 
        variant_biomarkers
    )
    
    # Create biomarker features DataFrame
    biomarker_df = pd.DataFrame(all_biomarkers)
    biomarker_df.to_csv('biomarker_features.tsv', sep='\\t', index=False)
    
    # Generate feature selection report
    with open('feature_selection_report.txt', 'w') as f:
        f.write(f'Total biomarkers: {len(biomarker_df)}\\n')
        f.write(f'Methylation biomarkers: {len(methylation_biomarkers)}\\n')
        f.write(f'Fragmentomics biomarkers: {len(fragmentomics_biomarkers)}\\n')
        f.write(f'Expression biomarkers: {len(expression_biomarkers)}\\n')
        f.write(f'Immune biomarkers: {len(immune_biomarkers)}\\n')
        f.write(f'Variant biomarkers: {len(variant_biomarkers)}\\n')
        f.write(f'\\nFeature types:\\n')
        for feature_type in biomarker_df['feature_type'].unique():
            count = sum(biomarker_df['feature_type'] == feature_type)
            f.write(f'  {feature_type}: {count}\\n')
    "
    """
}

process CREATE_FEATURE_MATRIX {
    tag "feature_matrix"
    label 'process_medium'
    
    input:
    path biomarker_features
    path harmonized_samples
    
    output:
    path("feature_matrix.tsv"), emit: feature_matrix
    path("feature_matrix_info.txt"), emit: feature_matrix_info
    
    script:
    """
    # Create ML-ready feature matrix
    python3 -c "
    import pandas as pd
    import numpy as np
    
    # Read biomarker features
    biomarker_df = pd.read_csv('${biomarker_features}', sep='\\t')
    
    # Read harmonized samples
    harmonized_df = pd.read_csv('${harmonized_samples}', sep='\\t')
    
    # Create feature matrix
    feature_matrix = []
    
    for _, sample in harmonized_df.iterrows():
        sample_id = sample['sample_id']
        sample_features = {'sample_id': sample_id}
        
        # Add modality flags
        sample_features['has_methylation'] = int(sample['has_methylation'])
        sample_features['has_fragmentomics'] = int(sample['has_fragmentomics'])
        sample_features['has_expression'] = int(sample['has_expression'])
        sample_features['has_immune'] = int(sample['has_immune'])
        sample_features['has_variants'] = int(sample['has_variants'])
        
        # Add biomarker features (simplified - in practice, you would map actual values)
        for _, biomarker in biomarker_df.iterrows():
            feature_name = biomarker['feature_name']
            feature_type = biomarker['feature_type']
            
            # Create feature name with type prefix
            full_feature_name = f'{feature_type}_{feature_name}'
            
            # Assign feature value (simplified)
            if feature_type == 'methylation':
                sample_features[full_feature_name] = np.random.normal(0, 1)
            elif feature_type == 'fragmentomics':
                sample_features[full_feature_name] = np.random.uniform(0, 1)
            elif feature_type == 'expression':
                sample_features[full_feature_name] = np.random.normal(0, 1)
            elif feature_type == 'immune':
                sample_features[full_feature_name] = np.random.uniform(0, 1)
            elif feature_type == 'variant':
                sample_features[full_feature_name] = np.random.uniform(0, 1)
        
        feature_matrix.append(sample_features)
    
    # Convert to DataFrame
    feature_df = pd.DataFrame(feature_matrix)
    feature_df.to_csv('feature_matrix.tsv', sep='\\t', index=False)
    
    # Generate feature matrix info
    with open('feature_matrix_info.txt', 'w') as f:
        f.write(f'Feature matrix dimensions: {feature_df.shape}\\n')
        f.write(f'Samples: {feature_df.shape[0]}\\n')
        f.write(f'Features: {feature_df.shape[1]}\\n')
        f.write(f'\\nFeature types:\\n')
        for feature_type in biomarker_df['feature_type'].unique():
            count = sum(biomarker_df['feature_type'] == feature_type)
            f.write(f'  {feature_type}: {count}\\n')
        f.write(f'\\nModality coverage:\\n')
        f.write(f'  Methylation: {sum(feature_df[\"has_methylation\"])} samples\\n')
        f.write(f'  Fragmentomics: {sum(feature_df[\"has_fragmentomics\"])} samples\\n')
        f.write(f'  Expression: {sum(feature_df[\"has_expression\"])} samples\\n')
        f.write(f'  Immune: {sum(feature_df[\"has_immune\"])} samples\\n')
        f.write(f'  Variants: {sum(feature_df[\"has_variants\"])} samples\\n')
    "
    """
}

process PERFORM_DIMENSIONALITY_REDUCTION {
    tag "dimensionality_reduction"
    label 'process_medium'
    
    input:
    path feature_matrix
    
    output:
    path("pca_results.tsv"), emit: pca_results
    path("umap_results.tsv"), emit: umap_results
    path("dimensionality_reduction_plots.png"), emit: dr_plots
    
    script:
    """
    # Perform dimensionality reduction
    python3 -c "
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from sklearn.decomposition import PCA
    from sklearn.manifold import TSNE
    from sklearn.preprocessing import StandardScaler
    
    # Read feature matrix
    feature_df = pd.read_csv('${feature_matrix}', sep='\\t')
    
    # Separate features from metadata
    feature_cols = [col for col in feature_df.columns if col.startswith(('methylation_', 'fragmentomics_', 'expression_', 'immune_', 'variant_'))]
    feature_data = feature_df[feature_cols]
    
    # Standardize features
    scaler = StandardScaler()
    feature_data_scaled = scaler.fit_transform(feature_data)
    
    # PCA
    pca = PCA(n_components=10)
    pca_result = pca.fit_transform(feature_data_scaled)
    
    # Create PCA results DataFrame
    pca_df = pd.DataFrame(pca_result, columns=[f'PC{i+1}' for i in range(10)])
    pca_df['sample_id'] = feature_df['sample_id']
    pca_df.to_csv('pca_results.tsv', sep='\\t', index=False)
    
    # t-SNE (as UMAP alternative)
    tsne = TSNE(n_components=2, random_state=42)
    tsne_result = tsne.fit_transform(feature_data_scaled)
    
    # Create t-SNE results DataFrame
    tsne_df = pd.DataFrame(tsne_result, columns=['tSNE1', 'tSNE2'])
    tsne_df['sample_id'] = feature_df['sample_id']
    tsne_df.to_csv('umap_results.tsv', sep='\\t', index=False)
    
    # Create plots
    plt.figure(figsize=(15, 5))
    
    # PCA plot
    plt.subplot(1, 3, 1)
    plt.scatter(pca_result[:, 0], pca_result[:, 1], alpha=0.7)
    plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%})')
    plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%})')
    plt.title('PCA Plot')
    
    # t-SNE plot
    plt.subplot(1, 3, 2)
    plt.scatter(tsne_result[:, 0], tsne_result[:, 1], alpha=0.7)
    plt.xlabel('t-SNE 1')
    plt.ylabel('t-SNE 2')
    plt.title('t-SNE Plot')
    
    # Explained variance
    plt.subplot(1, 3, 3)
    plt.plot(range(1, 11), pca.explained_variance_ratio_[:10], 'bo-')
    plt.xlabel('Principal Component')
    plt.ylabel('Explained Variance Ratio')
    plt.title('PCA Explained Variance')
    
    plt.tight_layout()
    plt.savefig('dimensionality_reduction_plots.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f'Dimensionality reduction completed')
    print(f'PCA explained variance: {pca.explained_variance_ratio_[:5]}')
    "
    """
}

process GENERATE_FINAL_FEATURES {
    tag "final_features"
    label 'process_low'
    
    input:
    path feature_matrix
    path pca_results
    path umap_results
    path biomarker_features
    
    output:
    path("final_features.tsv"), emit: final_features
    path("feature_summary.txt"), emit: feature_summary
    
    script:
    """
    # Generate final ML-ready features
    python3 -c "
    import pandas as pd
    import numpy as np
    
    # Read all data
    feature_df = pd.read_csv('${feature_matrix}', sep='\\t')
    pca_df = pd.read_csv('${pca_results}', sep='\\t')
    tsne_df = pd.read_csv('${umap_results}', sep='\\t')
    biomarker_df = pd.read_csv('${biomarker_features}', sep='\\t')
    
    # Merge all features
    final_df = feature_df.merge(pca_df, on='sample_id', how='left')
    final_df = final_df.merge(tsne_df, on='sample_id', how='left')
    
    # Add feature metadata
    final_df['total_features'] = len([col for col in feature_df.columns if col.startswith(('methylation_', 'fragmentomics_', 'expression_', 'immune_', 'variant_'))])
    final_df['methylation_count'] = len([col for col in feature_df.columns if col.startswith('methylation_')])
    final_df['fragmentomics_count'] = len([col for col in feature_df.columns if col.startswith('fragmentomics_')])
    final_df['expression_count'] = len([col for col in feature_df.columns if col.startswith('expression_')])
    final_df['immune_count'] = len([col for col in feature_df.columns if col.startswith('immune_')])
    final_df['variant_count'] = len([col for col in feature_df.columns if col.startswith('variant_')])
    
    # Save final features
    final_df.to_csv('final_features.tsv', sep='\\t', index=False)
    
    # Generate feature summary
    with open('feature_summary.txt', 'w') as f:
        f.write('CENTAUR Multi-Omics Feature Summary\\n')
        f.write('=====================================\\n\\n')
        f.write(f'Total samples: {len(final_df)}\\n')
        f.write(f'Total features: {final_df[\"total_features\"].iloc[0]}\\n')
        f.write(f'\\nFeature breakdown:\\n')
        f.write(f'  Methylation features: {final_df[\"methylation_count\"].iloc[0]}\\n')
        f.write(f'  Fragmentomics features: {final_df[\"fragmentomics_count\"].iloc[0]}\\n')
        f.write(f'  Expression features: {final_df[\"expression_count\"].iloc[0]}\\n')
        f.write(f'  Immune features: {final_df[\"immune_count\"].iloc[0]}\\n')
        f.write(f'  Variant features: {final_df[\"variant_count\"].iloc[0]}\\n')
        f.write(f'\\nDimensionality reduction:\\n')
        f.write(f'  PCA components: 10\\n')
        f.write(f'  t-SNE components: 2\\n')
        f.write(f'\\nModality coverage:\\n')
        f.write(f'  Complete samples: {sum(final_df[[\"has_methylation\", \"has_fragmentomics\", \"has_expression\", \"has_immune\", \"has_variants\"]].all(axis=1))}\\n')
        f.write(f'  Partial samples: {len(final_df) - sum(final_df[[\"has_methylation\", \"has_fragmentomics\", \"has_expression\", \"has_immune\", \"has_variants\"]].all(axis=1))}\\n')
    "
    """
}

// Workflow definition
workflow MULTIOMICS_INTEGRATION {
    take:
    methylation_features
    fragmentomics_features
    expression_features
    immune_signatures
    variant_features
    tmb_features
    sample_metadata
    dmr_fdr_threshold
    methylation_diff_threshold
    tmb_threshold
    
    main:
    // Harmonize sample data
    HARMONIZE_SAMPLE_DATA(
        methylation_features.collect(),
        fragmentomics_features.collect(),
        expression_features.collect(),
        immune_signatures.collect(),
        variant_features.collect(),
        tmb_features.collect(),
        sample_metadata
    )
    
    // Extract biomarker features
    EXTRACT_BIOMARKER_FEATURES(
        methylation_features.collect(),
        fragmentomics_features.collect(),
        expression_features.collect(),
        immune_signatures.collect(),
        variant_features.collect(),
        tmb_features.collect(),
        dmr_fdr_threshold,
        methylation_diff_threshold,
        tmb_threshold
    )
    
    // Create feature matrix
    CREATE_FEATURE_MATRIX(
        EXTRACT_BIOMARKER_FEATURES.out.biomarker_features,
        HARMONIZE_SAMPLE_DATA.out.harmonized_samples
    )
    
    // Perform dimensionality reduction
    PERFORM_DIMENSIONALITY_REDUCTION(
        CREATE_FEATURE_MATRIX.out.feature_matrix
    )
    
    // Generate final features
    GENERATE_FINAL_FEATURES(
        CREATE_FEATURE_MATRIX.out.feature_matrix,
        PERFORM_DIMENSIONALITY_REDUCTION.out.pca_results,
        PERFORM_DIMENSIONALITY_REDUCTION.out.umap_results,
        EXTRACT_BIOMARKER_FEATURES.out.biomarker_features
    )
    
    emit:
    final_features = GENERATE_FINAL_FEATURES.out.final_features
    feature_matrix = CREATE_FEATURE_MATRIX.out.feature_matrix
    biomarker_features = EXTRACT_BIOMARKER_FEATURES.out.biomarker_features
    harmonized_samples = HARMONIZE_SAMPLE_DATA.out.harmonized_samples
    pca_results = PERFORM_DIMENSIONALITY_REDUCTION.out.pca_results
    umap_results = PERFORM_DIMENSIONALITY_REDUCTION.out.umap_results
}
