#!/usr/bin/env python3
"""
Basic Bioinformatics Analysis Example

This script demonstrates basic bioinformatics analysis capabilities
for metagenomics and metabolomics data.

Author: Sina Berk Golech
Date: 2024
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Tuple, Optional
import logging

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Set style for plots
plt.style.use('seaborn-v0_8')
sns.set_palette("husl")


class MetagenomicsAnalyzer:
    """Class for basic metagenomics analysis."""
    
    def __init__(self, data: pd.DataFrame):
        """
        Initialize the metagenomics analyzer.
        
        Args:
            data: DataFrame with taxonomic abundance data
        """
        self.data = data
        self.sample_names = data.columns.tolist()
        self.taxa_names = data.index.tolist()
        
    def calculate_diversity(self, method: str = 'shannon') -> Dict[str, float]:
        """
        Calculate alpha diversity metrics.
        
        Args:
            method: Diversity metric ('shannon', 'simpson', 'richness')
            
        Returns:
            Dictionary with diversity values for each sample
        """
        diversity_scores = {}
        
        for sample in self.sample_names:
            abundances = self.data[sample].values
            abundances = abundances[abundances > 0]  # Remove zeros
            
            if method == 'shannon':
                # Shannon diversity index
                proportions = abundances / abundances.sum()
                diversity = -np.sum(proportions * np.log(proportions))
            elif method == 'simpson':
                # Simpson diversity index
                proportions = abundances / abundances.sum()
                diversity = 1 - np.sum(proportions ** 2)
            elif method == 'richness':
                # Species richness
                diversity = len(abundances)
            else:
                raise ValueError(f"Unknown diversity method: {method}")
                
            diversity_scores[sample] = diversity
            
        return diversity_scores
    
    def plot_abundance_heatmap(self, top_n: int = 20, figsize: Tuple[int, int] = (12, 8)):
        """
        Plot abundance heatmap for top taxa.
        
        Args:
            top_n: Number of top taxa to display
            figsize: Figure size
        """
        # Get top taxa by mean abundance
        mean_abundances = self.data.mean(axis=1)
        top_taxa = mean_abundances.nlargest(top_n).index
        
        # Create heatmap data
        heatmap_data = self.data.loc[top_taxa]
        
        # Create plot
        plt.figure(figsize=figsize)
        sns.heatmap(heatmap_data, 
                   cmap='YlOrRd', 
                   cbar_kws={'label': 'Relative Abundance (%)'},
                   xticklabels=True,
                   yticklabels=True)
        plt.title(f'Top {top_n} Taxa Abundance Heatmap')
        plt.xlabel('Samples')
        plt.ylabel('Taxa')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.show()
        
    def plot_diversity_comparison(self, methods: List[str] = None):
        """
        Plot diversity comparison across samples.
        
        Args:
            methods: List of diversity methods to compare
        """
        if methods is None:
            methods = ['shannon', 'simpson', 'richness']
            
        fig, axes = plt.subplots(1, len(methods), figsize=(15, 5))
        if len(methods) == 1:
            axes = [axes]
            
        for i, method in enumerate(methods):
            diversity_scores = self.calculate_diversity(method)
            
            samples = list(diversity_scores.keys())
            scores = list(diversity_scores.values())
            
            axes[i].bar(samples, scores, color='skyblue', alpha=0.7)
            axes[i].set_title(f'{method.title()} Diversity')
            axes[i].set_xlabel('Samples')
            axes[i].set_ylabel('Diversity Score')
            axes[i].tick_params(axis='x', rotation=45)
            
        plt.tight_layout()
        plt.show()


class MetabolomicsAnalyzer:
    """Class for basic metabolomics analysis."""
    
    def __init__(self, data: pd.DataFrame, metadata: Optional[pd.DataFrame] = None):
        """
        Initialize the metabolomics analyzer.
        
        Args:
            data: DataFrame with metabolite abundance data
            metadata: Optional metadata for samples
        """
        self.data = data
        self.metadata = metadata
        self.sample_names = data.columns.tolist()
        self.metabolite_names = data.index.tolist()
        
    def normalize_data(self, method: str = 'log2') -> pd.DataFrame:
        """
        Normalize metabolomics data.
        
        Args:
            method: Normalization method ('log2', 'zscore', 'minmax')
            
        Returns:
            Normalized data
        """
        if method == 'log2':
            # Log2 transformation with pseudocount
            normalized = np.log2(self.data + 1)
        elif method == 'zscore':
            # Z-score normalization
            normalized = (self.data - self.data.mean()) / self.data.std()
        elif method == 'minmax':
            # Min-max normalization
            normalized = (self.data - self.data.min()) / (self.data.max() - self.data.min())
        else:
            raise ValueError(f"Unknown normalization method: {method}")
            
        return normalized
    
    def find_differential_metabolites(self, 
                                    group1: List[str], 
                                    group2: List[str],
                                    threshold: float = 0.05) -> pd.DataFrame:
        """
        Find differentially abundant metabolites between two groups.
        
        Args:
            group1: Sample names for group 1
            group2: Sample names for group 2
            threshold: P-value threshold for significance
            
        Returns:
            DataFrame with differential analysis results
        """
        from scipy import stats
        
        results = []
        
        for metabolite in self.metabolite_names:
            # Get abundances for each group
            group1_data = self.data.loc[metabolite, group1].values
            group2_data = self.data.loc[metabolite, group2].values
            
            # Perform t-test
            t_stat, p_value = stats.ttest_ind(group1_data, group2_data)
            
            # Calculate fold change
            fold_change = np.mean(group2_data) / np.mean(group1_data)
            log2_fold_change = np.log2(fold_change)
            
            # Determine significance
            significant = p_value < threshold
            
            results.append({
                'metabolite': metabolite,
                'p_value': p_value,
                't_statistic': t_stat,
                'fold_change': fold_change,
                'log2_fold_change': log2_fold_change,
                'significant': significant,
                'group1_mean': np.mean(group1_data),
                'group2_mean': np.mean(group2_data)
            })
            
        return pd.DataFrame(results)
    
    def plot_volcano_plot(self, 
                         differential_results: pd.DataFrame,
                         figsize: Tuple[int, int] = (10, 8)):
        """
        Plot volcano plot for differential metabolites.
        
        Args:
            differential_results: Results from differential analysis
            figsize: Figure size
        """
        plt.figure(figsize=figsize)
        
        # Create scatter plot
        plt.scatter(differential_results['log2_fold_change'], 
                   -np.log10(differential_results['p_value']),
                   alpha=0.6, s=30)
        
        # Highlight significant metabolites
        significant = differential_results[differential_results['significant']]
        plt.scatter(significant['log2_fold_change'], 
                   -np.log10(significant['p_value']),
                   color='red', alpha=0.8, s=50, label='Significant')
        
        # Add threshold lines
        plt.axhline(-np.log10(0.05), color='gray', linestyle='--', alpha=0.5)
        plt.axvline(1, color='gray', linestyle='--', alpha=0.5)
        plt.axvline(-1, color='gray', linestyle='--', alpha=0.5)
        
        plt.xlabel('Log2 Fold Change')
        plt.ylabel('-Log10 P-value')
        plt.title('Volcano Plot: Differential Metabolites')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()


def generate_sample_data() -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Generate sample data for demonstration.
    
    Returns:
        Tuple of (metagenomics_data, metabolomics_data)
    """
    # Generate sample metagenomics data
    np.random.seed(42)
    n_samples = 10
    n_taxa = 50
    
    sample_names = [f'Sample_{i+1}' for i in range(n_samples)]
    taxa_names = [f'Taxon_{i+1}' for i in range(n_taxa)]
    
    # Generate random abundances
    metagenomics_data = pd.DataFrame(
        np.random.dirichlet(np.ones(n_taxa), size=n_samples).T * 100,
        index=taxa_names,
        columns=sample_names
    )
    
    # Generate sample metabolomics data
    n_metabolites = 100
    metabolite_names = [f'Metabolite_{i+1}' for i in range(n_metabolites)]
    
    # Generate random abundances with some structure
    metabolomics_data = pd.DataFrame(
        np.random.lognormal(mean=5, sigma=1, size=(n_metabolites, n_samples)),
        index=metabolite_names,
        columns=sample_names
    )
    
    return metagenomics_data, metabolomics_data


def main():
    """Main function to demonstrate the analysis capabilities."""
    logger.info("Starting bioinformatics analysis demonstration...")
    
    # Generate sample data
    logger.info("Generating sample data...")
    metagenomics_data, metabolomics_data = generate_sample_data()
    
    # Metagenomics analysis
    logger.info("Performing metagenomics analysis...")
    metagenomics_analyzer = MetagenomicsAnalyzer(metagenomics_data)
    
    # Calculate diversity metrics
    shannon_diversity = metagenomics_analyzer.calculate_diversity('shannon')
    logger.info(f"Shannon diversity scores: {shannon_diversity}")
    
    # Plot abundance heatmap
    metagenomics_analyzer.plot_abundance_heatmap(top_n=15)
    
    # Plot diversity comparison
    metagenomics_analyzer.plot_diversity_comparison()
    
    # Metabolomics analysis
    logger.info("Performing metabolomics analysis...")
    metabolomics_analyzer = MetabolomicsAnalyzer(metabolomics_data)
    
    # Normalize data
    normalized_data = metabolomics_analyzer.normalize_data('log2')
    logger.info("Data normalized using log2 transformation")
    
    # Create sample groups for differential analysis
    group1 = metagenomics_data.columns[:5].tolist()
    group2 = metagenomics_data.columns[5:].tolist()
    
    # Find differential metabolites
    differential_results = metabolomics_analyzer.find_differential_metabolites(
        group1, group2, threshold=0.05
    )
    
    significant_count = differential_results['significant'].sum()
    logger.info(f"Found {significant_count} significantly different metabolites")
    
    # Plot volcano plot
    metabolomics_analyzer.plot_volcano_plot(differential_results)
    
    logger.info("Analysis demonstration completed!")


if __name__ == "__main__":
    main() 