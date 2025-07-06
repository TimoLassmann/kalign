# Kalign Ecosystem Integration Guide

This guide demonstrates how to seamlessly integrate Kalign with the broader bioinformatics ecosystem, including Biopython, scikit-bio, pandas, matplotlib, and other popular tools.

## Table of Contents

- [Biopython Integration](#biopython-integration)
- [scikit-bio Integration](#scikit-bio-integration)
- [Pandas Integration](#pandas-integration)
- [Visualization with Matplotlib](#visualization-with-matplotlib)
- [Phylogenetic Analysis](#phylogenetic-analysis)
- [Machine Learning Workflows](#machine-learning-workflows)
- [Jupyter Notebook Integration](#jupyter-notebook-integration)
- [Web Application Integration](#web-application-integration)
- [High-Throughput Pipelines](#high-throughput-pipelines)

## Biopython Integration

### Basic Biopython Workflow

```python
import kalign
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Read sequences with Biopython
records = list(SeqIO.parse("sequences.fasta", "fasta"))
sequences = [str(record.seq) for record in records]
ids = [record.id for record in records]

# Align with Kalign and get Biopython object
alignment = kalign.align(sequences, fmt="biopython", ids=ids)

# Now use Biopython's rich functionality
print(f"Alignment length: {alignment.get_alignment_length()}")
print(f"Number of sequences: {len(alignment)}")

# Export to various formats
AlignIO.write(alignment, "output.clustal", "clustal")
AlignIO.write(alignment, "output.phylip", "phylip")
AlignIO.write(alignment, "output.nexus", "nexus")
```

### Advanced Biopython Analysis

```python
import kalign
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Align import AlignInfo
import numpy as np

def biopython_analysis_pipeline(fasta_file):
    """Complete analysis pipeline using Biopython."""
    
    # Step 1: Read and align sequences
    print("🧬 Reading sequences...")
    records = list(SeqIO.parse(fasta_file, "fasta"))
    sequences = [str(record.seq) for record in records]
    ids = [record.id for record in records]
    
    print(f"Found {len(sequences)} sequences")
    
    # Step 2: Align with Kalign
    print("🔗 Aligning sequences...")
    alignment = kalign.align(sequences, fmt="biopython", ids=ids)
    
    # Step 3: Create alignment summary
    print("📊 Analyzing alignment...")
    summary_align = AlignInfo.SummaryInfo(alignment)
    
    # Calculate consensus
    consensus = summary_align.gap_consensus(threshold=0.7, ambiguous='N')
    print(f"Consensus (70%): {consensus}")
    
    # Position specific scoring matrix
    pssm = summary_align.pos_specific_score_matrix(consensus)
    
    # Step 4: Calculate distance matrix
    print("📏 Calculating distances...")
    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(alignment)
    
    # Step 5: Build phylogenetic tree
    print("🌳 Building tree...")
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(distance_matrix)
    
    # Step 6: Save results
    AlignIO.write(alignment, "analysis_alignment.fasta", "fasta")
    
    print("✅ Analysis complete!")
    return alignment, consensus, distance_matrix, tree

# Usage
# alignment, consensus, distances, tree = biopython_analysis_pipeline("sequences.fasta")
```

### Sequence Feature Analysis

```python
import kalign
from Bio import SeqIO
from Bio.SeqUtils import GC, molecular_weight
from Bio.SeqUtils.ProtParam import ProteinAnalysis

def analyze_protein_features(sequences, ids=None):
    """Analyze protein features before and after alignment."""
    
    if ids is None:
        ids = [f"seq{i}" for i in range(len(sequences))]
    
    print("🔬 Protein Feature Analysis")
    
    # Analyze original sequences
    print("\n📋 Original Sequence Features:")
    original_features = []
    
    for i, (seq, seq_id) in enumerate(zip(sequences, ids)):
        try:
            analysis = ProteinAnalysis(seq)
            features = {
                'id': seq_id,
                'length': len(seq),
                'molecular_weight': analysis.molecular_weight(),
                'aromaticity': analysis.aromaticity(),
                'instability_index': analysis.instability_index(),
                'isoelectric_point': analysis.isoelectric_point(),
                'secondary_structure_fraction': analysis.secondary_structure_fraction()
            }
            original_features.append(features)
            
            print(f"  {seq_id}: {len(seq)} aa, MW={features['molecular_weight']:.1f}, "
                  f"pI={features['isoelectric_point']:.2f}")
            
        except Exception as e:
            print(f"  {seq_id}: Error analyzing - {e}")
    
    # Align sequences
    print("\n🔗 Aligning sequences...")
    alignment = kalign.align(sequences, seq_type="protein", fmt="biopython", ids=ids)
    
    # Analyze conserved regions
    print("\n🎯 Conservation Analysis:")
    conservation_scores = []
    
    for i in range(alignment.get_alignment_length()):
        column = alignment[:, i]
        unique_chars = set(str(column).replace('-', ''))
        if len(unique_chars) == 1 and '-' not in unique_chars:
            conservation_scores.append(1.0)  # Fully conserved
        elif len(unique_chars) == 0:
            conservation_scores.append(0.0)  # All gaps
        else:
            # Calculate conservation score
            non_gap_chars = [c for c in str(column) if c != '-']
            if non_gap_chars:
                most_common = max(set(non_gap_chars), key=non_gap_chars.count)
                score = non_gap_chars.count(most_common) / len(non_gap_chars)
                conservation_scores.append(score)
            else:
                conservation_scores.append(0.0)
    
    # Find highly conserved regions
    conserved_regions = []
    threshold = 0.8
    in_region = False
    region_start = 0
    
    for i, score in enumerate(conservation_scores):
        if score >= threshold and not in_region:
            region_start = i
            in_region = True
        elif score < threshold and in_region:
            conserved_regions.append((region_start, i-1))
            in_region = False
    
    if in_region:
        conserved_regions.append((region_start, len(conservation_scores)-1))
    
    print(f"Found {len(conserved_regions)} highly conserved regions (≥{threshold:.0%}):")
    for start, end in conserved_regions:
        print(f"  Position {start+1}-{end+1} (length: {end-start+1})")
    
    return alignment, original_features, conservation_scores, conserved_regions

# Usage
protein_sequences = [
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWUAIFRRVVSAEFQRQPVHQSYLNTVLGSQGKL",
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWUAIFRRVVSAEFQRQPVHQSYLNTVLGSQGKL"
]

alignment, features, conservation, regions = analyze_protein_features(protein_sequences)
```

## scikit-bio Integration

### Basic scikit-bio Workflow

```python
import kalign
import skbio
from skbio import DNA, Protein, TabularMSA
import numpy as np

# Align sequences and get scikit-bio object
sequences = ["ATCGATCGATCG", "ATCGTCGATCG", "ATCGATCATCG"]
alignment = kalign.align(sequences, fmt="skbio")

# Use scikit-bio's analysis tools
print(f"Consensus: {alignment.consensus()}")

# Calculate conservation
conservation = alignment.conservation()
print(f"Conservation scores: {conservation}")

# Calculate Shannon entropy
entropy = alignment.conservation(metric='shannon_uncertainty')
print(f"Shannon entropy: {entropy}")

# Export in different formats
alignment.write("output.fasta", format="fasta")
```

### Advanced scikit-bio Analysis

```python
import kalign
import skbio
import numpy as np
import pandas as pd
from scipy import stats
from skbio.stats.composition import ancom
from skbio.diversity import alpha_diversity, beta_diversity

def skbio_diversity_analysis(sequences, ids=None, sequence_type="dna"):
    """Comprehensive diversity analysis using scikit-bio."""
    
    if ids is None:
        ids = [f"seq{i}" for i in range(len(sequences))]
    
    print("🌈 Diversity Analysis with scikit-bio")
    
    # Step 1: Align sequences
    print("🔗 Aligning sequences...")
    alignment = kalign.align(sequences, fmt="skbio")
    
    # Step 2: Calculate conservation metrics
    print("📊 Calculating conservation metrics...")
    conservation = alignment.conservation()
    shannon_entropy = alignment.conservation(metric='shannon_uncertainty')
    
    print(f"Mean conservation: {np.mean(conservation):.3f}")
    print(f"Mean Shannon entropy: {np.mean(shannon_entropy):.3f}")
    
    # Step 3: Position-wise analysis
    print("📍 Position-wise analysis...")
    position_stats = []
    
    for i in range(len(alignment[0])):
        column = [str(seq)[i] for seq in alignment]
        
        # Count character frequencies
        char_counts = {}
        for char in column:
            if char != '-':  # Exclude gaps
                char_counts[char] = char_counts.get(char, 0) + 1
        
        if char_counts:
            # Calculate diversity metrics
            total_chars = sum(char_counts.values())
            frequencies = np.array(list(char_counts.values())) / total_chars
            
            # Shannon diversity
            shannon_div = -np.sum(frequencies * np.log(frequencies))
            
            # Simpson diversity
            simpson_div = 1 - np.sum(frequencies ** 2)
            
            position_stats.append({
                'position': i + 1,
                'gap_fraction': column.count('-') / len(column),
                'n_unique_chars': len(char_counts),
                'shannon_diversity': shannon_div,
                'simpson_diversity': simpson_div,
                'conservation': conservation[i],
                'entropy': shannon_entropy[i]
            })
    
    # Step 4: Create summary DataFrame
    df_positions = pd.DataFrame(position_stats)
    
    # Step 5: Identify variable regions
    variable_threshold = df_positions['shannon_diversity'].quantile(0.75)
    variable_positions = df_positions[
        df_positions['shannon_diversity'] > variable_threshold
    ]
    
    print(f"\n🎯 Variable regions (top 25% Shannon diversity):")
    print(f"Found {len(variable_positions)} variable positions")
    
    # Group consecutive variable positions into regions
    variable_regions = []
    if len(variable_positions) > 0:
        positions = variable_positions['position'].values
        
        # Find consecutive runs
        runs = []
        current_run = [positions[0]]
        
        for i in range(1, len(positions)):
            if positions[i] == positions[i-1] + 1:
                current_run.append(positions[i])
            else:
                runs.append(current_run)
                current_run = [positions[i]]
        runs.append(current_run)
        
        # Convert to regions
        for run in runs:
            if len(run) >= 3:  # Minimum region size
                variable_regions.append((min(run), max(run)))
        
        print(f"Variable regions (≥3 positions):")
        for start, end in variable_regions:
            print(f"  Positions {start}-{end} (length: {end-start+1})")
    
    # Step 6: Calculate pairwise distances
    print("\n📏 Calculating pairwise distances...")
    
    # Convert alignment to distance matrix
    distance_matrix = []
    for i in range(len(alignment)):
        row = []
        for j in range(len(alignment)):
            if i == j:
                row.append(0.0)
            else:
                seq1 = str(alignment[i]).replace('-', '')
                seq2 = str(alignment[j]).replace('-', '')
                
                # Calculate Hamming distance
                min_len = min(len(seq1), len(seq2))
                if min_len > 0:
                    differences = sum(c1 != c2 for c1, c2 in zip(seq1[:min_len], seq2[:min_len]))
                    distance = differences / min_len
                else:
                    distance = 1.0
                
                row.append(distance)
        distance_matrix.append(row)
    
    distance_matrix = np.array(distance_matrix)
    
    # Summary statistics
    print(f"\n📈 Summary Statistics:")
    print(f"  Alignment length: {len(alignment[0])}")
    print(f"  Number of sequences: {len(alignment)}")
    print(f"  Mean pairwise distance: {np.mean(distance_matrix[np.triu_indices_from(distance_matrix, k=1)]):.3f}")
    print(f"  Mean gap fraction: {df_positions['gap_fraction'].mean():.3f}")
    print(f"  Most variable position: {df_positions.loc[df_positions['shannon_diversity'].idxmax(), 'position']}")
    print(f"  Most conserved position: {df_positions.loc[df_positions['shannon_diversity'].idxmin(), 'position']}")
    
    return alignment, df_positions, variable_regions, distance_matrix

# Usage
dna_sequences = [
    "ATCGATCGATCGATCGATCG",
    "ATCGATCGTCGATCGATCG",
    "ATCGATCGATCATCGATCG",
    "ATCGATCGAGCTCGATCG",
    "ATCGATCGATCGATCAACG"
]

alignment, position_stats, variable_regions, distances = skbio_diversity_analysis(dna_sequences)
```

## Pandas Integration

### Alignment Data Analysis with Pandas

```python
import kalign
import pandas as pd
import numpy as np

def alignment_to_dataframe(sequences, ids=None):
    """Convert alignment to pandas DataFrame for analysis."""
    
    if ids is None:
        ids = [f"seq{i}" for i in range(len(sequences))]
    
    # Align sequences
    aligned = kalign.align(sequences)
    
    # Convert to DataFrame
    alignment_data = []
    for i, (seq_id, seq) in enumerate(zip(ids, aligned)):
        for j, char in enumerate(seq):
            alignment_data.append({
                'sequence_id': seq_id,
                'sequence_index': i,
                'position': j + 1,
                'character': char,
                'is_gap': char == '-'
            })
    
    df = pd.DataFrame(alignment_data)
    
    # Add summary statistics
    stats = kalign.utils.alignment_stats(aligned)
    
    print(f"📊 Alignment DataFrame created:")
    print(f"  Shape: {df.shape}")
    print(f"  Alignment length: {stats['length']}")
    print(f"  Number of sequences: {stats['n_sequences']}")
    print(f"  Total positions: {len(df)}")
    
    return df, aligned, stats

def analyze_alignment_dataframe(df):
    """Perform comprehensive analysis on alignment DataFrame."""
    
    print("🔬 Analyzing alignment with pandas...")
    
    # Position-wise analysis
    position_summary = df.groupby('position').agg({
        'is_gap': ['sum', 'mean'],
        'character': ['nunique', lambda x: x.value_counts().iloc[0] if len(x) > 0 else 0]
    }).round(3)
    
    position_summary.columns = ['gap_count', 'gap_fraction', 'unique_chars', 'most_common_count']
    
    # Sequence-wise analysis
    sequence_summary = df.groupby('sequence_id').agg({
        'is_gap': ['sum', 'mean'],
        'character': 'count'
    }).round(3)
    
    sequence_summary.columns = ['total_gaps', 'gap_fraction', 'total_positions']
    
    # Find conserved positions (no gaps, single character)
    conserved_positions = position_summary[
        (position_summary['gap_count'] == 0) & 
        (position_summary['unique_chars'] == 1)
    ]
    
    # Find variable positions (multiple characters, low gaps)
    variable_positions = position_summary[
        (position_summary['unique_chars'] > 2) & 
        (position_summary['gap_fraction'] < 0.5)
    ]
    
    print(f"\n📈 Analysis Results:")
    print(f"  Conserved positions: {len(conserved_positions)}")
    print(f"  Variable positions: {len(variable_positions)}")
    print(f"  Mean gap fraction per position: {position_summary['gap_fraction'].mean():.3f}")
    print(f"  Positions with gaps: {(position_summary['gap_count'] > 0).sum()}")
    
    return position_summary, sequence_summary, conserved_positions, variable_positions

# Example usage
sequences = [
    "ATCGATCGATCGATCGATCG",
    "ATCGATCGTCGATCGATCG", 
    "ATCGATCGATCATCGATCG",
    "ATCGATCGAGCTCGATCG"
]
ids = ["Human", "Mouse", "Rat", "Dog"]

df, aligned, stats = alignment_to_dataframe(sequences, ids)
pos_summary, seq_summary, conserved, variable = analyze_alignment_dataframe(df)

# Display results
print("\n🎯 Top 5 most variable positions:")
print(variable.head())

print("\n🔒 Sample conserved positions:")
print(conserved.head())
```

### Comparative Analysis with Pandas

```python
import kalign
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def compare_multiple_alignments(sequence_groups, group_names):
    """Compare multiple alignment groups using pandas."""
    
    results = []
    
    for group_name, sequences in zip(group_names, sequence_groups):
        print(f"Analyzing {group_name}...")
        
        # Align sequences
        aligned = kalign.align(sequences)
        stats = kalign.utils.alignment_stats(aligned)
        
        # Calculate additional metrics
        pairwise_identity = kalign.utils.pairwise_identity_matrix(aligned)
        mean_identity = np.mean(pairwise_identity[np.triu_indices_from(pairwise_identity, k=1)])
        
        # Sequence length statistics
        original_lengths = [len(seq.replace('-', '')) for seq in sequences]
        
        results.append({
            'group': group_name,
            'n_sequences': len(sequences),
            'alignment_length': stats['length'],
            'gap_fraction': stats['gap_fraction'],
            'conservation': stats['conservation'],
            'mean_pairwise_identity': mean_identity,
            'mean_original_length': np.mean(original_lengths),
            'std_original_length': np.std(original_lengths),
            'min_original_length': min(original_lengths),
            'max_original_length': max(original_lengths)
        })
    
    # Create comparison DataFrame
    df_comparison = pd.DataFrame(results)
    
    print("\n📊 Comparison Results:")
    print(df_comparison.round(3))
    
    # Statistical analysis
    print(f"\n📈 Statistical Summary:")
    print(f"Most conserved group: {df_comparison.loc[df_comparison['conservation'].idxmax(), 'group']}")
    print(f"Highest identity group: {df_comparison.loc[df_comparison['mean_pairwise_identity'].idxmax(), 'group']}")
    print(f"Most variable lengths: {df_comparison.loc[df_comparison['std_original_length'].idxmax(), 'group']}")
    
    return df_comparison

# Example: Compare different sequence families
protein_family_1 = [
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPIL",
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPIS",
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPIM"
]

protein_family_2 = [
    "ACDEFGHIKLMNPQRSTVWY",
    "ACDEFGHIKLMNPQRSTVWF",
    "ACDEFGHIKLMNPQRSTVWH"
]

dna_family = [
    "ATCGATCGATCGATCGATCG",
    "ATCGATCGTCGATCGATCG",
    "ATCGATCGATCATCGATCG"
]

sequence_groups = [protein_family_1, protein_family_2, dna_family]
group_names = ["Protein_Family_A", "Protein_Family_B", "DNA_Sequences"]

comparison_df = compare_multiple_alignments(sequence_groups, group_names)
```

## Visualization with Matplotlib

### Alignment Visualization

```python
import kalign
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.colors import ListedColormap

def visualize_alignment(sequences, ids=None, figsize=(12, 8)):
    """Create comprehensive alignment visualization."""
    
    if ids is None:
        ids = [f"Seq{i+1}" for i in range(len(sequences))]
    
    # Align sequences
    aligned = kalign.align(sequences)
    
    # Convert to numeric matrix for visualization
    char_to_num = {'-': 0}  # Gap
    all_chars = set(''.join(aligned).replace('-', ''))
    for i, char in enumerate(sorted(all_chars), 1):
        char_to_num[char] = i
    
    # Create numeric matrix
    matrix = np.array([[char_to_num[char] for char in seq] for seq in aligned])
    
    # Create figure with subplots
    fig, axes = plt.subplots(3, 2, figsize=figsize)
    fig.suptitle('Kalign Alignment Analysis', fontsize=16, fontweight='bold')
    
    # 1. Alignment heatmap
    ax1 = axes[0, 0]
    im1 = ax1.imshow(matrix, aspect='auto', cmap='tab20')
    ax1.set_title('Alignment Heatmap')
    ax1.set_xlabel('Position')
    ax1.set_ylabel('Sequence')
    ax1.set_yticks(range(len(ids)))
    ax1.set_yticklabels(ids)
    
    # 2. Conservation plot
    ax2 = axes[0, 1]
    conservation = []
    for j in range(matrix.shape[1]):
        column = matrix[:, j]
        non_gap = column[column != 0]
        if len(non_gap) > 0:
            unique_vals = len(set(non_gap))
            conservation.append(1.0 / unique_vals if unique_vals > 0 else 0)
        else:
            conservation.append(0)
    
    ax2.plot(range(1, len(conservation) + 1), conservation, 'b-', linewidth=2)
    ax2.set_title('Position-wise Conservation')
    ax2.set_xlabel('Position')
    ax2.set_ylabel('Conservation Score')
    ax2.grid(True, alpha=0.3)
    
    # 3. Gap analysis
    ax3 = axes[1, 0]
    gap_fraction = [(row == 0).sum() / len(row) for row in matrix]
    bars = ax3.bar(range(len(ids)), gap_fraction, color='orange', alpha=0.7)
    ax3.set_title('Gap Fraction per Sequence')
    ax3.set_xlabel('Sequence')
    ax3.set_ylabel('Gap Fraction')
    ax3.set_xticks(range(len(ids)))
    ax3.set_xticklabels(ids, rotation=45)
    
    # Add values on bars
    for i, (bar, value) in enumerate(zip(bars, gap_fraction)):
        ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                f'{value:.2f}', ha='center', va='bottom', fontsize=8)
    
    # 4. Position-wise gap count
    ax4 = axes[1, 1]
    gap_counts = [(matrix[:, j] == 0).sum() for j in range(matrix.shape[1])]
    ax4.plot(range(1, len(gap_counts) + 1), gap_counts, 'r-', linewidth=2)
    ax4.fill_between(range(1, len(gap_counts) + 1), gap_counts, alpha=0.3, color='red')
    ax4.set_title('Gaps per Position')
    ax4.set_xlabel('Position')
    ax4.set_ylabel('Number of Gaps')
    ax4.grid(True, alpha=0.3)
    
    # 5. Pairwise identity heatmap
    ax5 = axes[2, 0]
    identity_matrix = kalign.utils.pairwise_identity_matrix(aligned)
    im5 = ax5.imshow(identity_matrix, cmap='RdYlBu_r', vmin=0, vmax=1)
    ax5.set_title('Pairwise Identity Matrix')
    ax5.set_xticks(range(len(ids)))
    ax5.set_yticks(range(len(ids)))
    ax5.set_xticklabels(ids, rotation=45)
    ax5.set_yticklabels(ids)
    
    # Add colorbar
    cbar5 = plt.colorbar(im5, ax=ax5, shrink=0.8)
    cbar5.set_label('Identity')
    
    # Add identity values to cells
    for i in range(len(ids)):
        for j in range(len(ids)):
            ax5.text(j, i, f'{identity_matrix[i, j]:.2f}',
                    ha='center', va='center',
                    color='white' if identity_matrix[i, j] < 0.5 else 'black',
                    fontsize=8)
    
    # 6. Summary statistics
    ax6 = axes[2, 1]
    ax6.axis('off')
    
    stats = kalign.utils.alignment_stats(aligned)
    summary_text = f"""
    Alignment Statistics
    
    Sequences: {stats['n_sequences']}
    Length: {stats['length']}
    Gap Fraction: {stats['gap_fraction']:.2%}
    Conservation: {stats['conservation']:.2%}
    Avg Identity: {stats['identity']:.2%}
    
    Most Similar Pair:
    {np.unravel_index(np.argmax(identity_matrix + np.eye(len(ids)) * -1), identity_matrix.shape)}
    Identity: {np.max(identity_matrix - np.eye(len(ids))):.2%}
    
    Least Similar Pair:
    {np.unravel_index(np.argmin(identity_matrix + np.eye(len(ids)) * 2), identity_matrix.shape)}
    Identity: {np.min(identity_matrix):.2%}
    """
    
    ax6.text(0.1, 0.9, summary_text, transform=ax6.transAxes,
             fontsize=10, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.5))
    
    plt.tight_layout()
    return fig, axes

# Example usage
sequences = [
    "ATCGATCGATCGATCGATCG",
    "ATCGATCGTCGATCGATCG",
    "ATCGATCGATCATCGATCG",
    "ATCGATCGAGCTCGATCG",
    "ATCGATCGATCGATCAACG"
]
ids = ["Human", "Mouse", "Rat", "Dog", "Cat"]

fig, axes = visualize_alignment(sequences, ids)
plt.show()
```

### Interactive Alignment Viewer

```python
import kalign
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.widgets import Slider
import numpy as np

def interactive_alignment_viewer(sequences, ids=None):
    """Create interactive alignment viewer with zoom and scroll."""
    
    if ids is None:
        ids = [f"Seq{i+1}" for i in range(len(sequences))]
    
    # Align sequences
    aligned = kalign.align(sequences)
    
    # Color mapping for different characters
    unique_chars = sorted(set(''.join(aligned)))
    colors = plt.cm.tab20(np.linspace(0, 1, len(unique_chars)))
    char_colors = dict(zip(unique_chars, colors))
    char_colors['-'] = (0.9, 0.9, 0.9, 1.0)  # Light gray for gaps
    
    # Initial parameters
    window_size = min(50, len(aligned[0]))
    start_pos = 0
    
    # Create figure
    fig, ax = plt.subplots(figsize=(15, 8))
    plt.subplots_adjust(bottom=0.2)
    
    def update_display(start=0, window=50):
        """Update the alignment display."""
        ax.clear()
        
        end_pos = min(start + window, len(aligned[0]))
        
        # Display alignment segment
        for i, (seq_id, seq) in enumerate(zip(ids, aligned)):
            y_pos = len(ids) - i - 1
            
            for j, char in enumerate(seq[start:end_pos]):
                rect = patches.Rectangle(
                    (j, y_pos), 1, 1,
                    facecolor=char_colors[char],
                    edgecolor='black',
                    linewidth=0.1
                )
                ax.add_patch(rect)
                
                # Add character text
                ax.text(j + 0.5, y_pos + 0.5, char,
                       ha='center', va='center',
                       fontsize=8, fontweight='bold')
        
        # Formatting
        ax.set_xlim(0, end_pos - start)
        ax.set_ylim(0, len(ids))
        ax.set_xlabel(f'Position ({start + 1} - {end_pos})')
        ax.set_ylabel('Sequences')
        ax.set_title(f'Kalign Alignment Viewer (Window: {window} positions)')
        
        # Set y-tick labels
        ax.set_yticks([i + 0.5 for i in range(len(ids))])
        ax.set_yticklabels(reversed(ids))
        
        # Set x-tick labels
        x_ticks = range(0, end_pos - start, max(1, (end_pos - start) // 10))
        ax.set_xticks([x + 0.5 for x in x_ticks])
        ax.set_xticklabels([start + x + 1 for x in x_ticks])
        
        plt.draw()
    
    # Create sliders
    ax_pos = plt.axes([0.1, 0.1, 0.5, 0.03])
    ax_window = plt.axes([0.1, 0.05, 0.5, 0.03])
    
    slider_pos = Slider(ax_pos, 'Start Position', 0, 
                       max(0, len(aligned[0]) - window_size), 
                       valinit=start_pos, valfmt='%d')
    
    slider_window = Slider(ax_window, 'Window Size', 10, 
                          min(100, len(aligned[0])), 
                          valinit=window_size, valfmt='%d')
    
    def update_sliders(val):
        start = int(slider_pos.val)
        window = int(slider_window.val)
        
        # Update slider limits
        slider_pos.valmax = max(0, len(aligned[0]) - window)
        if start > slider_pos.valmax:
            start = slider_pos.valmax
            slider_pos.set_val(start)
        
        update_display(start, window)
    
    slider_pos.on_changed(update_sliders)
    slider_window.on_changed(update_sliders)
    
    # Initial display
    update_display(start_pos, window_size)
    
    # Add legend
    legend_elements = [patches.Patch(facecolor=color, label=char) 
                      for char, color in char_colors.items() if char in ''.join(aligned[:5])]
    ax.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.15, 1))
    
    return fig, ax, slider_pos, slider_window

# Example usage
sequences = [
    "ATCGATCGATCGATCGATCGATCGATCGATCGATCG",
    "ATCGATCGTCGATCGATCGATCGATCGTCGATCG",
    "ATCGATCGATCATCGATCGATCGATCGATCATCG",
    "ATCGATCGAGCTCGATCGATCGATCGAGCTCG"
]
ids = ["Sequence_A", "Sequence_B", "Sequence_C", "Sequence_D"]

fig, ax, pos_slider, window_slider = interactive_alignment_viewer(sequences, ids)
plt.show()
```

## Phylogenetic Analysis

### Tree Construction from Alignments

```python
import kalign
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Phylo import draw
import matplotlib.pyplot as plt

def phylogenetic_analysis(sequences, ids=None, method='upgma'):
    """Construct phylogenetic tree from aligned sequences."""
    
    if ids is None:
        ids = [f"seq{i}" for i in range(len(sequences))]
    
    print(f"🌳 Phylogenetic Analysis ({method.upper()})")
    
    # Step 1: Align sequences
    print("🔗 Aligning sequences...")
    alignment = kalign.align(sequences, fmt="biopython", ids=ids)
    
    # Step 2: Calculate distance matrix
    print("📏 Calculating distances...")
    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(alignment)
    
    print("Distance Matrix:")
    print(distance_matrix)
    
    # Step 3: Construct tree
    print(f"🌲 Building {method.upper()} tree...")
    constructor = DistanceTreeConstructor()
    
    if method.lower() == 'upgma':
        tree = constructor.upgma(distance_matrix)
    elif method.lower() == 'nj':
        tree = constructor.nj(distance_matrix)
    else:
        raise ValueError("Method must be 'upgma' or 'nj'")
    
    # Step 4: Visualize tree
    print("🎨 Creating visualization...")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Tree visualization
    Phylo.draw(tree, axes=ax1, do_show=False)
    ax1.set_title(f'{method.upper()} Phylogenetic Tree')
    
    # Distance matrix heatmap
    import numpy as np
    n = len(ids)
    matrix_array = np.zeros((n, n))
    
    for i in range(n):
        for j in range(n):
            if i == j:
                matrix_array[i, j] = 0
            else:
                # Get distance from matrix
                matrix_array[i, j] = distance_matrix[ids[i], ids[j]]
    
    im = ax2.imshow(matrix_array, cmap='RdYlBu_r')
    ax2.set_title('Distance Matrix Heatmap')
    ax2.set_xticks(range(n))
    ax2.set_yticks(range(n))
    ax2.set_xticklabels(ids, rotation=45)
    ax2.set_yticklabels(ids)
    
    # Add distance values
    for i in range(n):
        for j in range(n):
            ax2.text(j, i, f'{matrix_array[i, j]:.3f}',
                    ha='center', va='center',
                    color='white' if matrix_array[i, j] > 0.5 else 'black')
    
    plt.colorbar(im, ax=ax2, label='Distance')
    plt.tight_layout()
    
    # Step 5: Tree statistics
    print(f"\n📊 Tree Statistics:")
    print(f"  Total branch length: {tree.total_branch_length():.4f}")
    print(f"  Number of terminals: {len(tree.get_terminals())}")
    print(f"  Tree depth: {tree.depths().get(tree.root, 0):.4f}")
    
    return tree, distance_matrix, fig

# Example usage
protein_sequences = [
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQF",
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEV",
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALP",
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPD"
]
ids = ["Human", "Mouse", "Rat", "Dog"]

tree, distances, fig = phylogenetic_analysis(protein_sequences, ids, method='upgma')
plt.show()

# Save tree in different formats
from Bio import Phylo
Phylo.write(tree, "phylogenetic_tree.nwk", "newick")
Phylo.write(tree, "phylogenetic_tree.xml", "phyloxml")
```

## Machine Learning Workflows

### Feature Extraction for ML

```python
import kalign
import numpy as np
import pandas as pd
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

def extract_alignment_features(sequences, ids=None):
    """Extract features from aligned sequences for ML applications."""
    
    if ids is None:
        ids = [f"seq{i}" for i in range(len(sequences))]
    
    print("🤖 Extracting features for machine learning...")
    
    # Align sequences
    aligned = kalign.align(sequences)
    
    # Feature 1: Sequence composition
    features_df = pd.DataFrame()
    
    for i, (seq_id, seq) in enumerate(zip(ids, aligned)):
        # Basic composition features
        seq_no_gaps = seq.replace('-', '')
        total_length = len(seq_no_gaps)
        
        features = {
            'sequence_id': seq_id,
            'original_length': total_length,
            'alignment_length': len(seq),
            'gap_count': seq.count('-'),
            'gap_fraction': seq.count('-') / len(seq)
        }
        
        # Character frequency features
        unique_chars = set(seq_no_gaps)
        for char in unique_chars:
            features[f'freq_{char}'] = seq_no_gaps.count(char) / total_length if total_length > 0 else 0
        
        # Positional features
        features['gaps_start'] = len(seq) - len(seq.lstrip('-'))
        features['gaps_end'] = len(seq) - len(seq.rstrip('-'))
        features['internal_gaps'] = features['gap_count'] - features['gaps_start'] - features['gaps_end']
        
        features_df = pd.concat([features_df, pd.DataFrame([features])], ignore_index=True)
    
    # Feature 2: Pairwise similarity features
    identity_matrix = kalign.utils.pairwise_identity_matrix(aligned)
    
    for i, seq_id in enumerate(ids):
        # Average similarity to other sequences
        other_similarities = [identity_matrix[i, j] for j in range(len(ids)) if i != j]
        features_df.loc[i, 'avg_similarity'] = np.mean(other_similarities)
        features_df.loc[i, 'max_similarity'] = max(other_similarities) if other_similarities else 0
        features_df.loc[i, 'min_similarity'] = min(other_similarities) if other_similarities else 0
        features_df.loc[i, 'similarity_std'] = np.std(other_similarities)
    
    # Feature 3: Conservation-based features
    conservation_scores = []
    alignment_array = kalign.utils.to_array(aligned)
    
    for col in range(alignment_array.shape[1]):
        column = alignment_array[:, col]
        non_gap = column[column != '-']
        if len(non_gap) > 0:
            unique_chars = len(set(non_gap))
            conservation = 1.0 / unique_chars if unique_chars > 0 else 0
        else:
            conservation = 0
        conservation_scores.append(conservation)
    
    # Add conservation-based features for each sequence
    for i, seq_id in enumerate(ids):
        seq_positions = []
        for j, char in enumerate(aligned[i]):
            if char != '-':
                seq_positions.append(conservation_scores[j])
        
        if seq_positions:
            features_df.loc[i, 'avg_position_conservation'] = np.mean(seq_positions)
            features_df.loc[i, 'conserved_positions_fraction'] = sum(1 for x in seq_positions if x > 0.8) / len(seq_positions)
        else:
            features_df.loc[i, 'avg_position_conservation'] = 0
            features_df.loc[i, 'conserved_positions_fraction'] = 0
    
    # Fill NaN values
    features_df = features_df.fillna(0)
    
    print(f"✅ Extracted {features_df.shape[1]-1} features for {len(sequences)} sequences")
    
    return features_df, aligned, identity_matrix

def ml_sequence_clustering(features_df, n_clusters=3):
    """Perform clustering analysis on sequence features."""
    
    print(f"🎯 Performing K-means clustering (k={n_clusters})...")
    
    # Prepare features for clustering
    feature_columns = [col for col in features_df.columns if col != 'sequence_id']
    X = features_df[feature_columns].values
    
    # Standardize features
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # Apply PCA for dimensionality reduction
    pca = PCA(n_components=min(5, X_scaled.shape[1]))
    X_pca = pca.fit_transform(X_scaled)
    
    # Clustering
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    clusters = kmeans.fit_predict(X_scaled)
    
    # Add cluster labels to dataframe
    features_df['cluster'] = clusters
    
    # Visualization
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    # PCA plot
    if X_pca.shape[1] >= 2:
        scatter = axes[0].scatter(X_pca[:, 0], X_pca[:, 1], c=clusters, cmap='tab10')
        axes[0].set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)')
        axes[0].set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)')
        axes[0].set_title('PCA with Clusters')
        
        # Add sequence labels
        for i, seq_id in enumerate(features_df['sequence_id']):
            axes[0].annotate(seq_id, (X_pca[i, 0], X_pca[i, 1]), 
                           xytext=(5, 5), textcoords='offset points', fontsize=8)
    
    # Feature importance (PCA loadings)
    feature_importance = np.abs(pca.components_[0]) if pca.components_.shape[0] > 0 else np.zeros(len(feature_columns))
    feature_names = feature_columns
    
    # Sort features by importance
    importance_order = np.argsort(feature_importance)[::-1]
    top_features = [feature_names[i] for i in importance_order[:10]]
    top_importance = [feature_importance[i] for i in importance_order[:10]]
    
    axes[1].barh(range(len(top_features)), top_importance)
    axes[1].set_yticks(range(len(top_features)))
    axes[1].set_yticklabels(top_features)
    axes[1].set_xlabel('Importance (|PC1 loading|)')
    axes[1].set_title('Top 10 Feature Importance')
    
    # Cluster composition
    cluster_counts = pd.Series(clusters).value_counts().sort_index()
    axes[2].bar(cluster_counts.index, cluster_counts.values)
    axes[2].set_xlabel('Cluster')
    axes[2].set_ylabel('Number of Sequences')
    axes[2].set_title('Cluster Sizes')
    
    # Add cluster labels
    for i, (cluster, count) in enumerate(cluster_counts.items()):
        axes[2].text(cluster, count + 0.1, str(count), ha='center', va='bottom')
    
    plt.tight_layout()
    
    # Print cluster summary
    print(f"\n📊 Clustering Results:")
    for cluster in range(n_clusters):
        cluster_sequences = features_df[features_df['cluster'] == cluster]['sequence_id'].tolist()
        print(f"  Cluster {cluster}: {cluster_sequences}")
    
    # Feature statistics by cluster
    print(f"\n📈 Cluster Characteristics:")
    cluster_stats = features_df.groupby('cluster')[feature_columns].mean()
    print(cluster_stats.round(3))
    
    return features_df, X_pca, clusters, pca, scaler, fig

# Example usage
sequences = [
    "ATCGATCGATCGATCGATCG",  # DNA-like
    "ATCGATCGTCGATCGATCG",   # DNA-like with variation
    "ATCGATCGATCATCGATCG",   # DNA-like with variation
    "ACDEFGHIKLMNPQRSTVWY",  # Protein-like
    "ACDEFGHIKLMNPQRSTVWF",  # Protein-like with variation
    "UUUUUUUUUUUUUUUUUUUU"   # RNA-like, highly repetitive
]

ids = ["DNA_1", "DNA_2", "DNA_3", "Protein_1", "Protein_2", "RNA_1"]

# Extract features
features_df, aligned, identity_matrix = extract_alignment_features(sequences, ids)

# Perform clustering
features_df, pca_data, clusters, pca_model, scaler, fig = ml_sequence_clustering(features_df, n_clusters=3)

plt.show()

# Display feature matrix
print("\n🔬 Feature Matrix:")
print(features_df.round(3))
```

## Jupyter Notebook Integration

### Notebook-Friendly Utilities

```python
import kalign
import matplotlib.pyplot as plt
from IPython.display import display, HTML, Markdown
import pandas as pd

class KalignNotebook:
    """Jupyter notebook utilities for Kalign."""
    
    def __init__(self):
        self.last_alignment = None
        self.last_stats = None
    
    def align_and_display(self, sequences, ids=None, seq_type="auto", **kwargs):
        """Align sequences and display results in notebook-friendly format."""
        
        if ids is None:
            ids = [f"seq{i+1}" for i in range(len(sequences))]
        
        # Perform alignment
        aligned = kalign.align(sequences, seq_type=seq_type, **kwargs)
        stats = kalign.utils.alignment_stats(aligned)
        
        # Store for later use
        self.last_alignment = aligned
        self.last_stats = stats
        
        # Display markdown summary
        summary = f"""
## Alignment Results

**Sequences:** {stats['n_sequences']}  
**Alignment Length:** {stats['length']}  
**Gap Fraction:** {stats['gap_fraction']:.2%}  
**Conservation:** {stats['conservation']:.2%}  
**Average Identity:** {stats['identity']:.2%}
        """
        display(Markdown(summary))
        
        # Display alignment as HTML table
        self.display_alignment_table(aligned, ids)
        
        return aligned, stats
    
    def display_alignment_table(self, aligned, ids, max_width=80):
        """Display alignment as formatted HTML table."""
        
        # Split long alignments into chunks
        alignment_length = len(aligned[0])
        chunks = []
        
        for start in range(0, alignment_length, max_width):
            end = min(start + max_width, alignment_length)
            chunk_data = []
            
            # Header row with positions
            positions = "".join([str((i+1) % 10) for i in range(start, end)])
            chunk_data.append(f"<tr><th>Pos</th><td><tt>{positions}</tt></td></tr>")
            
            # Sequence rows
            for seq_id, seq in zip(ids, aligned):
                chunk = seq[start:end]
                # Color gaps differently
                colored_chunk = chunk.replace('-', '<span style="color: #ccc;">-</span>')
                chunk_data.append(f"<tr><th>{seq_id}</th><td><tt>{colored_chunk}</tt></td></tr>")
            
            chunks.append(f"""
            <table style="font-family: monospace; border-collapse: collapse; margin: 10px 0;">
                <caption><strong>Positions {start+1}-{end}</strong></caption>
                {''.join(chunk_data)}
            </table>
            """)
        
        # Display all chunks
        html_content = f"""
        <div style="max-height: 400px; overflow-y: auto; border: 1px solid #ddd; padding: 10px;">
            {''.join(chunks)}
        </div>
        """
        display(HTML(html_content))
    
    def plot_conservation(self, figsize=(12, 4)):
        """Plot conservation scores for the last alignment."""
        
        if self.last_alignment is None:
            print("No alignment available. Run align_and_display() first.")
            return
        
        # Calculate conservation
        alignment_array = kalign.utils.to_array(self.last_alignment)
        conservation_scores = []
        
        for col in range(alignment_array.shape[1]):
            column = alignment_array[:, col]
            non_gap = column[column != '-']
            if len(non_gap) > 0:
                unique_chars = len(set(non_gap))
                conservation = 1.0 / unique_chars if unique_chars > 0 else 0
            else:
                conservation = 0
            conservation_scores.append(conservation)
        
        # Plot
        fig, ax = plt.subplots(figsize=figsize)
        ax.plot(range(1, len(conservation_scores) + 1), conservation_scores, 'b-', linewidth=2)
        ax.fill_between(range(1, len(conservation_scores) + 1), conservation_scores, alpha=0.3)
        ax.set_xlabel('Position')
        ax.set_ylabel('Conservation Score')
        ax.set_title('Position-wise Conservation')
        ax.grid(True, alpha=0.3)
        
        # Highlight highly conserved regions
        threshold = 0.8
        for i, score in enumerate(conservation_scores):
            if score >= threshold:
                ax.axvline(i + 1, color='red', alpha=0.5, linewidth=0.5)
        
        plt.tight_layout()
        plt.show()
        
        return fig, conservation_scores
    
    def export_results(self, filename_prefix="kalign_results"):
        """Export alignment results in multiple formats."""
        
        if self.last_alignment is None:
            print("No alignment available. Run align_and_display() first.")
            return
        
        # Generate IDs if not available
        ids = [f"seq{i+1}" for i in range(len(self.last_alignment))]
        
        # Export in different formats
        formats = {
            'fasta': kalign.io.write_fasta,
            'clustal': kalign.io.write_clustal,
            'stockholm': kalign.io.write_stockholm
        }
        
        exported_files = []
        for fmt, write_func in formats.items():
            filename = f"{filename_prefix}.{fmt}"
            try:
                write_func(self.last_alignment, filename, ids=ids)
                exported_files.append(filename)
            except Exception as e:
                print(f"Error exporting {fmt}: {e}")
        
        display(Markdown(f"**Exported files:** {', '.join(exported_files)}"))
        return exported_files

# Example notebook usage
def notebook_example():
    """Example usage in Jupyter notebook."""
    
    # Initialize notebook utilities
    kn = KalignNotebook()
    
    # Example sequences
    sequences = [
        "ATCGATCGATCGATCGATCG",
        "ATCGATCGTCGATCGATCG",
        "ATCGATCGATCATCGATCG",
        "ATCGATCGAGCTCGATCG"
    ]
    
    ids = ["Human", "Mouse", "Rat", "Dog"]
    
    # Align and display
    aligned, stats = kn.align_and_display(sequences, ids=ids, seq_type="dna", n_threads=2)
    
    # Plot conservation
    fig, conservation = kn.plot_conservation()
    
    # Export results
    files = kn.export_results("example_alignment")
    
    return kn, aligned, stats

# Run example (uncomment in notebook)
# kn, aligned, stats = notebook_example()
```

### Interactive Widgets

```python
import kalign
import ipywidgets as widgets
from IPython.display import display, clear_output
import matplotlib.pyplot as plt

def create_alignment_widget():
    """Create interactive widget for sequence alignment."""
    
    # Widget components
    sequence_input = widgets.Textarea(
        value="ATCGATCGATCG\nATCGTCGATCG\nATCGATCATCG",
        placeholder="Enter sequences (one per line)",
        description="Sequences:",
        layout=widgets.Layout(width='100%', height='100px')
    )
    
    seq_type_dropdown = widgets.Dropdown(
        options=['auto', 'dna', 'rna', 'protein', 'divergent'],
        value='auto',
        description='Seq Type:'
    )
    
    threads_slider = widgets.IntSlider(
        value=1,
        min=1,
        max=8,
        step=1,
        description='Threads:'
    )
    
    format_dropdown = widgets.Dropdown(
        options=['plain', 'biopython', 'skbio'],
        value='plain',
        description='Format:'
    )
    
    align_button = widgets.Button(
        description="Align Sequences",
        button_style='primary'
    )
    
    output_area = widgets.Output()
    
    def on_align_clicked(b):
        with output_area:
            clear_output()
            
            try:
                # Parse sequences
                sequences = [seq.strip() for seq in sequence_input.value.split('\n') if seq.strip()]
                
                if len(sequences) < 2:
                    print("❌ Please provide at least 2 sequences")
                    return
                
                print(f"🧬 Aligning {len(sequences)} sequences...")
                
                # Perform alignment
                aligned = kalign.align(
                    sequences,
                    seq_type=seq_type_dropdown.value,
                    n_threads=threads_slider.value,
                    fmt=format_dropdown.value
                )
                
                # Display results
                if format_dropdown.value == 'plain':
                    print("✅ Alignment complete!")
                    print("\n📋 Aligned sequences:")
                    for i, seq in enumerate(aligned):
                        print(f"Seq {i+1}: {seq}")
                    
                    # Calculate and display statistics
                    stats = kalign.utils.alignment_stats(aligned)
                    print(f"\n📊 Statistics:")
                    print(f"  Length: {stats['length']}")
                    print(f"  Gap fraction: {stats['gap_fraction']:.2%}")
                    print(f"  Conservation: {stats['conservation']:.2%}")
                    print(f"  Average identity: {stats['identity']:.2%}")
                
                else:
                    print(f"✅ Alignment complete! (Format: {format_dropdown.value})")
                    print(f"Object type: {type(aligned)}")
                    
                    if hasattr(aligned, 'get_alignment_length'):
                        print(f"Alignment length: {aligned.get_alignment_length()}")
                    elif hasattr(aligned, '__len__'):
                        print(f"Number of sequences: {len(aligned)}")
                
            except Exception as e:
                print(f"❌ Error: {e}")
    
    align_button.on_click(on_align_clicked)
    
    # Layout
    controls = widgets.VBox([
        widgets.HBox([seq_type_dropdown, threads_slider, format_dropdown]),
        align_button
    ])
    
    widget = widgets.VBox([
        sequence_input,
        controls,
        output_area
    ])
    
    return widget

# Usage in notebook:
# alignment_widget = create_alignment_widget()
# display(alignment_widget)
```

This comprehensive ecosystem integration guide demonstrates how Kalign seamlessly integrates with the entire bioinformatics ecosystem. The documentation provides practical examples for every major use case, from basic Biopython workflows to advanced machine learning applications and interactive Jupyter notebook tools.

<function_calls>
<invoke name="TodoWrite">
<parameter name="todos">[{"content": "Create comprehensive Python API documentation (docs/python-api.md)", "status": "completed", "priority": "high", "id": "1"}, {"content": "Create quick start guide with examples (docs/python-quickstart.md)", "status": "completed", "priority": "high", "id": "2"}, {"content": "Create ecosystem integration guide (docs/python-ecosystem.md)", "status": "completed", "priority": "high", "id": "3"}, {"content": "Create performance tuning guide (docs/python-performance.md)", "status": "in_progress", "priority": "medium", "id": "4"}, {"content": "Create troubleshooting guide (docs/python-troubleshooting.md)", "status": "pending", "priority": "medium", "id": "5"}, {"content": "Update main README with Python section", "status": "pending", "priority": "high", "id": "6"}, {"content": "Create example scripts directory", "status": "pending", "priority": "medium", "id": "7"}]