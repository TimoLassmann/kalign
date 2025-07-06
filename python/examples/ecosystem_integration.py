#!/usr/bin/env python3
"""
Ecosystem Integration Examples

Demonstrates integration with Biopython, scikit-bio, pandas, and matplotlib.
"""

import kalign
import sys

def check_dependencies():
    """Check which optional dependencies are available."""
    available = {}
    
    try:
        import Bio
        available['biopython'] = Bio.__version__
    except ImportError:
        available['biopython'] = None
    
    try:
        import skbio
        available['skbio'] = skbio.__version__
    except ImportError:
        available['skbio'] = None
    
    try:
        import pandas
        available['pandas'] = pandas.__version__
    except ImportError:
        available['pandas'] = None
    
    try:
        import matplotlib
        available['matplotlib'] = matplotlib.__version__
    except ImportError:
        available['matplotlib'] = None
    
    return available


def example_biopython():
    """Example: Biopython integration."""
    print("=" * 50)
    print("Example: Biopython Integration")
    print("=" * 50)
    
    try:
        from Bio import AlignIO
        from Bio.Align import AlignInfo
        
        # Sample sequences
        sequences = [
            "ATCGATCGATCGATCGATCG",
            "ATCGATCGTCGATCGATCG",
            "ATCGATCGATCATCGATCG",
            "ATCGATCGAGCTCGATCG"
        ]
        ids = ["human", "mouse", "rat", "dog"]
        
        print("Aligning sequences with Biopython output...")
        
        # Align and get Biopython object
        alignment = kalign.align(sequences, fmt="biopython", ids=ids)
        
        print(f"Alignment type: {type(alignment)}")
        print(f"Alignment length: {alignment.get_alignment_length()}")
        print(f"Number of sequences: {len(alignment)}")
        
        # Use Biopython's modern consensus tools
        try:
            # Try modern approach first (Bio 1.81+)
            from Bio.motifs import Motif
            motif = Motif(alignment=alignment.alignment)
            consensus = motif.consensus
        except (ImportError, AttributeError):
            # Fallback to older method with warning suppression
            import warnings
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=DeprecationWarning)
                summary = AlignInfo.SummaryInfo(alignment)
                consensus = summary.gap_consensus(threshold=0.7, ambiguous='N')
        
        print(f"\nConsensus sequence: {consensus}")
        
        # Display alignment
        print("\nAlignment:")
        for record in alignment:
            print(f"{record.id:8}: {record.seq}")
        
        # Export to different formats
        print("\nExporting to formats...")
        AlignIO.write(alignment, "output.clustal", "clustal")
        AlignIO.write(alignment, "output.phylip", "phylip")
        print("Exported to output.clustal and output.phylip")
        
        # Clean up
        import os
        for f in ["output.clustal", "output.phylip"]:
            if os.path.exists(f):
                os.remove(f)
        
    except ImportError:
        print("❌ Biopython not available. Install with: pip install kalign[biopython]")


def example_skbio():
    """Example: scikit-bio integration."""
    print("\n" + "=" * 50)
    print("Example: scikit-bio Integration")
    print("=" * 50)
    
    try:
        import skbio
        
        # Sample sequences
        sequences = [
            "ATCGATCGATCGATCGATCG",
            "ATCGATCGTCGATCGATCG",
            "ATCGATCGATCATCGATCG",
            "ATCGATCGAGCTCGATCG"
        ]
        
        print("Aligning sequences with scikit-bio output...")
        
        # Align and get scikit-bio object
        alignment = kalign.align(sequences, fmt="skbio")
        
        print(f"Alignment type: {type(alignment)}")
        print(f"Number of sequences: {len(alignment)}")
        
        # Use scikit-bio's analysis tools
        consensus = alignment.consensus()
        print(f"\nConsensus sequence: {consensus}")
        
        # Calculate conservation
        conservation = alignment.conservation()
        print(f"Conservation scores (first 10): {conservation[:10]}")
        
        # Calculate gap frequencies
        gap_frequencies = []
        for i in range(len(alignment[0])):
            column = [str(seq)[i] for seq in alignment]
            gap_freq = column.count('-') / len(column)
            gap_frequencies.append(gap_freq)
        
        print(f"Gap frequencies (first 10): {[f'{x:.2f}' for x in gap_frequencies[:10]]}")
        
        # Export alignment
        alignment.write("output.fasta", format="fasta")
        print("Exported to output.fasta")
        
        # Clean up
        import os
        if os.path.exists("output.fasta"):
            os.remove("output.fasta")
        
    except ImportError:
        print("❌ scikit-bio not available. Install with: pip install kalign[skbio]")


def example_pandas():
    """Example: Pandas integration for analysis."""
    print("\n" + "=" * 50)
    print("Example: Pandas Integration")
    print("=" * 50)
    
    try:
        import pandas as pd
        import numpy as np
        
        # Create multiple alignments for comparison
        sequence_groups = {
            "Group_A": [
                "ATCGATCGATCGATCGATCG",
                "ATCGATCGTCGATCGATCG",
                "ATCGATCGATCATCGATCG"
            ],
            "Group_B": [
                "GCTAGCTAGCTAGCTAGCTA",
                "GCTAGCTACTAGCTAGCTA",
                "GCTAGCTAGCTAACTAGCTA"
            ],
            "Group_C": [
                "TTTTTTTTTTTTTTTTTTTT",
                "TTTTTTCTTTTTTTTTTTTT",
                "TTTTTTTTTTCTTTTTTTTT"
            ]
        }
        
        print("Comparing multiple sequence groups...")
        
        # Analyze each group
        results = []
        for group_name, sequences in sequence_groups.items():
            aligned = kalign.align(sequences)
            stats = kalign.utils.alignment_stats(aligned)
            
            # Calculate additional metrics
            identity_matrix = kalign.utils.pairwise_identity_matrix(aligned)
            mean_identity = np.mean(identity_matrix[np.triu_indices_from(identity_matrix, k=1)])
            
            results.append({
                'group': group_name,
                'n_sequences': len(sequences),
                'alignment_length': stats['length'],
                'gap_fraction': stats['gap_fraction'],
                'conservation': stats['conservation'],
                'mean_pairwise_identity': mean_identity
            })
        
        # Create DataFrame
        df = pd.DataFrame(results)
        
        print("\nComparison results:")
        print(df.round(3))
        
        # Find best and worst groups
        best_conservation = df.loc[df['conservation'].idxmax(), 'group']
        best_identity = df.loc[df['mean_pairwise_identity'].idxmax(), 'group']
        most_gaps = df.loc[df['gap_fraction'].idxmax(), 'group']
        
        print(f"\nSummary:")
        print(f"  Most conserved: {best_conservation}")
        print(f"  Highest identity: {best_identity}")
        print(f"  Most gaps: {most_gaps}")
        
    except ImportError:
        print("❌ Pandas not available. Install with: pip install pandas")


def example_visualization():
    """Example: Matplotlib visualization."""
    print("\n" + "=" * 50)
    print("Example: Matplotlib Visualization")
    print("=" * 50)
    
    try:
        import matplotlib.pyplot as plt
        import numpy as np
        
        # Create test sequences
        sequences = [
            "ATCGATCGATCGATCGATCGATCGATCGATCGATCG",
            "ATCGATCGTCGATCGATCGATCGATCGTCGATCG",
            "ATCGATCGATCATCGATCGATCGATCGATCATCG",
            "ATCGATCGAGCTCGATCGATCGATCGAGCTCG",
            "ATCGATCGATCGATCAACGATCGATCGATCAACG"
        ]
        
        print("Creating visualization of alignment...")
        
        # Align sequences
        aligned = kalign.align(sequences)
        
        # Convert to numeric matrix for visualization
        char_to_num = {'-': 0, 'A': 1, 'T': 2, 'C': 3, 'G': 4}
        matrix = np.array([[char_to_num.get(char, 0) for char in seq] for seq in aligned])
        
        # Create visualization
        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        fig.suptitle('Kalign Alignment Analysis')
        
        # 1. Alignment heatmap
        im1 = axes[0, 0].imshow(matrix, aspect='auto', cmap='viridis')
        axes[0, 0].set_title('Alignment Heatmap')
        axes[0, 0].set_xlabel('Position')
        axes[0, 0].set_ylabel('Sequence')
        
        # 2. Position-wise conservation
        conservation = []
        for j in range(matrix.shape[1]):
            column = matrix[:, j]
            non_gap = column[column != 0]
            if len(non_gap) > 0:
                unique_vals = len(set(non_gap))
                conservation.append(1.0 / unique_vals)
            else:
                conservation.append(0)
        
        axes[0, 1].plot(conservation)
        axes[0, 1].set_title('Conservation Score by Position')
        axes[0, 1].set_xlabel('Position')
        axes[0, 1].set_ylabel('Conservation')
        
        # 3. Gap analysis
        gap_counts = [(matrix[:, j] == 0).sum() for j in range(matrix.shape[1])]
        axes[1, 0].bar(range(len(gap_counts)), gap_counts)
        axes[1, 0].set_title('Gaps per Position')
        axes[1, 0].set_xlabel('Position')
        axes[1, 0].set_ylabel('Number of Gaps')
        
        # 4. Pairwise identity matrix
        identity_matrix = kalign.utils.pairwise_identity_matrix(aligned)
        im4 = axes[1, 1].imshow(identity_matrix, cmap='RdYlBu_r', vmin=0, vmax=1)
        axes[1, 1].set_title('Pairwise Identity Matrix')
        axes[1, 1].set_xlabel('Sequence')
        axes[1, 1].set_ylabel('Sequence')
        
        plt.tight_layout()
        plt.savefig('alignment_analysis.png', dpi=150, bbox_inches='tight')
        print("Saved visualization to: alignment_analysis.png")
        
        # Show plot if in interactive mode
        if hasattr(plt.get_backend(), 'show'):
            plt.show()
        
        plt.close()
        
        # Clean up
        import os
        if os.path.exists('alignment_analysis.png'):
            os.remove('alignment_analysis.png')
        
    except ImportError:
        print("❌ Matplotlib not available. Install with: pip install matplotlib")


def example_comprehensive_workflow():
    """Example: Comprehensive analysis workflow."""
    print("\n" + "=" * 50)
    print("Example: Comprehensive Workflow")
    print("=" * 50)
    
    # Sample protein sequences
    protein_sequences = [
        "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ",
        "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQF",
        "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALP",
        "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPD"
    ]
    
    ids = ["Human", "Mouse", "Rat", "Dog"]
    
    print(f"Analyzing {len(protein_sequences)} protein sequences...")
    
    # Step 1: Align sequences
    aligned = kalign.align(protein_sequences, seq_type="protein")
    
    # Step 2: Calculate basic statistics
    stats = kalign.utils.alignment_stats(aligned)
    print(f"\nBasic Statistics:")
    print(f"  Alignment length: {stats['length']}")
    print(f"  Gap fraction: {stats['gap_fraction']:.2%}")
    print(f"  Conservation: {stats['conservation']:.2%}")
    print(f"  Average identity: {stats['identity']:.2%}")
    
    # Step 3: Generate consensus
    consensus = kalign.utils.consensus_sequence(aligned, threshold=0.75)
    print(f"\nConsensus (75%): {consensus}")
    
    # Step 4: Pairwise analysis
    identity_matrix = kalign.utils.pairwise_identity_matrix(aligned)
    print(f"\nPairwise identities:")
    for i in range(len(ids)):
        for j in range(i+1, len(ids)):
            identity = identity_matrix[i, j]
            print(f"  {ids[i]} vs {ids[j]}: {identity:.2%}")
    
    # Step 5: Find conserved regions
    import numpy as np
    array = kalign.utils.to_array(aligned)
    conserved_positions = []
    
    for col in range(array.shape[1]):
        column = array[:, col]
        non_gap = column[column != '-']
        if len(non_gap) > 0 and len(set(non_gap)) == 1:
            conserved_positions.append(col + 1)  # 1-based position
    
    print(f"\nConserved positions: {conserved_positions[:10]}...")
    print(f"Total conserved positions: {len(conserved_positions)}")
    
    # Step 6: Export results
    print(f"\nExporting results...")
    
    # Save alignment
    with open("protein_alignment.fasta", 'w') as f:
        for seq_id, seq in zip(ids, aligned):
            f.write(f">{seq_id}\n{seq}\n")
    
    # Save analysis report
    with open("analysis_report.txt", 'w') as f:
        f.write("Protein Sequence Analysis Report\n")
        f.write("=" * 40 + "\n\n")
        f.write(f"Sequences analyzed: {len(protein_sequences)}\n")
        f.write(f"Alignment length: {stats['length']}\n")
        f.write(f"Gap fraction: {stats['gap_fraction']:.2%}\n")
        f.write(f"Conservation: {stats['conservation']:.2%}\n")
        f.write(f"Average identity: {stats['identity']:.2%}\n\n")
        f.write(f"Consensus sequence:\n{consensus}\n\n")
        f.write(f"Conserved positions ({len(conserved_positions)} total):\n")
        f.write(f"{conserved_positions}\n")
    
    print("Saved: protein_alignment.fasta, analysis_report.txt")
    
    # Clean up
    import os
    for f in ["protein_alignment.fasta", "analysis_report.txt"]:
        if os.path.exists(f):
            os.remove(f)
    
    print("✅ Comprehensive analysis complete!")


def main():
    """Run ecosystem integration examples."""
    print("🌐 Kalign Ecosystem Integration Examples")
    print("Demonstrates integration with popular bioinformatics tools")
    
    # Check available dependencies
    dependencies = check_dependencies()
    print("\nAvailable dependencies:")
    for dep, version in dependencies.items():
        status = f"✅ {version}" if version else "❌ Not installed"
        print(f"  {dep}: {status}")
    
    try:
        # Run available examples
        if dependencies['biopython']:
            try:
                example_biopython()
            except Exception as e:
                print(f"\n❌ Biopython example failed: {e}")
        
        if dependencies['skbio']:
            try:
                example_skbio()
            except Exception as e:
                print(f"\n❌ scikit-bio example failed: {e}")
        
        if dependencies['pandas']:
            try:
                example_pandas()
            except Exception as e:
                print(f"\n❌ Pandas example failed: {e}")
        
        if dependencies['matplotlib']:
            try:
                example_visualization()
            except Exception as e:
                print(f"\n❌ Matplotlib example failed: {e}")
        
        # Always run comprehensive workflow
        try:
            example_comprehensive_workflow()
        except Exception as e:
            print(f"\n❌ Comprehensive workflow failed: {e}")
        
        print("\n" + "=" * 50)
        print("✅ Ecosystem integration examples completed!")
        print("=" * 50)
        
        # Installation suggestions
        missing_deps = [dep for dep, version in dependencies.items() if not version]
        if missing_deps:
            print(f"\n💡 To unlock more features, install missing dependencies:")
            for dep in missing_deps:
                if dep == 'biopython':
                    print(f"   pip install kalign[biopython]")
                elif dep == 'skbio':
                    print(f"   pip install kalign[skbio]")
                else:
                    print(f"   pip install {dep}")
        
    except Exception as e:
        print(f"\n❌ Critical error running examples: {e}")
        print("Check your installations and try again.")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()