#!/usr/bin/env python3
"""
Basic Kalign Usage Examples

This script demonstrates the most common use cases for Kalign Python package.
"""

import kalign

def example_1_basic_alignment():
    """Example 1: Basic sequence alignment."""
    print("=" * 50)
    print("Example 1: Basic Sequence Alignment")
    print("=" * 50)
    
    # DNA sequences
    dna_sequences = [
        "ATCGATCGATCGATCG",
        "ATCGATCGTCGATCG",
        "ATCGATCGATCATCG",
        "ATCGATCGAGATCG"
    ]
    
    print("Input sequences:")
    for i, seq in enumerate(dna_sequences):
        print(f"  Seq {i+1}: {seq}")
    
    # Align sequences
    aligned = kalign.align(dna_sequences, seq_type="dna")
    
    print("\nAligned sequences:")
    for i, seq in enumerate(aligned):
        print(f"  Seq {i+1}: {seq}")
    
    # Calculate statistics
    stats = kalign.utils.alignment_stats(aligned)
    print(f"\nAlignment statistics:")
    print(f"  Length: {stats['length']}")
    print(f"  Gap fraction: {stats['gap_fraction']:.2%}")
    print(f"  Conservation: {stats['conservation']:.2%}")
    print(f"  Average identity: {stats['identity']:.2%}")


def example_2_protein_alignment():
    """Example 2: Protein sequence alignment."""
    print("\n" + "=" * 50)
    print("Example 2: Protein Sequence Alignment")
    print("=" * 50)
    
    # Protein sequences
    protein_sequences = [
        "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK",
        "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAV",
        "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSG"
    ]
    
    print("Input protein sequences:")
    for i, seq in enumerate(protein_sequences):
        print(f"  Protein {i+1}: {seq[:30]}...")
    
    # Align with custom gap penalties
    aligned = kalign.align(
        protein_sequences,
        seq_type="protein",
        gap_open=-10.0,
        gap_extend=-1.0
    )
    
    print("\nAligned sequences:")
    for i, seq in enumerate(aligned):
        print(f"  Protein {i+1}: {seq}")
    
    # Generate consensus
    consensus = kalign.utils.consensus_sequence(aligned, threshold=0.7)
    print(f"\nConsensus (70%): {consensus}")


def example_3_threading():
    """Example 3: Multi-threaded alignment."""
    print("\n" + "=" * 50)
    print("Example 3: Multi-threaded Alignment")
    print("=" * 50)
    
    # Create larger dataset for threading demo
    base_seq = "ATCGATCGATCG" * 50  # 600 bp sequence
    sequences = []
    
    import random
    random.seed(42)  # For reproducible results
    
    for i in range(20):
        seq = list(base_seq)
        # Add 2% mutations
        n_mutations = len(seq) // 50
        for _ in range(n_mutations):
            pos = random.randint(0, len(seq) - 1)
            seq[pos] = random.choice(['A', 'T', 'C', 'G'])
        sequences.append(''.join(seq))
    
    print(f"Created {len(sequences)} sequences of {len(sequences[0])} bp each")
    
    # Set global thread count
    kalign.set_num_threads(4)
    print(f"Using {kalign.get_num_threads()} threads")
    
    import time
    start_time = time.time()
    aligned = kalign.align(sequences, seq_type="dna")
    end_time = time.time()
    
    print(f"Alignment completed in {end_time - start_time:.2f} seconds")
    print(f"Aligned {len(aligned)} sequences to length {len(aligned[0])}")


def example_4_file_operations():
    """Example 4: File input/output operations."""
    print("\n" + "=" * 50)
    print("Example 4: File Input/Output")
    print("=" * 50)
    
    # Create sample FASTA file
    sample_file = "sample_sequences.fasta"
    sequences = [
        "ATCGATCGATCGATCG",
        "ATCGATCGTCGATCG",
        "ATCGATCGATCATCG"
    ]
    ids = ["seq1", "seq2", "seq3"]
    
    # Write sample file
    with open(sample_file, 'w') as f:
        for seq_id, seq in zip(ids, sequences):
            f.write(f">{seq_id}\n{seq}\n")
    
    print(f"Created sample file: {sample_file}")
    
    # Read and align from file
    aligned = kalign.align_from_file(sample_file, seq_type="dna")
    print(f"Aligned {len(aligned)} sequences from file")
    
    # Write aligned sequences
    output_file = "aligned_sequences.fasta"
    kalign.write_alignment(aligned, output_file, format="fasta")
    print(f"Wrote aligned sequences to: {output_file}")
    
    # Clean up
    import os
    os.remove(sample_file)
    os.remove(output_file)
    print("Cleaned up temporary files")


def example_5_analysis():
    """Example 5: Alignment analysis."""
    print("\n" + "=" * 50)
    print("Example 5: Alignment Analysis")
    print("=" * 50)
    
    # Sequences with varying similarity
    sequences = [
        "ATCGATCGATCGATCGATCG",
        "ATCGATCGTCGATCGATCG",  # 1 difference
        "ATCGATCGATCATCGATCG",  # 1 difference  
        "ATCGATCGAGCTCGATCG",   # 2 differences
        "GCTAGCTAGCTAGCTAGCTA"  # Very different
    ]
    
    # Align sequences
    aligned = kalign.align(sequences)
    
    print("Aligned sequences:")
    for i, seq in enumerate(aligned):
        print(f"  Seq {i+1}: {seq}")
    
    # Pairwise identity matrix
    identity_matrix = kalign.utils.pairwise_identity_matrix(aligned)
    
    print("\nPairwise identity matrix:")
    print("     ", end="")
    for i in range(len(sequences)):
        print(f"Seq{i+1:1}  ", end="")
    print()
    
    for i in range(len(sequences)):
        print(f"Seq{i+1}: ", end="")
        for j in range(len(sequences)):
            print(f"{identity_matrix[i, j]:.2f} ", end="")
        print()
    
    # Find most and least similar sequences
    import numpy as np
    
    # Mask diagonal for finding max/min
    masked_matrix = identity_matrix + np.eye(len(sequences)) * -1
    max_pos = np.unravel_index(np.argmax(masked_matrix), masked_matrix.shape)
    min_pos = np.unravel_index(np.argmin(identity_matrix), identity_matrix.shape)
    
    print(f"\nMost similar: Seq{max_pos[0]+1} vs Seq{max_pos[1]+1} "
          f"({identity_matrix[max_pos]:.2%} identity)")
    print(f"Least similar: Seq{min_pos[0]+1} vs Seq{min_pos[1]+1} "
          f"({identity_matrix[min_pos]:.2%} identity)")


def main():
    """Run all examples."""
    print("🧬 Kalign Python Examples")
    print("This script demonstrates common usage patterns")
    
    try:
        example_1_basic_alignment()
        example_2_protein_alignment()
        example_3_threading()
        example_4_file_operations()
        example_5_analysis()
        
        print("\n" + "=" * 50)
        print("✅ All examples completed successfully!")
        print("=" * 50)
        
    except Exception as e:
        print(f"\n❌ Error running examples: {e}")
        print("Check your Kalign installation and try again.")


if __name__ == "__main__":
    main()