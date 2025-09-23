import sys
import os
from Bio import SeqIO
from controller import run_full_pipeline
from tabulate import tabulate

def print_color(text, color):
    """Print colored text"""
    colors = {
        "red": "\033[91m",
        "green": "\033[92m",
        "yellow": "\033[93m",
        "reset": "\033[0m"
    }
    print(f"{colors.get(color, '')}{text}{colors['reset']}")

def test_sequence_analysis(sequence_file):
    """Run full analysis on a sequence file and display formatted results"""
    print("\n🧬 Running Genomic Analysis Pipeline...\n")
    
    try:
        # Run the full analysis pipeline
        results = run_full_pipeline(
            fasta_path=sequence_file,
            num_blast_hits=5,
            max_variants=10,
            max_pubmed_results=5
        )
        
        # Display summary if available
        if 'summary' in results:
            print("\n" + "=" * 50)
            print(results['summary'])
            print("=" * 50)
        
        # Display risk assessment details
        if 'risk_assessment' in results:
            risk = results['risk_assessment']
            print("\nDetailed Risk Assessment:")
            print("-" * 30)
            
            # Print risk score with color
            color = risk['color']
            score_text = f"Risk Score: {risk['risk_score']}% ({risk['risk_level']} Risk)"
            print_color(score_text, color)
            
            print(f"Confidence: {risk['confidence']}%")
            
            # Display analysis details in table format
            analysis = risk['analysis_details']
            if analysis['sequence_match']['pathogen_match']:
                match_data = [
                    ["Matching Pathogen", analysis['sequence_match']['matching_pathogen']],
                    ["Sequence Identity", f"{analysis['sequence_match']['highest_identity']:.1f}%"],
                    ["Mutations Found", analysis['mutations']['total_count']],
                    ["Critical Regions Affected", "Yes" if analysis['mutations']['critical_regions_affected'] else "No"]
                ]
                print("\nSequence Analysis:")
                print(tabulate(match_data, tablefmt="simple"))
            
            # Display BLAST hits
            if 'blast' in results and results['blast'].get('hits'):
                print("\nTop BLAST Hits:")
                blast_data = []
                for hit in results['blast']['hits'][:3]:  # Show top 3 hits
                    blast_data.append([
                        hit.get('title', 'Unknown'),
                        f"{hit.get('percent_identity', 0):.1f}%"
                    ])
                print(tabulate(blast_data, headers=['Organism', 'Identity'], tablefmt="simple"))
        
        # Display mutation details if available
        if 'mutation_top1' in results and results['mutation_top1'].get('mutations'):
            print("\nDetected Mutations:")
            mutations = results['mutation_top1']['mutations'][:5]  # Show top 5 mutations
            mutation_data = [[m['position'], f"{m['ref']} → {m['user']}", m.get('mutation', '')] for m in mutations]
            print(tabulate(mutation_data, headers=['Position', 'Change', 'Mutation'], tablefmt="simple"))
            
            if len(results['mutation_top1']['mutations']) > 5:
                print(f"... and {len(results['mutation_top1']['mutations']) - 5} more mutations")
        
    except Exception as e:
        print(f"❌ Error during analysis: {str(e)}")
        return

if __name__ == "__main__":
    # Test sequences
    test_sequences = [
        ("test_sequence.fasta", "ATGCGAGCTGACCCAGTACGTCGTACACGACTCGACGTCAGCTGACGTCAGCTAGCTAGACTGACT"),
        ("covid_spike.fasta", """ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAATCTTACAACCAGAACTCAAT
TACCCCCTGCATACACTAATTCTTTCACACGTGGTGTTTATTACCCTGACAAAGTTTTCAGATCC"""),
        ("random_sequence.fasta", "GGGGAAATTCCCCAAAATTTGGGGCCCCAAAATTTTGGGGAAAACCCCTTTTAAAAGGGCCCC")
    ]
    
    # Create and test each sequence
    for filename, sequence in test_sequences:
        print(f"\n📄 Testing sequence: {filename}")
        print("-" * 50)
        
        # Create temporary FASTA file
        with open(filename, "w") as f:
            f.write(f">test_sequence\n{sequence}\n")
        
        # Run analysis
        test_sequence_analysis(filename)
        
        # Clean up
        os.unlink(filename)
        
        input("\nPress Enter to continue to next sequence...")