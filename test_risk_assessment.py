"""
Test the Risk Assessment Framework on user_sequence.fasta
"""

import sys
import os
from Bio import SeqIO
from controller import run_full_pipeline, run_blast_controller, run_mutation_controller
from tabulate import tabulate
from risk_assessment import InfectionRiskAssessment, calculate_infection_risk
from app import load_user_sequence, detect_sequence_type, run_blast_top_hits

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

def test_diverse_sequences():
    """
    Test the framework with diverse sequence types:
    - Bacterial gene
    - Viral protein
    - Fungal toxin
    - Human sequence (negative control)
    """
    print("Testing Risk Assessment Framework with Diverse Sequences")
    print("-------------------------------------------------------")
    
    # Create an instance of the risk assessor
    assessor = InfectionRiskAssessment()
    
    # 1. E. coli toxin (bacterial)
    print("\n1. Testing bacterial sequence (E. coli toxin):")
    ecoli_blast = {"hits": [{"title": "Shiga toxin [Escherichia coli O157:H7]", "percent_identity": 98.5}]}
    ecoli_mutations = {"mutations": [{"position": 45, "ref": "A", "alt": "T"}]}
    ecoli_result = assessor.calculate_infection_risk(ecoli_blast, ecoli_mutations)
    print(f"Risk Level: {ecoli_result['risk_level']}")
    print(f"Risk Score: {ecoli_result['risk_score']}")
    print(f"Confidence: {ecoli_result['confidence']}%")
    print("Risk Factors:")
    for factor in ecoli_result['risk_factors']:
        print(f"- {factor}")
    print("Recommendations:")
    for rec in ecoli_result['recommendations']:
        print(f"- {rec}")
    
    # 2. Influenza hemagglutinin (viral)
    print("\n2. Testing viral sequence (Influenza hemagglutinin):")
    flu_blast = {"hits": [{"title": "Hemagglutinin [Influenza A virus (A/California/07/2009(H1N1))]", "percent_identity": 96.0}]}
    flu_mutations = {"mutations": [{"position": 222, "ref": "D", "alt": "G"}]}
    flu_result = assessor.calculate_infection_risk(flu_blast, flu_mutations)
    print(f"Risk Level: {flu_result['risk_level']}")
    print(f"Risk Score: {flu_result['risk_score']}")
    print(f"Confidence: {flu_result['confidence']}%")
    print("Risk Factors:")
    for factor in flu_result['risk_factors']:
        print(f"- {factor}")
    print("Recommendations:")
    for rec in flu_result['recommendations']:
        print(f"- {rec}")
    
    # 3. Aspergillus fumigatus toxin (fungal)
    print("\n3. Testing fungal sequence (Aspergillus toxin):")
    fungal_blast = {"hits": [{"title": "Gliotoxin synthetase [Aspergillus fumigatus]", "percent_identity": 92.0}]}
    fungal_mutations = {"mutations": []}
    fungal_result = assessor.calculate_infection_risk(fungal_blast, fungal_mutations)
    print(f"Risk Level: {fungal_result['risk_level']}")
    print(f"Risk Score: {fungal_result['risk_score']}")
    print(f"Confidence: {fungal_result['confidence']}%")
    print("Risk Factors:")
    for factor in fungal_result['risk_factors']:
        print(f"- {factor}")
    print("Recommendations:")
    for rec in fungal_result['recommendations']:
        print(f"- {rec}")
    
    # 4. Human actin (negative control)
    print("\n4. Testing human sequence (Beta-actin, negative control):")
    human_blast = {"hits": [{"title": "Beta-actin [Homo sapiens]", "percent_identity": 100.0}]}
    human_mutations = {"mutations": [{"position": 10, "ref": "G", "alt": "A"}]}
    human_result = assessor.calculate_infection_risk(human_blast, human_mutations)
    print(f"Risk Level: {human_result['risk_level']}")
    print(f"Risk Score: {human_result['risk_score']}")
    print(f"Confidence: {human_result['confidence']}%")
    print("Risk Factors:")
    for factor in human_result['risk_factors']:
        print(f"- {factor}")
    print("Recommendations:")
    for rec in human_result['recommendations']:
        print(f"- {rec}")
    
    # 5. Test with real-world SARS-CoV-2 data
    print("\n5. Testing with SARS-CoV-2 spike protein:")
    covid_blast = {
        "hits": [
            {
                "title": "Spike glycoprotein [Severe acute respiratory syndrome coronavirus 2]",
                "percent_identity": 99.8,
                "evalue": 0.0
            }
        ]
    }
    covid_mutations = {
        "mutations": [
            {"position": 484, "ref": "E", "alt": "K"},
            {"position": 501, "ref": "N", "alt": "Y"},
            {"position": 614, "ref": "D", "alt": "G"}
        ]
    }
    covid_literature = {
        "papers": [
            {"title": "SARS-CoV-2 variants of concern", "abstract": "Analysis of mutations in the spike protein and their impact on transmissibility", "year": 2023},
            {"title": "COVID-19 pandemic and evolution of variants", "abstract": "Study of viral variants and their pathogenic properties", "year": 2022}
        ]
    }
    covid_result = assessor.calculate_infection_risk(covid_blast, covid_mutations, covid_literature)
    print(f"Risk Level: {covid_result['risk_level']}")
    print(f"Risk Score: {covid_result['risk_score']}")
    print(f"Confidence: {covid_result['confidence']}%")
    print("Risk Factors:")
    for factor in covid_result['risk_factors']:
        print(f"- {factor}")
    print("Recommendations:")
    for rec in covid_result['recommendations']:
        print(f"- {rec}")
    
    print("\nRisk Assessment Testing Complete!")

def run_risk_assessment_on_user_sequence():
    """
    Run risk assessment on the user_sequence.fasta file
    """
    print("\n" + "="*80)
    print("FASTALYSIS RISK ASSESSMENT ON USER SEQUENCE".center(80))
    print("="*80 + "\n")
    
    # Check if the file exists
    fasta_path = "user_sequence.fasta"
    if not os.path.exists(fasta_path):
        print(f"Error: File {fasta_path} not found!")
        print(f"Current directory: {os.getcwd()}")
        return
    
    print(f"Running risk assessment on: {fasta_path}")
    
    # Load the user sequence
    try:
        seq, seq_id = load_user_sequence(fasta_path)
        print(f"Loaded sequence: {seq_id}")
        print(f"Sequence type: {detect_sequence_type(seq)}")
        print(f"Sequence length: {len(seq)} bases")
        print(f"Sequence preview: {seq[:50]}...")
    except Exception as e:
        print(f"Error loading sequence: {str(e)}")
        return
        
    # Step 1: Run BLAST analysis
    print("\n[1/3] Running BLAST analysis...")
    try:
        blast_results = run_blast_controller(fasta_path, num_hits=5)
        print(f"BLAST complete. Found {len(blast_results.get('hits', []))} hits")
        
        # Display top hit
        if blast_results.get('hits'):
            top_hit = blast_results['hits'][0]
            print(f"Top hit: {top_hit.get('title', 'Unknown')}")
            print(f"Identity: {top_hit.get('percent_identity', 0):.1f}%")
            print(f"E-value: {top_hit.get('evalue', 'N/A')}")
    except Exception as e:
        print(f"Error running BLAST: {str(e)}")
        # Create minimal BLAST results to continue
        blast_results = {"hits": []}
    
    # Step 2: Run mutation analysis
    print("\n[2/3] Running mutation analysis...")
    try:
        mutation_results = run_mutation_controller(fasta_path)
        if "error" in mutation_results:
            print(f"Mutation analysis warning: {mutation_results['error']}")
            # Create empty mutation results
            mutation_results = {"mutations": []}
        else:
            # Make sure we have a 'mutations' key with a list
            if "mutations" not in mutation_results:
                mutation_results["mutations"] = []
            
            mutations = mutation_results.get("mutations", [])
            print(f"Mutation analysis complete. Found {len(mutations)} mutations")
            
            # Display a few mutations if available
            if mutations:
                print("Sample mutations:")
                for i, mut in enumerate(mutations[:3]):
                    print(f"  - {mut.get('ref', '?')}{mut.get('position', '?')}{mut.get('alt', '?')}")
                if len(mutations) > 3:
                    print(f"  - ... and {len(mutations) - 3} more")
    except Exception as e:
        print(f"Error running mutation analysis: {str(e)}")
        # Create empty mutation results
        mutation_results = {"mutations": []}
    
    # Step 3: Run risk assessment
    print("\n[3/3] Running risk assessment...")
    try:
        # Initialize the enhanced risk assessment class
        try:
            risk_assessor = InfectionRiskAssessment()
            result = risk_assessor.calculate_infection_risk(blast_results, mutation_results)
        except (NameError, AttributeError):
            # Fall back to the simple function if the class isn't available
            print("Using basic risk assessment function")
            result = calculate_infection_risk(blast_results, mutation_results)
            
        # Print results
        print("\n" + "-"*80)
        print("RISK ASSESSMENT RESULTS".center(80))
        print("-"*80)
        
        # Display risk level with color indicator
        risk_level = result.get("risk_level", "Unknown")
        risk_score = result.get("risk_score", 0)
        color = result.get("color", "gray")
        
        print_color(f"\nRisk Level: {risk_level.upper()} ({risk_score:.1f}/100)", color)
        print(f"Confidence: {result.get('confidence', 0):.1f}%")
        
        # Display risk color indicator
        color_indicator = {
            "red": "🔴 HIGH RISK",
            "yellow": "🟡 MODERATE RISK",
            "green": "🟢 LOW RISK",
            "gray": "⚪ UNKNOWN RISK"
        }
        print_color(f"\n{color_indicator.get(color, '⚪ UNKNOWN RISK')}", color)
        
        # Display risk factors
        print("\nRisk Factors:")
        if result.get("risk_factors"):
            for factor in result["risk_factors"]:
                print(f"  - {factor}")
        else:
            print("  - No specific risk factors identified")
        
        # Display recommendations
        print("\nRecommendations:")
        if result.get("recommendations"):
            for rec in result["recommendations"]:
                print(f"  - {rec}")
        else:
            print("  - No specific recommendations available")
            
        # Display analysis details
        print("\nAnalysis Details:")
        details = result.get("analysis_details", {})
        
        # Sequence match details
        seq_match = details.get("sequence_match", {})
        if seq_match:
            print(f"  - Highest sequence identity: {seq_match.get('highest_identity', 0):.1f}%")
            
        # Taxonomy details
        taxonomy = details.get("taxonomy", {})
        if taxonomy:
            if taxonomy.get("genus") and taxonomy.get("species"):
                print(f"  - Organism: {taxonomy.get('genus', '')} {taxonomy.get('species', '')}")
            if taxonomy.get("family"):
                print(f"  - Family: {taxonomy.get('family', '')}")
            if taxonomy.get("risk_group"):
                print(f"  - Biosafety Risk Group: {taxonomy.get('risk_group', '')}")
                
        # Mutation details
        mutations = details.get("mutations", {})
        if mutations:
            print(f"  - Total mutations: {mutations.get('total_count', 0)}")
            if mutations.get("in_functional_sites", 0) > 0:
                print(f"  - Mutations in functional sites: {mutations.get('in_functional_sites', 0)}")
                
        # Literature evidence
        lit = details.get("literature", {})
        if lit:
            print(f"  - Literature evidence score: {lit.get('evidence_score', 0)}/10")
            
    except Exception as e:
        print(f"Error running risk assessment: {str(e)}")
        import traceback
        traceback.print_exc()
        
    print("\n" + "="*80)
    print("END OF RISK ASSESSMENT".center(80))
    print("="*80 + "\n")

# Run the test when this file is executed
if __name__ == "__main__":
    # Run the assessment on the user_sequence.fasta file
    run_risk_assessment_on_user_sequence()
    
    # Optional: Uncomment below to run the synthetic test sequences too
    """
    # Test sequences
    test_sequences = [
        ("test_sequence.fasta", "ATGCGAGCTGACCCAGTACGTCGTACACGACTCGACGTCAGCTGACGTCAGCTAGCTAGACTGACT"),
        ("covid_spike.fasta", "ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAATCTTACAACCAGAACTCAATTACCCCCTGCATACACTAATTCTTTCACACGTGGTGTTTATTACCCTGACAAAGTTTTCAGATCC"),
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
    """
    
    # Uncomment to run diverse test sequences
    # test_diverse_sequences()