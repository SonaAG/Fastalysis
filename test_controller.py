
from controller import (
    run_blast_controller,
    run_mutation_controller,
    run_variant_table_controller,
    run_pubmed_controller,
    run_full_pipeline
)
import json
import os

print("Current working directory:", os.getcwd())

try:
    # Example FASTA file path and accessions (update with your real data)
    fasta_path = "user_sequence.fasta"
    ref_accession = "NM_007294"
    gene_symbol = "BRCA1"

    results = {}

    # Test BLAST
    results['blast'] = run_blast_controller(fasta_path, num_hits=3)

    # Test Mutation Analysis
    results['mutation'] = run_mutation_controller(fasta_path, ref_accession)

    # Test Variant Table (by gene symbol)
    results['variant_table'] = run_variant_table_controller(gene_symbol, by_accession=False, max_variants=5)

    # Test PubMed Search
    results['pubmed'] = run_pubmed_controller(ref_accession, additional_terms="mutation", max_results=3)

    # Test Full Pipeline
    results['full_pipeline'] = run_full_pipeline(
        fasta_path,
        ref_accession=ref_accession,
        variant_accession_or_gene=gene_symbol,
        pubmed_accession=ref_accession,
        additional_pubmed_terms="mutation",
        num_blast_hits=3,
        max_variants=5,
        max_pubmed_results=3
    )

    # Write results to a file
    with open("controller_test_output.json", "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)

    print("Test results written to controller_test_output.json")
except Exception as e:
    print("An error occurred:", str(e))