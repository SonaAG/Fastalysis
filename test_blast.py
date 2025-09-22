# Simple test to verify BLAST functionality

from app import run_blast_top_hits, detect_sequence_type

# Test with a known sequence
test_sequence = "ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGACGCGTACAGGAAACACAGAAAAAAGCCCGCACCTGACAGTGCGGGCTTTTTTTTCGACCAAAGGTAACGAGGTAACAACCATGCGAGTGTTGAAGTTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCGTGTTGCCGATATTCTGGAAAGCAATGCCAGGCAGGGGCAGGTGGCCACCGTCCTCTCTGCCCCCGCCAAAATCACCAACCACCTGGTGGCGATGATTGAAAAAACCATTAGCGGCCAGGATGCTTTACCCAATATCAGCGATGCCGAACGTATTTTTGCCGAACTTTTGACGGGACTCGCCGCCGCCCAGCCGGGGTTCCCGCTGGCGCAATTGAAAACTTTCGTCGATCAGGAATTTGCCCAAATAAAACATGTCCTGCATGGCATTAGTTTGTTGGTTCTTCATCGGCCCGCAGGGGCCGTAAACACCTGCCTTCTTCCTTAAAACCATACCTACTGGTGCCGATCCTGTGCCTTTTCGTCATGCCATGTGACCCTGAGGGCTTAAACTTCCGAGATAAGGGGCGTTTGCTGAAAGGTCTGTATGGGAACCTTAAGGAAAAGGAAGGGGTAAAAAGGTGCATGAAG"

print(f"Testing BLAST with sequence length: {len(test_sequence)}")
print(f"Sequence type: {detect_sequence_type(test_sequence)}")

# Run BLAST
hits = run_blast_top_hits(test_sequence, "nucleotide", hitlist_size=3)

print(f"\nBLAST Results: Found {len(hits)} hits")

for i, hit in enumerate(hits):
    print(f"\nHit {i+1}:")
    print(f"  Accession: {hit.get('accession', 'N/A')}")
    print(f"  Identity: {hit.get('percent_identity', 0)}%")
    print(f"  E-value: {hit.get('evalue', 'N/A')}")
    print(f"  Title: {hit.get('title', 'N/A')[:100]}...")

if len(hits) == 0:
    print("❌ BLAST failed - no hits found")
else:
    print(f"✅ BLAST working - found {len(hits)} hits")