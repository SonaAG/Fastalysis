# --- Import all main features from each module ---
from app import (
	run_blast_top_hits,
	detect_sequence_type,
	load_user_sequence
)
from mutation import calculate_stats
from mu_table import (
	fetch_dbsnp_variants_for_gene,
	get_gene_symbol_from_refseq
)
from pub_med import (
	get_gene_symbol_from_accession,
	search_pubmed
)

# --- Controller functions for each feature ---
def run_blast_controller(fasta_path, num_hits=5):
	seq, seq_id = load_user_sequence(fasta_path)
	seq_type = detect_sequence_type(seq)
	hits = run_blast_top_hits(seq, seq_type, hitlist_size=num_hits)
	return {
		"sequence_id": seq_id,
		"sequence_type": seq_type,
		"hits": hits
	}

def run_mutation_controller(user_fasta_path, hit_index=0):

	from Bio import SeqIO, Entrez
	from Bio.Blast import NCBIWWW, NCBIXML
	import plotly.graph_objects as go
	import os
	Entrez.email = "your@email.com"
	user_seq, user_id = load_user_sequence(user_fasta_path)
	seq_type = detect_sequence_type(user_seq)

	# BLAST user sequence against RefSeq (as in mutation.py)
	program = "blastn" if seq_type == "nucleotide" else "blastp"
	db = "refseq_rna" if seq_type == "nucleotide" else "refseq_protein"

	result_handle = NCBIWWW.qblast(program, db, user_seq, hitlist_size=hit_index+1)
	blast_record = NCBIXML.read(result_handle)
	if not blast_record.alignments:
		return {"error": "No BLAST hits found for user sequence."}
	if hit_index >= len(blast_record.alignments):
		return {"error": f"Requested hit_index {hit_index} but only {len(blast_record.alignments)} hits found."}

	# Use the selected hit and its first HSP (matched region)
	selected_hit = blast_record.alignments[hit_index]
	hsp = selected_hit.hsps[0]
	ref_accession = selected_hit.accession
	ref_title = selected_hit.hit_def

	# Use only the aligned region (HSP) for analysis
	ref_aligned = hsp.sbjct
	user_aligned = hsp.query
	match_line = hsp.match
	ref_start = hsp.sbjct_start
	ref_end = hsp.sbjct_end
	score = hsp.score
	evalue = hsp.expect
	identities = hsp.identities
	align_len = hsp.align_length
	percent_identity = round((identities / align_len) * 100, 2)
	query_coverage = round((align_len / len(user_seq)) * 100, 2)

	# Calculate stats
	from mutation import calculate_stats
	length, identity, similarity, gaps = calculate_stats(ref_aligned, user_aligned)

	# Find mutations in the aligned region
	mutations = []
	for i, (r, u) in enumerate(zip(ref_aligned, user_aligned)):
		if r != u and r != '-' and u != '-':
			mutations.append({
				"position": i+1,
				"ref": r,
				"user": u,
				"mutation": f"c.{i+1}{r}>{u}"
			})

	# Plotly diagrams for the aligned region
	def plot_mutation_scatter(ref_seq, user_seq, ref_label="Ref", user_label="User"):
		positions = []
		hover_texts = []
		for i, (r, u) in enumerate(zip(ref_seq, user_seq)):
			if r != u and r != '-' and u != '-':
				positions.append(i+1)
				hover_texts.append(f"{ref_label}: {r} / {user_label}: {u}")
		fig = go.Figure()
		fig.add_trace(go.Scatter(
			x=positions,
			y=[1]*len(positions),
			mode='markers',
			marker=dict(color='purple', size=12, symbol='diamond'),
			name="Mutations",
			text=hover_texts,
			hovertemplate="Position: %{x}<br>%{text}<extra></extra>"
		))
		fig.update_layout(
			title="Mutation Scatter Plot (HSP region)",
			xaxis_title="Position",
			yaxis=dict(showticklabels=False),
			showlegend=False,
			height=400
		)
		return fig

	def plot_mutation_track(ref_seq, user_seq, ref_label="Ref", user_label="User"):
		positions = []
		mutation_types = []
		for i, (r, u) in enumerate(zip(ref_seq, user_seq)):
			if r != u and r != '-' and u != '-':
				positions.append(i+1)
				mutation_types.append('Substitution')
		fig = go.Figure()
		fig.add_trace(go.Bar(
			x=positions,
			y=[1]*len(positions),
			marker_color=['#f39c12']*len(positions),
			text=mutation_types,
			hovertemplate="Position: %{x}<br>Type: %{text}<extra></extra>",
			name="Mutations"
		))
		fig.update_layout(
			title="Mutation Track Diagram (HSP region)",
			xaxis_title="Position",
			yaxis=dict(showticklabels=False),
			showlegend=False,
			height=300
		)
		return fig

	# Save plots as HTML (for just the HSP region)
	scatter_html = f"mutation_scatter_{user_id}_{ref_accession}_hsp_{hit_index}.html"
	track_html = f"mutation_track_{user_id}_{ref_accession}_hsp_{hit_index}.html"
	plot_mutation_scatter(ref_aligned, user_aligned).write_html(scatter_html)
	plot_mutation_track(ref_aligned, user_aligned).write_html(track_html)

	return {
		"user_id": user_id,
		"ref_id": ref_accession,
		"ref_title": ref_title,
		"alignment": {
			"length": length,
			"identity": identity,
			"similarity": similarity,
			"gaps": gaps,
			"ref_aligned": ref_aligned,
			"user_aligned": user_aligned,
			"match_line": match_line,
			"ref_start": ref_start,
			"ref_end": ref_end,
			"score": score,
			"evalue": evalue,
			"percent_identity": percent_identity,
			"query_coverage": query_coverage
		},
		"mutations": mutations,
		"mutation_scatter_html": scatter_html,
		"mutation_track_html": track_html
	}

def run_variant_table_controller(accession_or_gene, by_accession=True, max_variants=20):
	if by_accession:
		gene_symbol = get_gene_symbol_from_refseq(accession_or_gene)
		if not gene_symbol:
			return {"error": "Could not extract gene symbol from accession."}
	else:
		gene_symbol = accession_or_gene
	variants = fetch_dbsnp_variants_for_gene(gene_symbol, max_records=max_variants)
	return {
		"gene_symbol": gene_symbol,
		"variants": variants
	}

def run_pubmed_controller(accession, additional_terms=None, max_results=10):
	gene_symbol = get_gene_symbol_from_accession(accession)
	if gene_symbol:
		search_term = gene_symbol
		if additional_terms:
			search_term += f" AND {additional_terms}"
	else:
		search_term = accession
		if additional_terms:
			search_term += f" AND {additional_terms}"
	papers = search_pubmed(search_term, max_results)
	return {
		"search_term": search_term,
		"papers": papers
	}

# --- Combined pipeline: run all steps and return results ---
def run_full_pipeline(fasta_path, ref_accession=None, variant_accession_or_gene=None, pubmed_accession=None, additional_pubmed_terms=None, num_blast_hits=5, max_variants=20, max_pubmed_results=10):
	"""
	Run all main Fastalysis features in sequence and return a dict with all results.
	Includes infection risk assessment and clinical recommendations.
	"""
	results = {}
	from risk_assessment import calculate_infection_risk
	# BLAST
	if fasta_path:
		results['blast'] = run_blast_controller(fasta_path, num_hits=num_blast_hits)
	# Mutation analysis (default to top hit)
	if fasta_path:
		results['mutation_top1'] = run_mutation_controller(fasta_path, hit_index=0)
		results['mutation_top2'] = run_mutation_controller(fasta_path, hit_index=1)
	# Variant table
	if variant_accession_or_gene:
		# If string looks like NM_ or NP_, treat as accession
		by_accession = variant_accession_or_gene.upper().startswith(('NM_', 'NP_', 'NC_', 'NR_'))
		results['variant_table'] = run_variant_table_controller(variant_accession_or_gene, by_accession=by_accession, max_variants=max_variants)
	# PubMed
	if pubmed_accession:
		results['pubmed'] = run_pubmed_controller(pubmed_accession, additional_terms=additional_pubmed_terms, max_results=max_pubmed_results)
	
	# Generate risk assessment if we have sequence analysis results
	if 'blast' in results and 'mutation_top1' in results:
		results['risk_assessment'] = calculate_infection_risk(
			blast_results=results['blast'],
			mutation_results=results['mutation_top1'],
			literature_data=results.get('pubmed')
		)
		
		# Add a human-readable summary
		risk_level = results['risk_assessment']['risk_level']
		risk_score = results['risk_assessment']['risk_score']
		confidence = results['risk_assessment']['confidence']
		
		summary = [
			"🧬 Genomic Analysis Summary",
			"------------------------",
			f"Risk Level: {risk_level} ({risk_score}%)",
			f"Confidence: {confidence}%",
			"",
			"Key Findings:",
		]
		
		for factor in results['risk_assessment']['risk_factors']:
			summary.append(f"• {factor}")
		
		summary.extend([
			"",
			"Recommendations:",
		])
		
		for rec in results['risk_assessment']['recommendations']:
			summary.append(f"• {rec}")
			
		results['summary'] = "\n".join(summary)
	
	return results

# --- Expose all main functions for external use ---
__all__ = [
	'run_blast_controller',
	'run_mutation_controller',
	'run_variant_table_controller',
	'run_pubmed_controller',
	'run_full_pipeline',
	'run_blast_top_hits',
	'detect_sequence_type',
	'load_user_sequence',
	'calculate_stats',
	'fetch_dbsnp_variants_for_gene',
	'get_gene_symbol_from_refseq',
	'get_gene_symbol_from_accession',
	'search_pubmed',
]
