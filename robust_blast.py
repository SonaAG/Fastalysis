"""
Helper module for robust BLAST searches with fallback mechanisms
"""

import os
import tempfile
from Bio import SeqIO, Entrez
from Bio.Blast import NCBIWWW, NCBIXML
import time
import requests
import json
from typing import Dict, List, Optional, Union, Any
import socket
import urllib.error

# Set a longer timeout for BLAST requests
socket.setdefaulttimeout(120)  # 2 minute timeout

class RobustBLASTSearch:
    """BLAST search with multiple fallback mechanisms"""
    
    def __init__(self, email: str = "your.email@example.com"):
        """Initialize with NCBI requirements"""
        self.email = email
        Entrez.email = email
        # Track which method succeeded
        self.search_method_used = None
    
    def run_blast_search(self, 
                         sequence: str, 
                         max_retries: int = 3, 
                         hitlist_size: int = 5) -> Dict[str, Any]:
        """
        Perform BLAST search with multiple fallback mechanisms
        
        Args:
            sequence: The sequence to search
            max_retries: Number of times to retry NCBI BLAST
            hitlist_size: Number of hits to return
            
        Returns:
            Dictionary with search results
        """
        # First detect sequence type
        seq_type = self._detect_sequence_type(sequence)
        
        # Try direct NCBI BLAST first (most reliable)
        for attempt in range(max_retries):
            try:
                print(f"BLAST attempt {attempt+1} via direct NCBI API")
                results = self._run_ncbi_blast(sequence, seq_type, hitlist_size)
                self.search_method_used = "direct_ncbi"
                return {
                    "status": "success", 
                    "data": {
                        "hits": results,
                        "method": "ncbi_direct",
                        "seq_type": seq_type
                    }
                }
            except (urllib.error.URLError, socket.timeout, ConnectionError) as e:
                print(f"NCBI BLAST connection error (attempt {attempt+1}): {str(e)}")
                time.sleep(2)  # Wait before retry
        
        # If direct BLAST failed, try using local API if available
        try:
            print("Trying local FastAPI BLAST endpoint")
            results = self._run_local_api_blast(sequence, hitlist_size)
            if results:
                self.search_method_used = "local_api"
                return {
                    "status": "success", 
                    "data": {
                        "hits": results,
                        "method": "local_api",
                        "seq_type": seq_type
                    }
                }
        except Exception as e:
            print(f"Local API BLAST failed: {str(e)}")
        
        # If all else fails, use alternative BLAST provider
        try:
            print("Trying alternative BLAST provider")
            results = self._run_alternative_blast(sequence, seq_type, hitlist_size)
            if results:
                self.search_method_used = "alternative_provider"
                return {
                    "status": "success", 
                    "data": {
                        "hits": results,
                        "method": "alternative_provider",
                        "seq_type": seq_type
                    }
                }
        except Exception as e:
            print(f"Alternative BLAST provider failed: {str(e)}")
        
        # Last resort - generate placeholder result based on sequence
        self.search_method_used = "fallback_analysis"
        return {
            "status": "partial", 
            "data": {
                "hits": self._generate_fallback_results(sequence, seq_type),
                "method": "fallback_analysis",
                "seq_type": seq_type
            }
        }
    
    def _detect_sequence_type(self, seq: str) -> str:
        """Determine if sequence is nucleotide or protein"""
        protein_set = set("ACDEFGHIKLMNPQRSTVWY")
        nucleotide_set = set("ACGTU")
        seq_set = set(seq.upper())
        
        if seq_set.issubset(nucleotide_set) or len(seq_set & nucleotide_set) / len(seq_set) > 0.95:
            return "nucleotide"
        elif seq_set.issubset(protein_set):
            return "protein"
        else:
            # Default to nucleotide if ambiguous
            return "nucleotide"
    
    def _run_ncbi_blast(self, seq: str, seq_type: str, hitlist_size: int) -> List[Dict]:
        """Run BLAST search via NCBI API directly"""
        program = "blastn" if seq_type == "nucleotide" else "blastp"
        db = "nt" if seq_type == "nucleotide" else "nr"
        
        print(f"Running {program} against {db} database")
        result_handle = NCBIWWW.qblast(program, db, seq, hitlist_size=hitlist_size)
        blast_record = NCBIXML.read(result_handle)
        
        hits = []
        if blast_record.alignments:
            for i, alignment in enumerate(blast_record.alignments):
                if not alignment.hsps:
                    continue
                
                hsp = alignment.hsps[0]  # Best HSP
                
                # Calculate metrics
                identities = hsp.identities
                align_len = hsp.align_length
                percent_identity = round((identities / align_len) * 100, 2) if align_len > 0 else 0
                query_coverage = round((align_len / len(seq)) * 100, 2) if len(seq) > 0 else 0
                
                hit_data = {
                    "accession": alignment.accession,
                    "title": alignment.hit_def,
                    "percent_identity": percent_identity,
                    "e_value": hsp.expect,
                    "score": hsp.score,
                    "query_coverage": query_coverage,
                }
                hits.append(hit_data)
                
        return hits
    
    def _run_local_api_blast(self, seq: str, hitlist_size: int) -> List[Dict]:
        """Try running BLAST through local API if available"""
        try:
            response = requests.post(
                "http://localhost:8000/blast",
                json={"sequence": seq, "num_hits": hitlist_size},
                timeout=30
            )
            
            if response.status_code == 200:
                result = response.json()
                if "data" in result and "hits" in result["data"]:
                    return result["data"]["hits"]
            return []
        except requests.exceptions.RequestException:
            return []
    
    def _run_alternative_blast(self, seq: str, seq_type: str, hitlist_size: int) -> List[Dict]:
        """
        Try an alternative BLAST provider
        This is a placeholder - in a real implementation you might use
        EBI's BLAST service or another alternative
        """
        # Currently returns empty results as this is just a placeholder
        # In a real implementation, you would integrate with an alternative service
        return []
    
    def _generate_fallback_results(self, seq: str, seq_type: str) -> List[Dict]:
        """
        Generate fallback results when all BLAST methods fail
        Performs basic sequence analysis to provide some useful information
        """
        # Calculate basic sequence properties
        seq_len = len(seq)
        gc_content = (seq.upper().count('G') + seq.upper().count('C')) / seq_len * 100 if seq_type == "nucleotide" else 0
        
        # Create a single fallback "hit" with basic sequence info
        fallback_hit = {
            "accession": "OFFLINE_ANALYSIS",
            "title": "Network Error - Basic Sequence Analysis Only",
            "percent_identity": None,
            "e_value": None,
            "score": None,
            "query_coverage": 100,
            "sequence_length": seq_len,
            "gc_content": round(gc_content, 2) if seq_type == "nucleotide" else None,
            "is_fallback": True,
            "message": "Unable to connect to BLAST servers. This is an offline analysis only."
        }
        
        return [fallback_hit]