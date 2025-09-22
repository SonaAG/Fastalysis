"""
RAG (Retrieval-Augmented Generation) System for Genomics Knowledge
Provides context-aware knowledge retrieval for genes, diseases, mutations
"""

import chromadb
from chromadb.config import Settings
import json
from typing import List, Dict, Optional
from sentence_transformers import SentenceTransformer
import requests
from pathlib import Path

class GenomicsRAG:
    def __init__(self, db_path: str = "./genomics_db"):
        """Initialize RAG system with vector database"""
        self.db_path = Path(db_path)
        self.db_path.mkdir(exist_ok=True)
        
        # Initialize ChromaDB
        self.client = chromadb.PersistentClient(path=str(self.db_path))
        
        # Initialize embedding model
        self.embedding_model = SentenceTransformer('all-MiniLM-L6-v2')
        
        # Create collections
        self.setup_collections()
        
    def setup_collections(self):
        """Create vector database collections"""
        
        # Gene information collection
        self.genes_collection = self.client.get_or_create_collection(
            name="genes_info",
            metadata={"description": "Gene functions, pathways, diseases"}
        )
        
        # Disease associations
        self.diseases_collection = self.client.get_or_create_collection(
            name="disease_associations", 
            metadata={"description": "Gene-disease relationships"}
        )
        
        # Mutation knowledge
        self.mutations_collection = self.client.get_or_create_collection(
            name="mutation_knowledge",
            metadata={"description": "Mutation types, effects, significance"}
        )
        
        # Literature abstracts
        self.literature_collection = self.client.get_or_create_collection(
            name="literature",
            metadata={"description": "Research paper abstracts and findings"}
        )
        
    def populate_gene_knowledge(self):
        """Populate database with essential gene information"""
        
        # Core genes important in genomics research
        gene_knowledge = [
            {
                "gene": "BRCA1",
                "function": "DNA repair, tumor suppression",
                "diseases": ["Breast cancer", "Ovarian cancer"],
                "pathway": "Homologous recombination DNA repair",
                "clinical_significance": "High penetrance cancer susceptibility gene"
            },
            {
                "gene": "BRCA2", 
                "function": "DNA repair, homologous recombination",
                "diseases": ["Breast cancer", "Ovarian cancer", "Prostate cancer"],
                "pathway": "Homologous recombination DNA repair",
                "clinical_significance": "High penetrance cancer susceptibility gene"
            },
            {
                "gene": "TP53",
                "function": "Cell cycle regulation, apoptosis",
                "diseases": ["Li-Fraumeni syndrome", "Various cancers"],
                "pathway": "p53 signaling pathway",
                "clinical_significance": "Guardian of the genome, most mutated gene in cancer"
            },
            {
                "gene": "CFTR",
                "function": "Chloride channel",
                "diseases": ["Cystic fibrosis"],
                "pathway": "Ion transport",
                "clinical_significance": "Autosomal recessive disease gene"
            }
        ]
        
        # Add to vector database
        for i, gene_info in enumerate(gene_knowledge):
            text = f"Gene {gene_info['gene']}: {gene_info['function']}. Associated with {', '.join(gene_info['diseases'])}. Pathway: {gene_info['pathway']}. {gene_info['clinical_significance']}"
            
            self.genes_collection.add(
                documents=[text],
                metadatas=[gene_info],
                ids=[f"gene_{i}"]
            )
            
    def populate_mutation_knowledge(self):
        """Populate with mutation type information"""
        
        mutation_types = [
            {
                "type": "Missense mutation",
                "description": "Single nucleotide change resulting in amino acid substitution",
                "effect": "May alter protein function depending on location and amino acid properties",
                "examples": "BRCA1 p.Cys61Gly"
            },
            {
                "type": "Nonsense mutation", 
                "description": "Creates premature stop codon",
                "effect": "Usually results in truncated, non-functional protein",
                "examples": "BRCA1 p.Arg1751*"
            },
            {
                "type": "Frameshift mutation",
                "description": "Insertion or deletion not divisible by 3, shifts reading frame",
                "effect": "Usually results in truncated protein",
                "examples": "BRCA1 c.185delAG"
            }
        ]
        
        for i, mut_info in enumerate(mutation_types):
            text = f"{mut_info['type']}: {mut_info['description']}. Effect: {mut_info['effect']}. Example: {mut_info['examples']}"
            
            self.mutations_collection.add(
                documents=[text],
                metadatas=[mut_info], 
                ids=[f"mutation_{i}"]
            )
    
    def add_pubmed_abstracts(self, abstracts: List[Dict]):
        """Add PubMed abstracts to literature collection"""
        
        for i, abstract in enumerate(abstracts):
            if 'abstract' in abstract and 'title' in abstract:
                text = f"{abstract['title']} {abstract['abstract']}"
                
                self.literature_collection.add(
                    documents=[text],
                    metadatas=[abstract],
                    ids=[f"pubmed_{abstract.get('pmid', i)}"]
                )
                
    def query_knowledge(self, query: str, collection_name: str = "genes_info", n_results: int = 3) -> Dict:
        """Query the knowledge base"""
        
        collection = getattr(self, f"{collection_name}_collection")
        
        results = collection.query(
            query_texts=[query],
            n_results=n_results
        )
        
        return {
            "query": query,
            "results": results['documents'][0],
            "metadata": results['metadatas'][0],
            "distances": results['distances'][0] if results['distances'] else []
        }
        
    def get_gene_context(self, gene: str) -> str:
        """Get comprehensive context about a gene"""
        
        contexts = []
        
        # Query gene information
        gene_info = self.query_knowledge(f"gene {gene} function", "genes_info")
        if gene_info['results']:
            contexts.append(f"Gene Information: {gene_info['results'][0]}")
            
        # Query disease associations  
        disease_info = self.query_knowledge(f"{gene} disease", "disease_associations")
        if disease_info['results']:
            contexts.append(f"Disease Associations: {disease_info['results'][0]}")
            
        # Query recent literature
        lit_info = self.query_knowledge(f"{gene} research", "literature")
        if lit_info['results']:
            contexts.append(f"Recent Research: {lit_info['results'][0]}")
            
        return "\n\n".join(contexts)
        
    def get_mutation_context(self, mutation_desc: str) -> str:
        """Get context about mutation types and effects"""
        
        results = self.query_knowledge(mutation_desc, "mutation_knowledge")
        
        if results['results']:
            return f"Mutation Context: {results['results'][0]}"
        else:
            return "No specific mutation information found in knowledge base."
            
    def initialize_knowledge_base(self):
        """Initialize with core genomics knowledge"""
        print("Initializing genomics knowledge base...")
        
        # Check if already populated
        if self.genes_collection.count() == 0:
            self.populate_gene_knowledge()
            print(f"Added {self.genes_collection.count()} gene records")
            
        if self.mutations_collection.count() == 0:
            self.populate_mutation_knowledge() 
            print(f"Added {self.mutations_collection.count()} mutation records")
            
        print("Knowledge base initialization complete!")

# Integration with your existing system
class EnhancedGenomicsController:
    """Enhanced controller that combines your existing analysis with RAG"""
    
    def __init__(self):
        self.rag = GenomicsRAG()
        self.rag.initialize_knowledge_base()
        
    def enhanced_blast_analysis(self, sequence: str, num_hits: int = 5) -> Dict:
        """BLAST analysis with contextual knowledge"""
        from controller import run_blast_controller
        import tempfile
        import os
        
        # Run standard BLAST
        temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
        temp_file.write(f">user_sequence\n{sequence}")
        temp_file.close()
        
        blast_result = run_blast_controller(temp_file.name, num_hits=num_hits)
        os.unlink(temp_file.name)
        
        # Add contextual knowledge
        enhanced_result = blast_result.copy()
        
        if blast_result.get('hits'):
            for hit in blast_result['hits']:
                # Extract gene info from hit title if possible
                if 'accession' in hit:
                    gene_context = self.rag.get_gene_context(hit.get('title', ''))
                    hit['knowledge_context'] = gene_context
                    
        return enhanced_result
        
    def enhanced_mutation_analysis(self, sequence: str, hit_index: int = 0) -> Dict:
        """Mutation analysis with biological context"""
        from controller import run_mutation_controller
        import tempfile
        import os
        
        # Run standard mutation analysis
        temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
        temp_file.write(f">user_sequence\n{sequence}")
        temp_file.close()
        
        mutation_result = run_mutation_controller(temp_file.name, hit_index=hit_index)
        os.unlink(temp_file.name)
        
        # Add mutation context
        enhanced_result = mutation_result.copy()
        
        if mutation_result.get('mutations'):
            mutation_contexts = []
            for mutation in mutation_result['mutations']:
                context = self.rag.get_mutation_context(mutation.get('mutation', ''))
                mutation_contexts.append(context)
            enhanced_result['mutation_contexts'] = mutation_contexts
            
        return enhanced_result

# Example usage
if __name__ == "__main__":
    # Test the RAG system
    rag = GenomicsRAG()
    rag.initialize_knowledge_base()
    
    # Test queries
    gene_info = rag.query_knowledge("BRCA1 cancer function")
    print("Gene query result:", json.dumps(gene_info, indent=2))
    
    mutation_info = rag.query_knowledge("missense mutation effect")
    print("Mutation query result:", json.dumps(mutation_info, indent=2))