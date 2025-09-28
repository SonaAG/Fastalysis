"""
RAG (Retrieval-Augmented Generation) System for Genomics Knowledge
Provides context-aware knowledge retrieval for genes, diseases, mutations
"""

import chromadb
from chromadb.config import Settings
import json
import re
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
                "function": "Cystic Fibrosis Transmembrane Conductance Regulator, functions as a chloride channel and regulator of other channels",
                "diseases": ["Cystic fibrosis", "Congenital bilateral absence of vas deferens", "Pancreatitis"],
                "pathway": "Ion transport, epithelial fluid regulation",
                "clinical_significance": "Mutations cause cystic fibrosis, an autosomal recessive disease affecting the lungs, pancreas, and other organs. The most common mutation is F508del (deletion of phenylalanine at position 508)."
            },
            {
                "gene": "Dengue virus polyprotein",
                "function": "Encodes all viral proteins needed for replication",
                "diseases": ["Dengue fever", "Dengue hemorrhagic fever", "Dengue shock syndrome"],
                "pathway": "Viral replication and host immune evasion",
                "clinical_significance": "Single polyprotein is cleaved into structural proteins (C, prM, E) and non-structural proteins (NS1, NS2A, NS2B, NS3, NS4A, NS4B, NS5). Key target for antiviral drug development and vaccine design. Variations contribute to different serotypes (DENV-1, DENV-2, DENV-3, DENV-4)."
            }
        ]
        
        # Add to vector database
        for i, gene_info in enumerate(gene_knowledge):
            # Convert list fields to strings for ChromaDB compatibility
            metadata_safe = {k: (', '.join(v) if isinstance(v, list) else v) 
                            for k, v in gene_info.items()}
            
            text = f"Gene {gene_info['gene']}: {gene_info['function']}. Associated with {', '.join(gene_info['diseases'])}. Pathway: {gene_info['pathway']}. {gene_info['clinical_significance']}"
            
            self.genes_collection.add(
                documents=[text],
                metadatas=[metadata_safe],
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
            },
            {
                "type": "Silent mutation",
                "description": "Single nucleotide change that doesn't alter the amino acid sequence",
                "effect": "Usually no impact on protein function, though may affect mRNA splicing in some cases",
                "examples": "TP53 c.215C>G (p.Pro72Pro)"
            },
            {
                "type": "Splice site mutation",
                "description": "Affects the splice donor/acceptor sites, altering mRNA processing",
                "effect": "May lead to exon skipping, intron retention, or use of cryptic splice sites",
                "examples": "CFTR c.2657+5G>A, affecting exon 16 splicing"
            }
        ]
        
        for i, mut_info in enumerate(mutation_types):
            # Convert any list fields to strings for ChromaDB compatibility
            metadata_safe = {k: (', '.join(v) if isinstance(v, list) else v) 
                            for k, v in mut_info.items()}
                            
            text = f"{mut_info['type']}: {mut_info['description']}. Effect: {mut_info['effect']}. Example: {mut_info['examples']}"
            
            self.mutations_collection.add(
                documents=[text],
                metadatas=[metadata_safe], 
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
                
    def query_knowledge(self, query: str, collection_name: str = "genes", n_results: int = 3) -> Dict:
        """Query the knowledge base"""
        
        # Map collection name aliases to actual collection names
        collection_map = {
            "genes_info": "genes",
            "genes": "genes",
            "disease_associations": "diseases",
            "diseases": "diseases",
            "mutation_knowledge": "mutations",
            "mutations": "mutations",
            "literature": "literature"
        }
        
        try:
            # Get the correct collection name
            actual_collection_name = collection_map.get(collection_name, collection_name)
            
            # Correct known collection mapping issues
            if actual_collection_name == "genes_info":
                actual_collection_name = "genes"
                
            # Collection attribute name with "_collection" suffix
            collection_attr = f"{actual_collection_name}_collection"
            
            if not hasattr(self, collection_attr):
                print(f"Collection attribute '{collection_attr}' not found, trying default")
                collection_attr = "genes_collection"  # Fallback to genes collection
                
            collection = getattr(self, collection_attr)
            
            # Verify collection exists and has items
            if collection.count() == 0:
                print(f"Warning: Collection '{collection_attr}' is empty")
                return {"results": []}
                
            results = collection.query(
                query_texts=[query],
                n_results=min(n_results, collection.count())  # Prevent querying more than available
            )
            
            # Validate results structure before returning
            if results and 'documents' in results and len(results['documents']) > 0:
                return {
                    "query": query,
                    "results": results['documents'][0],
                    "metadata": results['metadatas'][0] if 'metadatas' in results and results['metadatas'] else [],
                    "distances": results['distances'][0] if 'distances' in results and results['distances'] else []
                }
            else:
                return {"query": query, "results": [], "metadata": [], "distances": []}
                
        except Exception as e:
            import traceback
            print(f"Error in query_knowledge for collection '{collection_name}': {str(e)}")
            print(traceback.format_exc())
            return {"query": query, "results": [], "metadata": [], "distances": []}
        
    def get_gene_context(self, gene: str) -> str:
        """Get comprehensive context about a gene"""
        
        try:
            contexts = []
            
            if not gene or not isinstance(gene, str) or len(gene.strip()) < 2:
                return "Insufficient information to retrieve context"
            
            # Special case handling for common genes and pathogens
            gene_lower = gene.lower()
            
            # Handle dengue virus specifically
            if "dengue" in gene_lower or "denv" in gene_lower:
                return """Dengue virus (DENV) is a single-stranded RNA virus of the Flaviviridae family, transmitted primarily by Aedes aegypti mosquitoes. 
                
There are four distinct serotypes (DENV-1, DENV-2, DENV-3, and DENV-4), each capable of causing the full spectrum of disease. 

The viral genome (~11kb) encodes a single polyprotein that is cleaved into:
- Three structural proteins: capsid (C), membrane precursor (prM), and envelope (E)
- Seven non-structural proteins (NS1, NS2A, NS2B, NS3, NS4A, NS4B, NS5) involved in viral replication

Primary infection typically causes self-limiting dengue fever, while secondary infection with a different serotype increases risk of severe dengue through antibody-dependent enhancement (ADE). DENV NS1 protein plays a critical role in viral replication and vascular leakage associated with severe dengue.

Research focuses on antiviral development targeting NS3 protease and NS5 polymerase, and vaccines must provide protection against all four serotypes."""
            
            # Handle CFTR gene specifically
            if gene_lower == "cftr" or "cystic fibrosis" in gene_lower:
                return """Gene CFTR (Cystic Fibrosis Transmembrane Conductance Regulator): Functions as a chloride channel and regulator of other channels in cell membranes. 

Associated with Cystic fibrosis, Congenital bilateral absence of vas deferens, and Pancreatitis.

Pathway: Ion transport, epithelial fluid regulation

Clinical significance: Mutations cause cystic fibrosis, an autosomal recessive disease affecting the lungs, pancreas, and other organs. The most common mutation is F508del (deletion of phenylalanine at position 508), which affects protein folding and trafficking.

Located on chromosome 7, the CFTR protein has 1,480 amino acids and contains two membrane-spanning domains, two nucleotide-binding domains, and a regulatory domain. Over 2,000 mutations have been identified, categorized into six classes based on their effect on the protein."""
                
            # Clean up the gene name/description to improve search
            cleaned_gene = re.sub(r'[^a-zA-Z0-9\s]', ' ', gene).strip()
            search_term = cleaned_gene[:50]  # Limit search term length
            
            # Try multiple collection names to ensure we get results
            # First try genes collection
            gene_info = self.query_knowledge(f"gene {search_term} function", "genes")
            if gene_info and 'results' in gene_info and gene_info['results']:
                contexts.append(f"Gene Information: {gene_info['results'][0]}")
            
            # Then try diseases collection    
            disease_info = self.query_knowledge(f"{search_term} disease", "diseases")
            if disease_info and 'results' in disease_info and disease_info['results']:
                contexts.append(f"Disease Associations: {disease_info['results'][0]}")
            
            # Then try literature collection
            lit_info = self.query_knowledge(f"{search_term} research", "literature")
            if lit_info and 'results' in lit_info and lit_info['results']:
                contexts.append(f"Recent Research: {lit_info['results'][0]}")
            
            # If we still have no results, try with different collection or variations of the name
            if not contexts:
                # Try searching all collections with slight variations
                variations = [search_term, search_term.split()[0] if ' ' in search_term else search_term]
                collection_names = ["genes", "diseases", "mutations", "literature"]
                
                for variation in variations:
                    for coll_name in collection_names:
                        info = self.query_knowledge(variation, coll_name)
                        if info and 'results' in info and info['results']:
                            contexts.append(f"Information from {coll_name}: {info['results'][0]}")
                            break
                    
                    if contexts:  # If we found something, stop trying variations
                        break
            
            # Return combined results or default message
            if contexts:
                return "\n\n".join(contexts)
            else:
                # Return a generic message if we couldn't find anything
                if len(search_term) >= 3:  # If it looks like a gene or protein name
                    return f"No specific information found for {gene} in the knowledge base. This appears to be a genetic sequence or protein that may require further analysis."
                else:
                    return f"No specific information found for {gene}. Please provide more details for a comprehensive analysis."
            
        except Exception as e:
            import traceback
            print(f"Error in get_gene_context: {str(e)}")
            print(traceback.format_exc())
            return f"Error retrieving context information. The knowledge base may need to be refreshed."
        
    def get_mutation_context(self, mutation_desc: str) -> str:
        """Get context about mutation types and effects"""
        
        try:
            if not mutation_desc or not isinstance(mutation_desc, str) or len(mutation_desc.strip()) < 2:
                return "Insufficient information to retrieve mutation context"
                
            # Handle common mutation patterns
            mutation_lower = mutation_desc.lower()
            
            # Special handling for common mutation types
            if "deletion" in mutation_lower or "del" in mutation_lower:
                special_context = """
Deletion: Removal of one or more nucleotides from the DNA sequence.
Impact: Can cause frameshift if not divisible by 3, potentially resulting in a premature stop codon and truncated protein. In-frame deletions may result in a protein missing specific amino acids.
Clinical significance: Depends on location and size, from benign to pathogenic. Deletions in critical domains often disrupt protein function.
Examples: CFTR F508del (loss of phenylalanine at position 508), common in cystic fibrosis.
                """
                return special_context.strip()
                
            elif "insertion" in mutation_lower or "ins" in mutation_lower:
                special_context = """
Insertion: Addition of one or more nucleotides to the DNA sequence.
Impact: Can cause frameshift if not divisible by 3, potentially resulting in a premature stop codon and truncated protein. In-frame insertions add extra amino acids that may disrupt protein folding.
Clinical significance: Depends on location and size, from benign to pathogenic. Insertions in critical domains often disrupt protein function.
Examples: BRCA1 c.5382insC (insertion of C after position 5382), associated with breast and ovarian cancer.
                """
                return special_context.strip()
                
            # Clean up the mutation description
            cleaned_desc = re.sub(r'[^a-zA-Z0-9\s]', ' ', mutation_desc).strip()
            search_term = cleaned_desc[:50]  # Limit search term length
            
            # Try multiple collections to ensure we get results
            # First try mutations collection
            results = self.query_knowledge(search_term, "mutations")
            if results and 'results' in results and results['results']:
                return f"Mutation Context: {results['results'][0]}"
                
            # If no results, try genes collection (may have mutation information)
            gene_results = self.query_knowledge(search_term, "genes")
            if gene_results and 'results' in gene_results and gene_results['results']:
                return f"Related Gene Context: {gene_results['results'][0]}"
                
            # If still no results, try literature collection
            lit_results = self.query_knowledge(search_term, "literature")
            if lit_results and 'results' in lit_results and lit_results['results']:
                return f"Research Context: {lit_results['results'][0]}"
                
            # Try to extract potential gene names from the mutation description
            # Common formats: GENE_NAME[space/colon/dot]MUTATION (e.g., "BRCA1 C61G" or "TP53.R175H")
            potential_gene = re.match(r'^([A-Za-z0-9_-]+)[:\.\s]', cleaned_desc)
            if potential_gene:
                gene_name = potential_gene.group(1)
                gene_context = self.get_gene_context(gene_name)
                if gene_context and "No specific information" not in gene_context:
                    return f"Related gene context for mutation in {gene_name}: {gene_context}"
            
            return "No specific mutation information found in knowledge base. Consider providing more details about the mutation or associated gene for better results."
        except Exception as e:
            import traceback
            print(f"Error in get_mutation_context: {str(e)}")
            print(traceback.format_exc())
            return "Error retrieving mutation context. The knowledge base may need to be refreshed."
            
    def initialize_knowledge_base(self):
        """Initialize with core genomics knowledge"""
        print("Initializing genomics knowledge base...")
        
        try:
            # Check if already populated
            if self.genes_collection.count() == 0:
                self.populate_gene_knowledge()
                print(f"Added {self.genes_collection.count()} gene records")
                
            if self.mutations_collection.count() == 0:
                self.populate_mutation_knowledge() 
                print(f"Added {self.mutations_collection.count()} mutation records")
                
            if self.diseases_collection.count() == 0:
                self.populate_disease_knowledge()
                print(f"Added {self.diseases_collection.count()} disease records")
                
            if self.literature_collection.count() == 0:
                self.populate_literature_knowledge()
                print(f"Added {self.literature_collection.count()} literature records")
                
            print("Knowledge base initialization complete!")
        except Exception as e:
            import traceback
            print(f"Error in initialize_knowledge_base: {str(e)}")
            print(traceback.format_exc())
            print("Will continue with limited knowledge base functionality")

    def populate_disease_knowledge(self):
        """Populate with disease information"""
        
        disease_knowledge = [
            {
                "disease": "Cystic Fibrosis",
                "description": "Genetic disorder affecting the lungs, pancreas, and other organs",
                "genetic_cause": "Mutations in CFTR gene",
                "inheritance": "Autosomal recessive",
                "symptoms": "Thick mucus in lungs and digestive tract, recurrent infections, poor growth",
                "treatment": "Airway clearance, antibiotics, CFTR modulators"
            },
            {
                "disease": "Dengue Fever",
                "description": "Mosquito-borne viral infection causing flu-like illness that can develop into severe dengue",
                "genetic_cause": "Dengue virus (DENV) - four distinct serotypes",
                "inheritance": "Not inherited, vector-borne",
                "symptoms": "High fever, severe headache, pain behind the eyes, muscle and joint pain, rash, mild bleeding",
                "treatment": "Supportive care, fluid management, no specific antiviral treatment"
            },
            {
                "disease": "Severe Dengue",
                "description": "Life-threatening complication of dengue infection",
                "genetic_cause": "Dengue virus, often from secondary infection with different serotype",
                "inheritance": "Not inherited, immune-mediated pathology",
                "symptoms": "Severe abdominal pain, persistent vomiting, rapid breathing, bleeding gums, fatigue, restlessness, blood in vomit, plasma leakage",
                "treatment": "Close monitoring, intravenous fluids, blood transfusion if needed"
            },
            {
                "disease": "Breast Cancer",
                "description": "Cancer that forms in the cells of the breasts",
                "genetic_cause": "Multiple genes including BRCA1, BRCA2, TP53",
                "inheritance": "Often sporadic, can be hereditary",
                "symptoms": "Breast lump, change in size or shape, skin changes, nipple discharge",
                "treatment": "Surgery, radiation, chemotherapy, hormone therapy, targeted therapy"
            }
        ]
        
        # Add to vector database
        for i, disease_info in enumerate(disease_knowledge):
            # Convert list fields to strings for ChromaDB compatibility
            metadata_safe = {k: (', '.join(v) if isinstance(v, list) else v) 
                            for k, v in disease_info.items()}
            
            text = f"Disease: {disease_info['disease']}. {disease_info['description']}. Genetic cause: {disease_info['genetic_cause']}. Inheritance: {disease_info['inheritance']}. Symptoms include {disease_info['symptoms']}. Treatment options: {disease_info['treatment']}."
            
            self.diseases_collection.add(
                documents=[text],
                metadatas=[metadata_safe],
                ids=[f"disease_{i}"]
            )
    
    def populate_literature_knowledge(self):
        """Populate with literature summaries"""
        
        literature_entries = [
            {
                "title": "Cystic fibrosis: a clinical view",
                "authors": "Davies JC, Alton EW, Bush A",
                "journal": "Lancet",
                "year": "2007",
                "summary": "Comprehensive review of cystic fibrosis pathophysiology, clinical manifestations, and treatment approaches. Highlights the role of CFTR mutations and potential therapies targeting the underlying genetic defect."
            },
            {
                "title": "Dengue virus: molecular basis of cell entry and pathogenesis",
                "authors": "Uno N, Ross TM",
                "journal": "Virology Journal",
                "year": "2018",
                "summary": "Review of dengue virus structure, mechanisms of cellular entry, and pathogenesis. Discusses how the virus evades immune responses and the differences between primary and secondary infections that may lead to severe dengue."
            },
            {
                "title": "Current status of dengue therapeutics research and development",
                "authors": "Low JG, Ooi EE, Vasudevan SG",
                "journal": "Journal of Infectious Diseases",
                "year": "2017",
                "summary": "Overview of dengue therapeutic approaches including antivirals, host-targeted therapies, and immune modulators. Discusses challenges in drug development due to the four different serotypes and complex pathophysiology."
            }
        ]
        
        # Add to vector database
        for i, lit_info in enumerate(literature_entries):
            # Convert list fields to strings for ChromaDB compatibility
            metadata_safe = {k: (', '.join(v) if isinstance(v, list) else v) 
                            for k, v in lit_info.items()}
            
            text = f"{lit_info['title']} by {lit_info['authors']} ({lit_info['journal']}, {lit_info['year']}). {lit_info['summary']}"
            
            self.literature_collection.add(
                documents=[text],
                metadatas=[metadata_safe],
                ids=[f"literature_{i}"]
            )

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

    def reload_knowledge_base(self):
        """Force reload of the knowledge base by clearing collections and repopulating them"""
        try:
            print("Reloading genomics knowledge base...")
            
            # Clear existing collections
            try:
                self.genes_collection.delete(ids=self.genes_collection.get()["ids"])
                self.diseases_collection.delete(ids=self.diseases_collection.get()["ids"])
                self.mutations_collection.delete(ids=self.mutations_collection.get()["ids"])
                self.literature_collection.delete(ids=self.literature_collection.get()["ids"])
            except Exception as e:
                print(f"Error clearing collections: {str(e)}")
                # In case of failure, recreate collections
                self.setup_collections()
            
            # Repopulate all collections
            self.populate_gene_knowledge()
            self.populate_mutation_knowledge()
            self.populate_disease_knowledge()
            self.populate_literature_knowledge()
            
            print("Knowledge base reload complete!")
            return True
        except Exception as e:
            import traceback
            print(f"Error in reload_knowledge_base: {str(e)}")
            print(traceback.format_exc())
            return False

# Example usage
if __name__ == "__main__":
    print("Testing GenomicsRAG system...")
    
    # Test the RAG system
    rag = GenomicsRAG()
    print("RAG system initialized")
    
    # Print collection status
    collections = {
        "genes_collection": rag.genes_collection,
        "diseases_collection": rag.diseases_collection,
        "mutations_collection": rag.mutations_collection,
        "literature_collection": rag.literature_collection
    }
    
    print("\nCollection status before initialization:")
    for name, collection in collections.items():
        try:
            count = collection.count()
            print(f"- {name}: {count} items")
        except Exception as e:
            print(f"- {name}: Error accessing ({str(e)})")
    
    # Initialize knowledge base
    print("\nInitializing knowledge base...")
    rag.initialize_knowledge_base()
    
    # Print collection status after initialization
    print("\nCollection status after initialization:")
    for name, collection in collections.items():
        try:
            count = collection.count()
            print(f"- {name}: {count} items")
        except Exception as e:
            print(f"- {name}: Error accessing ({str(e)})")
    
    # Test queries
    print("\nTesting queries:")
    
    print("\n1. Gene query (CFTR):")
    cftr_info = rag.get_gene_context("CFTR")
    print(cftr_info)
    
    print("\n2. Dengue virus query:")
    dengue_info = rag.get_gene_context("dengue virus")
    print(dengue_info)
    
    print("\n3. Mutation query (missense mutation):")
    mutation_info = rag.get_mutation_context("missense mutation")
    print(mutation_info)
    
    print("\n4. Mutation query (deletion):")
    deletion_info = rag.get_mutation_context("deletion")
    print(deletion_info)
    
    print("\nTest complete!")