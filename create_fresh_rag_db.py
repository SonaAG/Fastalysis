"""
Script to create a fresh RAG system database in a new location
"""

import os
import traceback
from datetime import datetime

print("Starting fresh RAG database creation...")
print(f"Current directory: {os.getcwd()}")

try:
    # Create a new database directory with timestamp
    db_path = f"./genomics_db_new_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    os.makedirs(db_path, exist_ok=True)
    print(f"Created fresh database directory: {db_path}")
    
    # Import and initialize RAG system
    from genomics_rag import GenomicsRAG
    print("Successfully imported GenomicsRAG")
    
    # Initialize RAG system with new database
    print(f"Initializing RAG system with new database at {db_path}...")
    rag = GenomicsRAG(db_path=db_path)
    print("RAG system initialized")
    
    # Initialize knowledge base
    print("Populating knowledge base...")
    rag.initialize_knowledge_base()
    print("Knowledge base populated")
    
    # Test the system
    print("\nTesting gene context retrieval...")
    cftr_info = rag.get_gene_context("CFTR")
    print(f"CFTR info successfully retrieved ({len(cftr_info)} chars)")
    
    print("\nTesting mutation context retrieval...")
    mutation_info = rag.get_mutation_context("missense mutation")
    print(f"Mutation info successfully retrieved ({len(mutation_info)} chars)")
    
    # Check final collection status
    collections = {
        "genes_collection": rag.genes_collection,
        "diseases_collection": rag.diseases_collection,
        "mutations_collection": rag.mutations_collection,
        "literature_collection": rag.literature_collection
    }
    
    print("\nFinal collection status:")
    for name, collection in collections.items():
        try:
            count = collection.count()
            print(f"- {name}: {count} items")
        except Exception as e:
            print(f"- {name}: Error accessing ({str(e)})")
    
    print("\nNew RAG database created successfully!")
    print(f"\nIMPORTANT: To use this database, update the db_path in enhanced_streamlit_app.py:")
    print(f"rag_system = GenomicsRAG(db_path=\"{db_path}\")")
    
except Exception as e:
    print(f"Error during database creation: {str(e)}")
    print(traceback.format_exc())