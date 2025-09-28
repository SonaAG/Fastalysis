"""
Script to rebuild the RAG system database from scratch
"""

import os
import shutil
import traceback
from datetime import datetime

print("Starting RAG database rebuild...")
print(f"Current directory: {os.getcwd()}")

# Backup existing database
try:
    db_path = "./genomics_db"
    if os.path.exists(db_path):
        backup_path = f"./genomics_db_backup_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        print(f"Backing up existing database to {backup_path}")
        shutil.copytree(db_path, backup_path)
        print("Backup completed")
        
        # Delete existing database
        print("Removing existing database")
        shutil.rmtree(db_path)
        print("Existing database removed")
    
    # Create fresh database directory
    os.makedirs(db_path, exist_ok=True)
    print("Created fresh database directory")
    
    # Import and initialize RAG system
    from genomics_rag import GenomicsRAG
    print("Successfully imported GenomicsRAG")
    
    # Initialize RAG system with new database
    print("Initializing RAG system with new database...")
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
    
    print("\nRAG database rebuild completed successfully!")
    
except Exception as e:
    print(f"Error during database rebuild: {str(e)}")
    print(traceback.format_exc())