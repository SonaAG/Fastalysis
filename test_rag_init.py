"""
Simple test script to check RAG system initialization
"""

import sys
import os
import traceback

print("Testing RAG system initialization...")
print(f"Current directory: {os.getcwd()}")

try:
    from genomics_rag import GenomicsRAG
    print("Successfully imported GenomicsRAG")
    
    # Create DB directory if it doesn't exist
    os.makedirs("./genomics_db", exist_ok=True)
    print("Database directory created/verified")
    
    # Initialize RAG system
    rag = GenomicsRAG(db_path="./genomics_db")
    print("RAG system initialized")
    
    # Check collections
    collections = {
        "genes_collection": rag.genes_collection,
        "diseases_collection": rag.diseases_collection,
        "mutations_collection": rag.mutations_collection,
        "literature_collection": rag.literature_collection
    }
    
    print("\nCollection status:")
    for name, collection in collections.items():
        try:
            count = collection.count()
            print(f"- {name}: {count} items")
        except Exception as e:
            print(f"- {name}: Error accessing ({str(e)})")
    
    print("\nTest complete!")
except Exception as e:
    print(f"Error: {str(e)}")
    print(traceback.format_exc())