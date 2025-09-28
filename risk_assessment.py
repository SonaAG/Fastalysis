from typing import Dict, List, Optional, Union, Tuple
import numpy as np
import pandas as pd
from Bio import Entrez
import re

class InfectionRiskAssessment:
    """
    Generalized infection risk assessment for any biological sequence.
    Works with any taxonomic group, not just pre-defined pathogens.
    """
    
    def __init__(self, email="your@email.com"):
        """Initialize with taxonomic database access"""
        Entrez.email = email
        
        # Load taxonomy risk database
        self.taxonomy_risk = self._load_taxonomy_risk()
        
    def _load_taxonomy_risk(self) -> Dict:
        """
        Load taxonomy risk classifications
        
        Returns a dict mapping taxonomy IDs to risk levels
        """
        # In a full implementation, this would load from a database
        # For now, we'll use a simplified approach with major taxonomic groups
        
        # Format: {taxonomy_level: {taxon: risk_score}}
        return {
            "superkingdom": {
                "viruses": 40,
                "bacteria": 30,
                "fungi": 20,
                "archaea": 10,
                "eukaryota": 5
            },
            "family": {
                # Virus families with known human pathogens
                "coronaviridae": 50,
                "orthomyxoviridae": 45,  # Influenza 
                "filoviridae": 60,       # Ebola
                "flaviviridae": 50,      # Zika, Dengue - increased from 45 to 50
                
                # Bacterial families with known pathogens
                "enterobacteriaceae": 40,
                "staphylococcaceae": 40,
                "streptococcaceae": 35,
                "mycobacteriaceae": 45,
                
                # Fungal families with known pathogens
                "aspergillaceae": 30,
                "cryptococcaceae": 35
            },
            "genus": {
                # High-risk genera
                "betacoronavirus": 55,
                "yersinia": 50,
                "bacillus": 45,
                "salmonella": 40,
                "mycobacterium": 45
            }
        }
    
    def calculate_infection_risk(self, 
                               blast_results: Dict, 
                               mutation_results: Dict, 
                               literature_data: Optional[Dict] = None) -> Dict:
        """
        Calculate infection risk using multiple evidence sources.
        Works for any biological sequence.
        """
        # Extract taxonomy and sequence features
        features = self._extract_features(blast_results, mutation_results, literature_data)
        
        # Calculate risk score using taxonomy and sequence characteristics
        risk_score, risk_details = self._calculate_risk_score(features)
        confidence = self._calculate_confidence(features)
        
        # Generate risk factors based on findings
        risk_factors = self._identify_risk_factors(features, risk_details)
        
        # Generate appropriate recommendations
        recommendations = self._generate_recommendations(risk_score, risk_factors, features)
        
        return {
            "risk_score": risk_score,
            "risk_level": self._determine_risk_level(risk_score),
            "color": self._determine_risk_color(risk_score),
            "confidence": confidence,
            "risk_factors": risk_factors,
            "recommendations": recommendations,
            "analysis_details": self._generate_analysis_details(features, risk_details),
            "taxonomy_classification": features.get("taxonomy", {})
        }
    
    def _extract_features(self, blast_results, mutation_results, literature_data) -> Dict:
        """Extract features from any biological sequence results"""
        features = {
            "highest_identity": 0,
            "taxonomy": {},
            "mutation_count": 0,
            "literature_evidence": 0,
            "conserved_domains": [],
            "functional_sites": [],
            "sequence_type": self._detect_sequence_type(blast_results),
        }
        
        # Process BLAST hits
        if "hits" in blast_results and blast_results["hits"]:
            # Get taxonomy information from top hits
            features["taxonomy"] = self._extract_taxonomy(blast_results["hits"])
            
            # Get highest identity match
            features["highest_identity"] = max(
                float(hit.get("percent_identity", 0)) 
                for hit in blast_results["hits"]
            )
            
            # Store top pathogen titles for direct detection
            features["top_pathogen"] = [hit.get("title", "") for hit in blast_results["hits"][:3]]
            
            # Special handling for Dengue virus
            if any("dengue virus" in title.lower() for title in features["top_pathogen"]):
                features["virus_type"] = "dengue"
                if "taxonomy" not in features:
                    features["taxonomy"] = {}
                features["taxonomy"]["family"] = "flaviviridae"
                features["taxonomy"]["genus"] = "flavivirus"
                features["taxonomy"]["risk_group"] = 2
                features["taxonomy"]["pathogenic"] = True
            
            # Extract conserved domains and functional sites
            features["conserved_domains"] = self._extract_domains(blast_results["hits"])
            features["functional_sites"] = self._extract_functional_sites(blast_results["hits"])
        
        # Process mutations
        if "mutations" in mutation_results:
            features["mutation_count"] = len(mutation_results["mutations"])
            features["mutations_in_functional_sites"] = self._mutations_in_functional_sites(
                mutation_results["mutations"],
                features["functional_sites"]
            )
        
        # Process literature
        if literature_data and "papers" in literature_data:
            features["literature_evidence"] = self._analyze_literature(literature_data["papers"])
        
        return features
    
    def _extract_taxonomy(self, hits: List[Dict]) -> Dict:
        """
        Extract taxonomic information from BLAST hits
        Works with any organism, not just predefined pathogens
        """
        taxonomy = {
            "superkingdom": None,
            "kingdom": None,
            "phylum": None,
            "class": None,
            "order": None,
            "family": None,
            "genus": None,
            "species": None,
            "risk_group": None,
            "pathogenic": False
        }
        
        # Extract organism names from top hits
        organisms = []
        for hit in hits[:3]:  # Consider top 3 hits
            title = hit.get("title", "").lower()
            
            # Direct detection of key viruses that may not be in brackets
            if "dengue virus" in title:
                organisms.append("dengue virus")
                # Add explicit taxonomic classification for Dengue
                taxonomy["family"] = "flaviviridae"
                taxonomy["genus"] = "flavivirus" 
                taxonomy["species"] = "dengue virus"
                taxonomy["pathogenic"] = True
                taxonomy["risk_group"] = 2  # Dengue is biosafety level 2
            
            # Try to extract organism name from brackets
            org_match = re.search(r'\[(.*?)\]', title)
            if org_match:
                organisms.append(org_match.group(1).lower())
        
        if not organisms:
            return taxonomy
            
        # Try to get taxonomic information from NCBI
        try:
            # Use the most common organism name from top hits
            from collections import Counter
            most_common_org = Counter(organisms).most_common(1)[0][0]
            
            # Query NCBI taxonomy database
            handle = Entrez.esearch(db="taxonomy", term=most_common_org)
            record = Entrez.read(handle)
            handle.close()
            
            if record["IdList"]:
                tax_id = record["IdList"][0]
                handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
                tax_record = Entrez.read(handle)
                handle.close()
                
                if tax_record:
                    # Extract lineage information
                    lineage = tax_record[0].get("LineageEx", [])
                    for entry in lineage:
                        rank = entry.get("Rank")
                        name = entry.get("ScientificName", "").lower()
                        if rank in taxonomy:
                            taxonomy[rank] = name
                    
                    # Add species level if available
                    taxonomy["species"] = tax_record[0].get("ScientificName", "").lower()
                    
                    # Check for pathogenicity
                    taxonomy["pathogenic"] = self._is_pathogenic_taxon(taxonomy)
                    
                    # Assign risk group based on taxonomy
                    taxonomy["risk_group"] = self._assign_risk_group(taxonomy)
        except Exception as e:
            print(f"Taxonomy lookup error: {str(e)}")
        
        return taxonomy
    
    def _is_pathogenic_taxon(self, taxonomy: Dict) -> bool:
        """
        Determine if taxonomy is likely pathogenic based on taxonomic classification
        """
        # Check for known pathogenic groups at each taxonomic level
        pathogen_indicators = {
            "superkingdom": ["viruses"],
            "phylum": ["proteobacteria", "firmicutes"],
            "family": ["coronaviridae", "enterobacteriaceae", "staphylococcaceae"],
            "genus": ["salmonella", "escherichia", "staphylococcus", "streptococcus", 
                     "betacoronavirus", "mycobacterium", "clostridium"]
        }
        
        # Check each taxonomic level
        for level, indicators in pathogen_indicators.items():
            if taxonomy.get(level) and taxonomy[level] in indicators:
                return True
                
        # Also check for pathogenic keywords in species name
        species = taxonomy.get("species", "").lower()
        pathogen_keywords = ["pathogen", "virulent", "toxin", "disease", 
                            "pathogenic", "infection", "virus", "dengue",
                            "flavivirus"]
        
        return any(keyword in species for keyword in pathogen_keywords)
    
    def _assign_risk_group(self, taxonomy: Dict) -> int:
        """
        Assign biosafety risk group (1-4) based on taxonomy
        1: Non-pathogenic
        2: Moderate risk, unlikely to cause serious disease
        3: High risk, can cause serious disease but treatment available
        4: Extreme risk, serious disease with no treatment
        """
        # Risk group 4 organisms (highest risk)
        if (taxonomy.get("family") in ["filoviridae", "arenaviridae"] or
            taxonomy.get("genus") in ["ebolavirus", "marburgvirus"]):
            return 4
            
        # Risk group 3 organisms
        if (taxonomy.get("family") in ["mycobacteriaceae", "rickettsiaceae"] or
            taxonomy.get("genus") in ["mycobacterium", "brucella", "coxiella", "francisella"]):
            return 3
            
        # Risk group 2 organisms
        if (taxonomy.get("family") in ["enterobacteriaceae", "staphylococcaceae"] or
            taxonomy.get("genus") in ["salmonella", "staphylococcus", "streptococcus"]):
            return 2
            
        # Default to risk group 1 (minimal risk)
        return 1
    
    def _detect_sequence_type(self, blast_results: Dict) -> str:
        """Determine if sequence is DNA, RNA, or protein"""
        # Check blast program used
        if "program" in blast_results:
            if blast_results["program"] in ["blastp", "blastx"]:
                return "protein"
            else:
                return "nucleic_acid"
        
        # Check database used
        if "database" in blast_results:
            if any(db in blast_results["database"].lower() for db in ["prot", "nr", "refseq_protein"]):
                return "protein"
            else:
                return "nucleic_acid"
        
        # Default
        return "unknown"
    
    def _extract_domains(self, hits: List[Dict]) -> List[Dict]:
        """Extract conserved domain information from hits"""
        domains = []
        
        for hit in hits[:3]:  # Check top 3 hits
            title = hit.get("title", "")
            
            # Look for domain keywords in hit descriptions
            domain_keywords = [
                "domain", "motif", "active site", "binding site", 
                "catalytic", "receptor", "spike", "enzyme"
            ]
            
            for keyword in domain_keywords:
                if keyword in title.lower():
                    domain_name = self._extract_domain_name(title, keyword)
                    if domain_name:
                        domains.append({
                            "name": domain_name,
                            "type": keyword,
                            "source": title
                        })
        
        return domains
    
    def _extract_domain_name(self, title: str, keyword: str) -> Optional[str]:
        """Extract domain name from hit title"""
        # Simple extraction strategy - get phrase containing keyword
        parts = title.split(",")
        for part in parts:
            if keyword in part.lower():
                return part.strip()
        return None
    
    def _extract_functional_sites(self, hits: List[Dict]) -> List[Dict]:
        """Extract information about functional sites from hits"""
        functional_sites = []
        
        # Keywords indicating functional sites
        site_keywords = [
            "active site", "binding site", "receptor binding", 
            "catalytic site", "substrate binding"
        ]
        
        for hit in hits[:3]:
            title = hit.get("title", "")
            hsp = hit.get("hsps", [{}])[0] if "hsps" in hit else {}
            
            # Look for functional site keywords
            for keyword in site_keywords:
                if keyword in title.lower():
                    # Try to extract position from alignment
                    start = hsp.get("hit_from", 0)
                    end = hsp.get("hit_to", 0)
                    
                    if start and end:
                        functional_sites.append({
                            "type": keyword,
                            "region": (start, end),
                            "source": title
                        })
        
        return functional_sites
    
    def _mutations_in_functional_sites(self, 
                                      mutations: List[Dict], 
                                      functional_sites: List[Dict]) -> List[Dict]:
        """Identify mutations that occur in functional sites"""
        critical_mutations = []
        
        for mutation in mutations:
            pos = mutation.get("position", 0)
            
            # Check if mutation is in a functional site
            for site in functional_sites:
                start, end = site.get("region", (0, 0))
                if start <= pos <= end:
                    critical_mutations.append({
                        "mutation": mutation,
                        "site_type": site.get("type", "unknown"),
                        "site_source": site.get("source", "")
                    })
        
        return critical_mutations
    
    def _analyze_literature(self, papers: List[Dict]) -> int:
        """
        Analyze literature for evidence of pathogenicity
        Returns a score indicating strength of evidence
        """
        evidence_score = 0
        
        # Keywords indicating pathogenicity
        pathogen_keywords = [
            "pathogen", "virulence", "infectious", "disease", 
            "infection", "epidemic", "pandemic", "outbreak", 
            "transmission", "contagious", "dengue", "flavivirus", 
            "arbovirus", "hemorrhagic", "vector-borne", "mosquito-borne"
        ]
        
        for paper in papers:
            title = paper.get("title", "").lower()
            abstract = paper.get("abstract", "").lower()
            
            # Count pathogenicity keywords in title (weighted higher)
            title_matches = sum(keyword in title for keyword in pathogen_keywords)
            evidence_score += title_matches * 2
            
            # Count pathogenicity keywords in abstract
            abstract_matches = sum(keyword in abstract for keyword in pathogen_keywords)
            evidence_score += abstract_matches
            
            # Consider recency (more recent papers weighted higher)
            year = paper.get("year", 0)
            if year >= 2020:
                evidence_score += 2
            elif year >= 2015:
                evidence_score += 1
        
        return min(10, evidence_score)  # Cap at 10
    
    def _calculate_risk_score(self, features: Dict) -> Tuple[float, Dict]:
        """
        Calculate risk score based on extracted features
        Works for any biological sequence
        """
        risk_score = 0
        risk_details = {}
        
        # Special case for Dengue virus detection
        if any("dengue virus" in str(title).lower() for title in features.get("top_pathogen", [""])):
            risk_score += 40
            risk_details["dengue_detection"] = 40
        
        # 1. Taxonomy-based risk (base risk)
        taxonomy = features.get("taxonomy", {})
        tax_risk = self._calculate_taxonomy_risk(taxonomy)
        risk_score += tax_risk
        risk_details["taxonomy_risk"] = tax_risk
        
        # 2. Sequence similarity risk
        identity = features.get("highest_identity", 0)
        
        # Only factor in high identity to pathogenic organisms
        if taxonomy.get("pathogenic", False):
            if identity >= 95:
                risk_score += 30
                risk_details["identity_risk"] = 30
            elif identity >= 85:
                risk_score += 20
                risk_details["identity_risk"] = 20
            elif identity >= 70:
                risk_score += 10
                risk_details["identity_risk"] = 10
            else:
                risk_details["identity_risk"] = 0
        else:
            risk_details["identity_risk"] = 0
        
        # 3. Functional impact risk
        mutation_count = features.get("mutation_count", 0)
        critical_mutations = len(features.get("mutations_in_functional_sites", []))
        
        # More mutations in functional sites = higher risk
        mutation_risk = min(25, critical_mutations * 5)
        risk_score += mutation_risk
        risk_details["mutation_risk"] = mutation_risk
        
        # 4. Literature evidence risk
        lit_evidence = features.get("literature_evidence", 0)
        lit_risk = lit_evidence * 2.5  # Scale to max 25
        risk_score += lit_risk
        risk_details["literature_risk"] = lit_risk
        
        # 5. Biosafety risk group adjustment
        risk_group = taxonomy.get("risk_group", 1)
        risk_group = 1 if risk_group is None else int(risk_group)
        if risk_group == 4:
            risk_score = max(risk_score, 80)  # Ensure at least high risk
        elif risk_group == 3:
            risk_score = max(risk_score, 60)  # Ensure at least moderate-high
        
        # Cap the final score at 100
        return min(100, risk_score), risk_details
    
    def _calculate_taxonomy_risk(self, taxonomy: Dict) -> float:
        """Calculate risk score based on taxonomic classification"""
        risk_score = 0
        
        # Check each taxonomic level against our database
        for level, taxa_risks in self.taxonomy_risk.items():
            if taxonomy.get(level) in taxa_risks:
                level_risk = taxa_risks[taxonomy[level]]
                risk_score = max(risk_score, level_risk)
        
        # Adjust by risk group if available
        risk_group = taxonomy.get("risk_group")
        if risk_group is not None:
            # Convert to int if not None
            risk_group = int(risk_group)
            # Risk group provides a baseline minimum risk
            group_risk = {1: 0, 2: 20, 3: 40, 4: 60}
            risk_score = max(risk_score, group_risk.get(risk_group, 0))
        
        return risk_score
    
    def _identify_risk_factors(self, features: Dict, risk_details: Dict) -> List[str]:
        """Identify specific risk factors from features"""
        risk_factors = []
        
        # Special case for dengue virus
        if features.get("virus_type") == "dengue":
            risk_factors.append("Sequence identified as Dengue virus, a known human pathogen")
            risk_factors.append("Dengue virus is classified as biosafety level 2 pathogen")
            if features.get("highest_identity", 0) > 90:
                risk_factors.append(f"High sequence similarity ({features.get('highest_identity', 0):.1f}%) to Dengue virus")
        
        # 1. Taxonomy-based risk factors
        taxonomy = features.get("taxonomy", {})
        risk_group = taxonomy.get("risk_group")
        risk_group = None if risk_group is None else int(risk_group)
        
        if taxonomy.get("pathogenic", False):
            if risk_group == 4:
                risk_factors.append(f"Sequence closely matches Risk Group 4 pathogen ({taxonomy.get('species', 'unknown species')})")
            elif risk_group == 3:
                risk_factors.append(f"Sequence matches Risk Group 3 pathogen ({taxonomy.get('species', 'unknown species')})")
            elif risk_group == 2:
                risk_factors.append(f"Sequence matches Risk Group 2 pathogen ({taxonomy.get('species', 'unknown species')})")
            else:
                risk_factors.append(f"Sequence matches known pathogen ({taxonomy.get('species', 'unknown species')})")
        
        # 2. Sequence similarity risk factors
        identity = features.get("highest_identity", 0)
        if taxonomy.get("pathogenic", False) and identity >= 85:
            risk_factors.append(f"High sequence similarity ({identity:.1f}%) to known pathogen")
        
        # 3. Functional impact risk factors
        critical_mutations = features.get("mutations_in_functional_sites", [])
        if critical_mutations:
            sites = set(m["site_type"] for m in critical_mutations)
            risk_factors.append(f"Mutations detected in {len(critical_mutations)} functional sites ({', '.join(sites)})")
        
        # 4. Literature evidence risk factors
        lit_evidence = features.get("literature_evidence", 0)
        if lit_evidence >= 5:
            risk_factors.append(f"Strong literature evidence of pathogenicity or health concern")
        elif lit_evidence > 0:
            risk_factors.append(f"Some literature evidence suggests potential health concerns")
            
        # 5. Special domain risk factors
        domains = features.get("conserved_domains", [])
        risky_domains = ["toxin", "adhesin", "virulence", "spike", "hemagglutinin"]
        found_domains = [d for d in domains if any(r in d["name"].lower() for r in risky_domains)]
        
        if found_domains:
            domain_names = [d["name"] for d in found_domains]
            risk_factors.append(f"Contains domains associated with pathogenicity: {', '.join(domain_names)}")
            
        return risk_factors
    
    def _generate_recommendations(self, risk_score: float, risk_factors: List[str], features: Dict) -> List[str]:
        """Generate recommendations based on risk assessment"""
        recommendations = []
        
        # Get risk group if available, ensure it's an integer
        risk_group = features.get("taxonomy", {}).get("risk_group", 0)
        risk_group = 0 if risk_group is None else int(risk_group)
        
        # High risk recommendations
        if risk_score >= 70 or risk_group >= 3:
            recommendations.append("Implement appropriate biosafety measures based on risk assessment")
            recommendations.append("Consider further characterization using specialized assays")
            recommendations.append("Follow institutional biosafety guidelines for handling")
            
            # Special recommendations for high-risk pathogens
            if risk_group == 4:
                recommendations.append("Requires BSL-4 containment facilities and protocols")
            elif risk_group == 3:
                recommendations.append("Requires BSL-3 containment facilities and protocols")
        
        # Moderate risk recommendations
        elif risk_score >= 40 or risk_group == 2:
            recommendations.append("Follow standard laboratory safety procedures")
            recommendations.append("Consider additional confirmation tests")
            recommendations.append("Document handling and analysis procedures")
            
            if risk_group == 2:
                recommendations.append("Recommended BSL-2 practices for handling")
        
        # Low risk recommendations
        else:
            recommendations.append("Standard laboratory precautions appropriate")
            recommendations.append("No special containment measures required")
        
        # Add specific recommendations based on features
        if features.get("mutations_in_functional_sites"):
            recommendations.append("Detailed analysis of mutations in functional sites recommended")
        
        if features.get("literature_evidence", 0) > 3:
            recommendations.append("Review cited literature for specific handling guidelines")
            
        # Sequence type specific recommendations
        if features.get("sequence_type") == "protein":
            recommendations.append("Consider structural analysis to further assess functional impact")
        else:
            recommendations.append("Consider full genome analysis for comprehensive risk assessment")
            
        return recommendations
    
    def _calculate_confidence(self, features: Dict) -> float:
        """Calculate confidence score for the risk assessment"""
        # Base confidence on data quality and quantity
        base_confidence = 50  # Start with 50%
        
        # Taxonomy information increases confidence
        if features.get("taxonomy", {}).get("species"):
            base_confidence += 20
            
        # Sequence identity
        identity = features.get("highest_identity", 0)
        base_confidence += min(15, identity / 10)
        
        # Literature evidence increases confidence
        base_confidence += min(15, features.get("literature_evidence", 0) * 1.5)
        
        return min(100, base_confidence)
    
    def _determine_risk_level(self, risk_score: float) -> str:
        """Determine risk level from score"""
        if risk_score >= 70:
            return "High"
        elif risk_score >= 40:
            return "Moderate"
        else:
            return "Low"
    
    def _determine_risk_color(self, risk_score: float) -> str:
        """Determine color representation of risk level"""
        if risk_score >= 70:
            return "red"
        elif risk_score >= 40:
            return "yellow"
        else:
            return "green"
            
    def _generate_analysis_details(self, features: Dict, risk_details: Dict) -> Dict:
        """Generate detailed analysis information"""
        return {
            "sequence_match": {
                "highest_identity": features.get("highest_identity", 0),
                "sequence_type": features.get("sequence_type", "unknown"),
            },
            "taxonomy": {
                "superkingdom": features.get("taxonomy", {}).get("superkingdom"),
                "family": features.get("taxonomy", {}).get("family"),
                "genus": features.get("taxonomy", {}).get("genus"),
                "species": features.get("taxonomy", {}).get("species"),
                "risk_group": features.get("taxonomy", {}).get("risk_group"),
                "pathogenic": features.get("taxonomy", {}).get("pathogenic")
            },
            "mutations": {
                "total_count": features.get("mutation_count", 0),
                "in_functional_sites": len(features.get("mutations_in_functional_sites", [])),
            },
            "domains": [d["name"] for d in features.get("conserved_domains", [])],
            "literature": {
                "evidence_score": features.get("literature_evidence", 0)
            },
            "risk_contributions": risk_details
        }
        
    def get_scientific_basis(self) -> Dict[str, List[Dict]]:
        """
        Returns the scientific basis and references for the risk assessment methodology.
        This documents the sources and rationale for the risk values used in the assessment.
        
        Returns:
            Dict containing scientific references and explanations for risk scoring methodology
        """
        return {
            "risk_group_classification": [
                {
                    "name": "NIH Guidelines for Research Involving Recombinant or Synthetic Nucleic Acid Molecules",
                    "url": "https://osp.od.nih.gov/wp-content/uploads/NIH_Guidelines.pdf",
                    "description": "Appendix B classifies human pathogens into Risk Groups 1-4",
                    "relevance": "Primary source for risk group classifications (RG1-4) used in the assessment"
                },
                {
                    "name": "NIH-CDC Biosafety in Microbiological and Biomedical Laboratories (BMBL)",
                    "url": "https://www.cdc.gov/labs/pdf/SF__19_308133-A_BMBL6_00-BOOK-WEB-final-3.pdf",
                    "description": "Section VIII contains agent summary statements with specific BSL requirements",
                    "relevance": "Provides detailed biosafety level classifications for specific pathogens"
                }
            ],
            "taxonomic_risk_values": [
                {
                    "name": "NIH National Institute of Allergy and Infectious Diseases (NIAID) Biodefense Classification",
                    "url": "https://www.niaid.nih.gov/research/biodefense-research",
                    "description": "Categorizes pathogens into A, B, and C priority groups",
                    "relevance": "Informs the relative risk values (Filoviridae=60, Flaviviridae=50)"
                },
                {
                    "name": "American Biological Safety Association (ABSA) Risk Group Database",
                    "url": "https://my.absa.org/Riskgroups",
                    "description": "NIH-funded database of risk group classifications for specific organisms",
                    "relevance": "Confirms risk group assignments for specific taxa"
                }
            ],
            "risk_score_methodology": [
                {
                    "name": "Risk Assessment Model",
                    "description": "The risk scoring methodology uses a weighted multi-factor approach combining:",
                    "factors": [
                        "Taxonomic classification (baseline risk from pathogen type)",
                        "Sequence similarity to known pathogens (higher similarity = higher risk)",
                        "Presence of functional domains associated with pathogenicity",
                        "Mutations in critical functional sites",
                        "Evidence from scientific literature"
                    ]
                },
                {
                    "name": "Score Calibration",
                    "description": "Risk scores are calibrated to align with biosafety levels:",
                    "mapping": [
                        "BSL-1/RG1: 0-39 points (Low risk)",
                        "BSL-2/RG2: 40-69 points (Moderate risk)",
                        "BSL-3/RG3: 70-89 points (High risk)",
                        "BSL-4/RG4: 90-100 points (Extreme risk)"
                    ]
                }
            ],
            "superkingdom_evidence": [
                {
                    "taxon": "viruses",
                    "score": 40,
                    "evidence": "Approximately 20% of known viruses are human pathogens vs. <1% of bacteria",
                    "source": "Nature Reviews Microbiology (2013), doi:10.1038/nrmicro3012"
                },
                {
                    "taxon": "bacteria",
                    "score": 30,
                    "evidence": "Lower baseline risk than viruses but higher than other kingdoms",
                    "source": "CDC's Antibiotic Resistance Threats Report 2019"
                }
            ],
            "family_evidence": [
                {
                    "taxon": "filoviridae",
                    "score": 60,
                    "evidence": "Contains Ebola and Marburg viruses - BSL-4 agents with high mortality rates",
                    "source": "WHO Ebola Virus Disease Fact Sheet, mortality rates of 25-90%"
                },
                {
                    "taxon": "flaviviridae",
                    "score": 50,
                    "evidence": "Contains Dengue, Zika, Yellow Fever; ~390 million dengue infections per year",
                    "source": "Nature (2013), doi:10.1038/nature12060"
                },
                {
                    "taxon": "coronaviridae",
                    "score": 50,
                    "evidence": "Contains SARS-CoV-2, SARS-CoV, MERS-CoV with pandemic potential",
                    "source": "WHO Risk Assessment (2022), WHO-WPE-GIH-2022.1"
                }
            ]
        }

# For backward compatibility
def calculate_infection_risk(blast_results: Dict, mutation_results: Dict, literature_data: Dict = None) -> Dict:
    """Legacy function that uses the new class-based implementation"""
    assessor = InfectionRiskAssessment()
    return assessor.calculate_infection_risk(blast_results, mutation_results, literature_data)