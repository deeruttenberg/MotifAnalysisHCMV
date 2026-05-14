#!/usr/bin/env python3
"""
Extract CDS and specific genes from an aligned FASTA file using GFF3 annotations.
Maintains indels (dashes) during extraction.
"""

import re
from typing import Dict, List, Tuple

class GFF3Parser:
    """Parse GFF3 annotation files."""
    
    def __init__(self, gff3_file: str):
        self.features = []
        self.parse(gff3_file)
    
    def parse(self, gff3_file: str):
        """Parse GFF3 file and store features."""
        with open(gff3_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                parts = line.split('\t')
                if len(parts) < 9:
                    continue
                
                feature = {
                    'seqname': parts[0],
                    'source': parts[1],
                    'feature_type': parts[2],
                    'start': int(parts[3]) - 1,  # Convert to 0-based
                    'end': int(parts[4]),
                    'score': parts[5],
                    'strand': parts[6],
                    'frame': parts[7],
                    'attributes': self._parse_attributes(parts[8])
                }
                self.features.append(feature)
    
    def _parse_attributes(self, attr_str: str) -> Dict:
        """Parse the attributes field of a GFF3 feature."""
        attributes = {}
        for item in attr_str.split(';'):
            if '=' in item:
                key, value = item.split('=', 1)
                attributes[key.strip()] = value.strip()
        return attributes
    
    def get_cds_features(self) -> List[Dict]:
        """Get all CDS features."""
        return [f for f in self.features if f['feature_type'] == 'CDS']
    
    def get_gene_by_name(self, gene_name: str) -> List[Dict]:
        """Get all features for a specific gene by name."""
        results = []
        for f in self.features:
            if 'gene' in f['attributes']:
                if gene_name.lower() in f['attributes']['gene'].lower():
                    results.append(f)
            # Also check in Name attribute
            if 'Name' in f['attributes']:
                if gene_name.lower() in f['attributes']['Name'].lower():
                    results.append(f)
        return results

class FASTAParser:
    """Parse and manipulate FASTA files."""
    
    def __init__(self, fasta_file: str):
        self.sequences = {}
        self.parse(fasta_file)
    
    def parse(self, fasta_file: str):
        """Parse FASTA file into dictionary."""
        current_id = None
        current_seq = []
        
        with open(fasta_file, 'r') as f:
            for line in f:
                line = line.rstrip('\n')
                if line.startswith('>'):
                    if current_id is not None:
                        self.sequences[current_id] = ''.join(current_seq)
                    current_id = line[1:].split()[0]  # Take only first part as ID
                    current_seq = []
                else:
                    current_seq.append(line)
            
            # Don't forget the last sequence
            if current_id is not None:
                self.sequences[current_id] = ''.join(current_seq)
    
    def extract_region(self, seq_id: str, start: int, end: int, 
                      maintain_indels: bool = True) -> str:
        """
        Extract a region from sequence, maintaining indels if requested.
        
        Args:
            seq_id: Sequence identifier
            start: 0-based start position
            end: End position (exclusive)
            maintain_indels: If True, preserve dashes; if False, remove them
        
        Returns:
            Extracted sequence
        """
        if seq_id not in self.sequences:
            raise ValueError(f"Sequence ID '{seq_id}' not found")
        
        seq = self.sequences[seq_id][start:end]
        
        if not maintain_indels:
            seq = seq.replace('-', '')
        
        return seq
    
    def write_fasta(self, output_file: str, sequences: Dict[str, str]):
        """Write sequences to FASTA file."""
        with open(output_file, 'w') as f:
            for seq_id, seq in sequences.items():
                f.write(f">{seq_id}\n")
                # Write sequence in 80-character lines
                for i in range(0, len(seq), 80):
                    f.write(seq[i:i+80] + '\n')

def extract_features_from_gff(fasta_parser: FASTAParser, gff_parser: GFF3Parser,
                              output_prefix: str):
    """
    Extract CDS and specific genes from aligned sequences.
    
    Args:
        fasta_parser: FASTAParser instance
        gff_parser: GFF3Parser instance
        output_prefix: Prefix for output files
    """
    
    # Get primary sequence ID (usually the reference)
    seq_id = list(fasta_parser.sequences.keys())[0]
    print(f"Using reference sequence: {seq_id}")
    
    # Extract CDS
    print("\nExtracting CDS features...")
    cds_features = gff_parser.get_cds_features()
    cds_sequences = {}
    
    for i, cds in enumerate(cds_features):
        feature_id = cds['attributes'].get('ID', f"CDS_{i}")
        try:
            seq = fasta_parser.extract_region(seq_id, cds['start'], cds['end'])
            cds_sequences[feature_id] = seq
        except ValueError as e:
            print(f"Warning: Could not extract {feature_id}: {e}")
    
    # Write CDS FASTA
    cds_output = f"{output_prefix}_CDS.fasta"
    fasta_parser.write_fasta(cds_output, cds_sequences)
    print(f"Wrote {len(cds_sequences)} CDS features to {cds_output}")
    
    # Extract specific genes
    target_genes = ['UL133', 'UL135', 'UL136', 'UL138']
    gene_sequences = {}
    
    for gene in target_genes:
        print(f"\nExtracting {gene}...")
        features = gff_parser.get_gene_by_name(gene)
        
        if not features:
            print(f"  Warning: No features found for {gene}")
            continue
        
        for j, feature in enumerate(features):
            feature_id = f"{gene}_{feature['feature_type']}"
            if j > 0:
                feature_id += f"_{j}"
            
            try:
                seq = fasta_parser.extract_region(seq_id, feature['start'], 
                                                 feature['end'])
                gene_sequences[feature_id] = seq
                print(f"  Extracted {feature_id} ({feature['feature_type']}): "
                      f"{len(seq)} bp")
            except ValueError as e:
                print(f"  Warning: Could not extract {feature_id}: {e}")
    
    # Write genes FASTA
    genes_output = f"{output_prefix}_target_genes.fasta"
    fasta_parser.write_fasta(genes_output, gene_sequences)
    print(f"\nWrote {len(gene_sequences)} features to {genes_output}")
    
    return cds_sequences, gene_sequences

def main():
    """Main execution."""
    # Configuration
    fasta_file = "Full_cmv_259sequences_Aln.fasta"
    gff_file = "hcmv.sorted.gff3"
    output_prefix = "hcmv_extracted"
    
    print("=" * 60)
    print("Extracting features from aligned FASTA with GFF3 annotations")
    print("=" * 60)
    
    # Parse input files
    print("\nParsing FASTA file...")
    fasta_parser = FASTAParser(fasta_file)
    print(f"Found {len(fasta_parser.sequences)} sequences")
    
    print("\nParsing GFF3 file...")
    gff_parser = GFF3Parser(gff_file)
    print(f"Found {len(gff_parser.features)} features")
    
    # Extract features
    cds_seqs, gene_seqs = extract_features_from_gff(fasta_parser, gff_parser,
                                                     output_prefix)
    
    print("\n" + "=" * 60)
    print("Feature extraction complete!")
    print("=" * 60)
    print(f"\nOutput files:")
    print(f"  - {output_prefix}_CDS.fasta: {len(cds_seqs)} CDS features")
    print(f"  - {output_prefix}_target_genes.fasta: {len(gene_seqs)} target gene features")

if __name__ == "__main__":
    main()
