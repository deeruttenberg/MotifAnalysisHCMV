#!/usr/bin/env python3
"""
Extract CDS and specific genes from an aligned FASTA file using GFF3 annotations.
Uses the NC_006273.2 (Merlin strain) as reference and removes indels for proper alignment.
Maintains indels in all other sequences during extraction.
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

class AlignedFASTAParser:
    """Parse and manipulate aligned FASTA files with gap tracking."""
    
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
    
    def remove_gaps(self, seq: str) -> str:
        """Remove gap characters (dashes) from sequence."""
        return seq.replace('-', '')
    
    def get_gap_adjusted_coords(self, seq: str, start: int, end: int) -> Tuple[int, int]:
        """
        Convert coordinates from ungapped reference to gapped alignment positions.
        
        Args:
            seq: The aligned sequence (with gaps)
            start: Start position in ungapped sequence (0-based)
            end: End position in ungapped sequence (0-based)
        
        Returns:
            (aligned_start, aligned_end) positions in the gapped sequence
        """
        ungapped_pos = 0
        aligned_start = None
        aligned_end = None
        
        for i, char in enumerate(seq):
            if char != '-':
                if ungapped_pos == start:
                    aligned_start = i
                if ungapped_pos == end - 1:
                    aligned_end = i + 1
                    break
                ungapped_pos += 1
            else:
                # For gaps, we need to know if we're before, at, or after our target
                if ungapped_pos >= start and aligned_start is None:
                    aligned_start = i
        
        # Handle edge cases
        if aligned_start is None:
            aligned_start = len(seq)
        if aligned_end is None:
            aligned_end = len(seq)
        
        return aligned_start, aligned_end
    
    def extract_region_from_aligned(self, seq_id: str, ungapped_start: int, 
                                   ungapped_end: int, ref_ungapped_seq: str) -> str:
        """
        Extract a region from an aligned sequence using coordinates from ungapped reference.
        
        Args:
            seq_id: Sequence identifier
            ungapped_start: Start in ungapped coordinates (0-based)
            ungapped_end: End in ungapped coordinates (0-based)
            ref_ungapped_seq: The ungapped reference sequence
        
        Returns:
            Extracted aligned sequence
        """
        if seq_id not in self.sequences:
            raise ValueError(f"Sequence ID '{seq_id}' not found")
        
        aligned_seq = self.sequences[seq_id]
        
        # Convert ungapped coordinates to aligned coordinates
        # by mapping through the reference
        ungapped_pos = 0
        aligned_start = None
        aligned_end = None
        
        for i, char in enumerate(aligned_seq):
            if char != '-':
                if ungapped_pos == ungapped_start:
                    aligned_start = i
                if ungapped_pos == ungapped_end:
                    aligned_end = i
                    break
                ungapped_pos += 1
        
        # Handle edge case for end of sequence
        if aligned_end is None:
            aligned_end = len(aligned_seq)
        if aligned_start is None:
            aligned_start = aligned_end
        
        return aligned_seq[aligned_start:aligned_end]
    
    def write_fasta(self, output_file: str, sequences: Dict[str, str]):
        """Write sequences to FASTA file."""
        with open(output_file, 'w') as f:
            for seq_id, seq in sequences.items():
                f.write(f">{seq_id}\n")
                # Write sequence in 80-character lines
                for i in range(0, len(seq), 80):
                    f.write(seq[i:i+80] + '\n')

def extract_features_from_gff(fasta_parser: AlignedFASTAParser, gff_parser: GFF3Parser,
                              reference_id: str, output_prefix: str):
    """
    Extract CDS and specific genes from aligned sequences using ungapped reference coordinates.
    
    Args:
        fasta_parser: AlignedFASTAParser instance
        gff_parser: GFF3Parser instance
        reference_id: ID of reference sequence in FASTA (without gaps)
        output_prefix: Prefix for output files
    """
    
    print(f"Using reference sequence: {reference_id}")
    
    # Get and ungap the reference sequence
    if reference_id not in fasta_parser.sequences:
        raise ValueError(f"Reference sequence '{reference_id}' not found in FASTA")
    
    ref_aligned_seq = fasta_parser.sequences[reference_id]
    ref_ungapped_seq = fasta_parser.remove_gaps(ref_aligned_seq)
    print(f"Reference length (aligned): {len(ref_aligned_seq)} bp")
    print(f"Reference length (ungapped): {len(ref_ungapped_seq)} bp")
    
    # Extract CDS
    print("\nExtracting CDS features...")
    cds_features = gff_parser.get_cds_features()
    cds_sequences = {}
    
    for i, cds in enumerate(cds_features):
        feature_id = cds['attributes'].get('ID', f"CDS_{i}")
        try:
            # Get the sequence from all sequences in the alignment
            seq = fasta_parser.extract_region_from_aligned(
                reference_id, cds['start'], cds['end'], ref_ungapped_seq
            )
            cds_sequences[feature_id] = seq
        except (ValueError, IndexError) as e:
            print(f"  Warning: Could not extract {feature_id}: {e}")
    
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
                seq = fasta_parser.extract_region_from_aligned(
                    reference_id, feature['start'], feature['end'], ref_ungapped_seq
                )
                gene_sequences[feature_id] = seq
                print(f"  Extracted {feature_id} ({feature['feature_type']}): {len(seq)} bp (aligned)")
            except (ValueError, IndexError) as e:
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
    # Reference sequence ID from NC_006273.2 (Merlin strain)
    reference_id = "Europe_United Kingdom_NC_006273.2_Human_herpesvirus_5_strain_Merlin_complete_genome"
    output_prefix = "hcmv_extracted"
    
    print("=" * 70)
    print("Extracting features from aligned FASTA with GFF3 annotations")
    print("=" * 70)
    
    # Parse input files
    print("\nParsing aligned FASTA file...")
    fasta_parser = AlignedFASTAParser(fasta_file)
    print(f"Found {len(fasta_parser.sequences)} sequences")
    
    # List available sequences (first 5)
    seq_ids = list(fasta_parser.sequences.keys())[:5]
    print(f"Sample sequence IDs: {', '.join(seq_ids)}")
    
    print("\nParsing GFF3 file...")
    gff_parser = GFF3Parser(gff_file)
    print(f"Found {len(gff_parser.features)} features")
    
    # Extract features
    cds_seqs, gene_seqs = extract_features_from_gff(
        fasta_parser, gff_parser, reference_id, output_prefix
    )
    
    print("\n" + "=" * 70)
    print("Feature extraction complete!")
    print("=" * 70)
    print(f"\nOutput files:")
    print(f"  - {output_prefix}_CDS.fasta: {len(cds_seqs)} CDS features")
    print(f"  - {output_prefix}_target_genes.fasta: {len(gene_seqs)} target gene features")
    print(f"\nNote: All sequences are in aligned format (with gaps)")

if __name__ == "__main__":
    main()
