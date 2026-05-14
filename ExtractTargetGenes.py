#!/usr/bin/env python3

import re
from collections import defaultdict
from Bio import SeqIO, SeqRecord
import gffutils

# Target genes
TARGET_GENES = {'UL133', 'UL135', 'UL136', 'UL138'}

def read_gff3(gff3_file):
    """Parse GFF3 file and organize features by type and gene name"""
    features = defaultdict(list)
    cds_features = []
    target_gene_features = defaultdict(list)
    
    # Create a database from GFF3
    db = gffutils.from_gff3(gff3_file)
    
    for feature in db.all_features():
        if feature.featuretype == 'CDS':
            cds_features.append(feature)
        
        # Check if this is one of our target genes
        if 'gene' in feature.attributes:
            gene_name = feature.attributes['gene'][0]
            if gene_name in TARGET_GENES:
                target_gene_features[gene_name].append(feature)
    
    return cds_features, target_gene_features, db

def extract_sequences(fasta_file, gff3_file):
    """Extract CDS and target gene sequences from aligned FASTA"""
    
    # Read the reference genome sequence
    sequences = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}
    
    # Parse GFF3
    cds_features, target_genes, db = read_gff3(gff3_file)
    
    # Extract all CDS sequences
    cds_sequences = []
    for i, feature in enumerate(cds_features, 1):
        seq_id = feature.seqname
        start = feature.start - 1  # GFF3 is 1-based
        end = feature.end
        strand = feature.strand
        
        if seq_id in sequences:
            seq = sequences[seq_id][start:end]
            if strand == '-':
                seq = reverse_complement(seq)
            
            feature_id = feature.attributes.get('ID', [f'CDS_{i}'])[0]
            cds_sequences.append(SeqRecord.SeqRecord(
                seq=seq,
                id=feature_id,
                description=f"CDS {i} - Position {start+1}:{end}"
            ))
    
    # Extract target gene sequences (maintaining indels from alignment)
    target_sequences = {}
    for gene_name in TARGET_GENES:
        if gene_name in target_genes:
            gene_features = target_genes[gene_name]
            # Extract the full gene region (includes introns if present)
            if gene_features:
                feature = gene_features[0]
                seq_id = feature.seqname
                start = feature.start - 1
                end = feature.end
                strand = feature.strand
                
                if seq_id in sequences:
                    seq = sequences[seq_id][start:end]
                    if strand == '-':
                        seq = reverse_complement(seq)
                    
                    target_sequences[gene_name] = SeqRecord.SeqRecord(
                        seq=seq,
                        id=gene_name,
                        description=f"Gene {gene_name} - Position {start+1}:{end}"
                    )
    
    return cds_sequences, target_sequences

def reverse_complement(seq):
    """Get reverse complement of DNA sequence"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'a': 't', 't': 'a', 'g': 'c', 'c': 'g', '-': '-'}
    return ''.join(complement.get(base, base) for base in reversed(seq))

def main():
    # File paths - adjust as needed
    fasta_file = "lab7/Full_cmv_259sequences_Aln.fasta"
    gff3_file = "lab7/hcmv.sorted.gff3"
    
    print("Extracting sequences...")
    cds_sequences, target_sequences = extract_sequences(fasta_file, gff3_file)
    
    # Write CDS sequences
    cds_output = "lab7/all_CDS_extracted.fasta"
    SeqIO.write(cds_sequences, cds_output, "fasta")
    print(f"Extracted {len(cds_sequences)} CDS regions to {cds_output}")
    
    # Write target gene sequences
    target_output = "lab7/target_genes_UL133_UL135_UL136_UL138.fasta"
    SeqIO.write(list(target_sequences.values()), target_output, "fasta")
    print(f"Extracted {len(target_sequences)} target genes to {target_output}")
    
    # Print summary
    print("\n--- Extraction Summary ---")
    print(f"Total CDS regions extracted: {len(cds_sequences)}")
    print(f"Target genes found: {len(target_sequences)}")
    for gene_name, seq_record in target_sequences.items():
        print(f"  - {gene_name}: {len(seq_record.seq)} bp")

if __name__ == "__main__":
    main()