#!/usr/bin/env python

"""
Author: Amanda Dougherty
Date: October 2025
Description: Deduplication program that removes PCR duplicates from a sorted SAM file based on UMI, strand, chromosome, and adjusted 5' start position.
Usage: python ./dougherty_deduper.py -f <input_sorted_sam_file> -o <output_deduped_sam_file> -u <umi_file> -h <help_message>
Requires: bioinfo.py module with parse_line, stranded, and adjust_pos functions; known list of umis; samtools installed and available in PATH.
"""

import argparse
import subprocess
from bioinfo import parse_line, stranded, adjust_pos

def get_args():
    parser = argparse.ArgumentParser(description="Deduplication program for a SAM file by UMI, strand, chromosome, and adjusted 5' start position, after sorting with samtools by chromosome. Assumes single-end reads.")
    parser.add_argument("-f", help="Designates absolute file path to sam file = input_sam", type=str, required =True)
    parser.add_argument("-o", help="Designates output for deduped file = deduped_sam", type=str, required=True)
    parser.add_argument("-s", help="Designates output file for sorted SAM = sorted_sam", type=str, required=True)
    parser.add_argument("-u", help="Designates file containing the list of UMIs", type=str, required=True)

    return parser.parse_args()

# A function that stores unique identifiers against which to compare each line's UMI
def load_umis(umi_file:str) -> set:
    with open(umi_file, "r") as uf:
        return set(line.strip() for line in uf)

# A function to call samtools-sort within the script 
def sort_by_chromo(input_sam:str, output_sam:str)-> str:
    cmd = ["samtools", "sort", "-O", "SAM", "-o", output_sam, input_sam]
    subprocess.run(cmd, check=True)
    print(f"Sorted SAM file created at {output_sam}")
    return output_sam  # <--- return it

# Main deduper function
def main():
    args=get_args()
    input_sam=args.f 
    deduped_sam=args.o
    umis=args.u
    known_umis = load_umis(umis)
    sorted_sam = args.s

    # Option to sort the SAM file
    sorted_sam = sort_by_chromo(input_sam, sorted_sam)
    # print(f"Sorted SAM file created at {sorted_sam}")

    # Keeping track of unknown UMIs, current line and chromosome, as well as counts
    # for the total number of reads, the unique reads we are writing out, as well as the 
    # header count and total number of chromosomes. 

    invalid_umi_count = 0
    seen = set()
    current_chrom = None
    total_reads = 0
    unique_reads = 0
    header_count = 0
    num_chroms = {}

    # Read in the sorted sam file and open the output file for writing
    with open (sorted_sam, 'r') as infile, open (deduped_sam, 'w') as outfile: 

        for line in infile:
            # Write out and count number of headers
            if line.startswith('@'):
                header_count += 1
                outfile.write(line)
                continue
            
            # Add to total read count for every line that isn't a header
            total_reads += 1

            # Collect line metadata
            umi, flag, chrom, pos, cig, strand, og_line = parse_line(line)

            # Track unique chromosomes
            if chrom not in num_chroms:
                num_chroms[chrom] = 0
        
            # Reset for every chromosome and clear
            if current_chrom != chrom:
                seen.clear()
                current_chrom = chrom 

            # Skip and track invalid UMIs
            if umi not in known_umis:
                invalid_umi_count += 1
                continue
            # Get true 5' position using adjust_pos() and generate key
            adj_pos = adjust_pos(pos, cig, flag)
            key = (umi, chrom, adj_pos, strand)

            # Write out the line only if it is first occurance; update unique read count and the number of each unique chromosome
            if key not in seen: 
                seen.add(key)
                outfile.write(og_line)
                unique_reads += 1
                num_chroms[chrom] += 1
 
    # Print out confirmation and generate summary data

    print(f"Skipped {invalid_umi_count} reads with unknown UMIs")
    print(f"Deduplication complete. Check {deduped_sam} for results.")

    print("\n===== Deduplication Summary =====")
    print(f"Input SAM: {input_sam}")
    print(f"Output SAM: {deduped_sam}")
    print(f"Headers: {header_count}")
    print(f"Total reads processed: {total_reads:,}")
    print(f"Unique reads written: {unique_reads:,}")
    print(f"Duplicates removed: {total_reads - invalid_umi_count - unique_reads:,}")
    print(f"Invalid UMIs skipped: {invalid_umi_count:,}")
    print("\nReads per chromosome:")
    for chrom, count in num_chroms.items():
        print(f"{chrom}\t{count}")
    print("=================================\n")
   
    
if __name__ == "__main__":
    main()

