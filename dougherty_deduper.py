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
import tempfile
import os

def get_args():
    parser = argparse.ArgumentParser(description="Deduplication program for a SAM file by UMI, strand, chromosome, and adjusted 5' start position, after sorting with samtools by chromosome. Assumes single-end reads.")
    parser.add_argument("-f", help="Designates absolute file path to sam file", type=str, required =True)
    parser.add_argument("-o", help="Designates output for deduped file", type=str, required=True)
    parser.add_argument("-s", help="Designates output file for sorted SAM", type=str, required=True)
    parser.add_argument("-u", help="Designates file containing the list of UMIs", type=str, required=True)

    return parser.parse_args()


def load_umis(umi_file:str) -> set:
    with open(umi_file, "r") as uf:
        return set(line.strip() for line in uf)

def sort_by_chromo(input_sam:str, output_sam:str)-> str:
    cmd = ["samtools", "sort", "-O", "SAM", "-o", output_sam, input_sam]
    subprocess.run(cmd, check=True)
    print(f"Sorted SAM file created at {output_sam}")
    return output_sam  # <--- return it


def process_chromo(chrom, reads, seen, known_umis, outfile):
    # seen.clear() --> should this be here? No, we want to keep seen across chromosomes right? 

    for line in reads: 
        umi, flag, chrom, pos, cig, strand, line = parse_line(line)
        if umi not in known_umis:
            continue

        adj_pos = adjust_pos(pos, cig, flag)
        key = (umi, flag, adj_pos, strand)

        if key not in seen: 
            seen[key] = True
            outfile.write(line)
    return

def main():
    args=get_args()
    f=args.f 
    out=args.o
    umi=args.u
    known_umis = load_umis(umi)
    output_sam = args.s
    sorted_sam = sort_by_chromo(f, output_sam)

    print(f"Sorted SAM file created at {sorted_sam}")

    with open (sorted_sam, 'r') as infile, open (out, 'w') as outfile: 
        current_chrom = None
        chrom_reads = []
        seen = {}


        for line in infile:
            if line.startswith('@'):
                outfile.write(line)
                continue
            else:   
                umi, flag, chrom, pos, cig, strand, line = parse_line(line)
                if current_chrom is None:
                    current_chrom = chrom
                # if chromosome changes, process new 
                if chrom != current_chrom: 
                    process_chromo(current_chrom, chrom_reads, seen, known_umis, outfile)
                    seen = {}

                    chrom_reads = []
                    current_chrom = chrom

                chrom_reads.append(line)
        # Process the last chromosome
            if current_chrom and chrom_reads:
                process_chromo(current_chrom, chrom_reads, seen, known_umis, outfile)

    os.remove(sorted_sam)

    print("Deduplication complete. Check output file for results.")
if __name__ == "__main__":
    main()
