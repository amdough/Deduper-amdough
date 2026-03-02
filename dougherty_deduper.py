#!/usr/bin/env python

"""
Author: Amanda Dougherty
Date: October 2025
Updated: February 2026
Description: Deduplication program that removes PCR duplicates from a sorted SAM file based on UMI, strand, chromosome, and adjusted 5' start position.
Usage: 

    If input SAM is NOT sorted: 
python dougherty_deduper.py -f input.sam -o deduped.sam -u <umi_file> -s sorted.sam

    If input SAM is already sorted: (or if test file lacks @HD header line, which is required for samtools sorting and it fails)
python dougherty_deduper.py -f input_sorted.sam -o deduped.sam -u <umi_file> --no-sort

    
Requires:

- bioinfo.py module with parse_line, stranded, and adjust_pos functions
- known list of umis -
- pre-sorted SAM file OR samtools installed and available in PATH

"""

import argparse
from ast import Dict
import os
import os
import subprocess
import tempfile
import tempfile
import sys
from typing import Optional
from bioinfo import parse_line, stranded, adjust_pos

def get_args():
    parser = argparse.ArgumentParser(description="Deduplication program for a SAM file by UMI, strand, chromosome, and adjusted 5' start position, after sorting with samtools by chromosome. Assumes single-end reads.")
    parser.add_argument("-f", help="Input SAM file (sorted or unsorted depending on --no-sort)", type=str, required =True)
    parser.add_argument("-o", help="Output SAM file (deduplicated). Headers are preserved.", type=str, required=True)
    parser.add_argument("-u", help="File containing list of UMIs (format: one per line)", type=str, required=True)
    # Sorting related arguments
    parser.add_argument("-s", help="OPTIONAL: path to write sorted SAM produced from samtools. Ignored if --no-sort is set", type=str, required=False, default = None)
    parser.add_argument("--no-sort", help = "Skip samtools sorting step (use if input SAM is already sorted. Note: input SAM must have @HD header line for samtools sorting to work, so if your test file lacks this, use this flag to skip sorting)", action="store_true")
    
    # Optional output stats file
    parser.add_argument("--stats", help="OPTIONAL: path to write deduplication summary statistics (text file)", type=str, required=False, default=None)
    return parser.parse_args()

# A function that stores unique identifiers against which to compare each line's UMI
def load_umis(umi_file:str) -> set:
    with open(umi_file, "r") as uf:
        return set(line.strip() for line in uf if line.strip())  # Skip empty lines

def ensure_samtools_available() -> None:
    """Exit with an error if samtools not available."""
    try:
        subprocess.run(["samtools", "--version"], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except Exception:
        sys.stderr.write("ERROR: samtools is not available in PATH. Load module or install samtools.\n")
        sys.exit(1)
        
# An OPTIONAL function to call samtools-sort within the script 
def sort_sam(input_sam:str, output_sam:str)-> str:
    """Sort the input SAM file by chromosome using samtools and write to output_sam. Returns path to sorted SAM."""
    ensure_samtools_available()
    cmd = ["samtools", "sort", "-O", "SAM", "-o", output_sam, input_sam]
    try: 
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"ERROR: samtools sort failed with error: {e}\n")
        sys.exit(1)
        raise e
    print(f"Sorted SAM file created at {output_sam}")
    return output_sam  

def chrom_sort(chrom: str) -> tuple[int, int, str]:
    """Helper function to sort chromosomes in natural order (e.g., chr1, chr2, ..., chrX, chrY)"""
    c = chrom
    if c.lower().startswith("chr"):
        c = c[3:]
    c_upper = c.upper()
    if c_upper.isdigit():
        return (0, int(c_upper), "")
    if c_upper == "X":
        return (1, 23, "")
    if c_upper == "Y":
        return (1, 24, "")
    if c_upper in {"M", "MT"}:
        return (1, 25, "")
    return (2, 999, chrom)


def write_stats_tsv(path: str, stats: dict[str, int]) -> None:
    """Write summary stats to a TSV file."""
    with open(path, "w") as out:
        out.write("metric\tvalue\n")
        for k, v in stats.items():
            out.write(f"{k}\t{v}\n")

# Main deduper function
def main():
    args=get_args()
    input_sam=args.f 
    deduped_sam=args.o
    umis=args.u
    known_umis = load_umis(umis)
    sorted_sam = args.s
    stats_file = args.stats

  # Decide which sam file to read (sorted_sam_path)
    temp_sorted_path: Optional[str] = None

    if args.no_sort:
        sorted_sam_path = input_sam
    else:
        # If user gave -s, use it; otherwise create a temp sorted file next to output
        if args.s:
            sorted_sam_path = args.s
        else:
            # temp file in same directory as output (helps on cluster filesystems)
            out_dir = os.path.dirname(os.path.abspath(deduped_sam)) or "."
            fd, temp_sorted_path = tempfile.mkstemp(prefix="deduper_sorted_", suffix=".sam", dir=out_dir)
            os.close(fd)
            sorted_sam_path = temp_sorted_path

        sort_sam(input_sam, sorted_sam_path)

    # Keeping track of unknown UMIs, current line and chromosome, as well as counts
    # for the total number of reads, the unique reads we are writing out, as well as the 
    # header count and total number of chromosomes. 
    invalid_umi_count = 0
    total_reads = 0
    unique_reads = 0
    header_count = 0
    current_chrom = None
    seen = set()
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
    duplicates_removed = total_reads - invalid_umi_count - unique_reads

    # Print out confirmation and generate summary data

    print(f"Skipped {invalid_umi_count} reads with unknown UMIs")
    print(f"Deduplication complete. Check {deduped_sam} for results.")

    print("\n===== Deduplication Summary =====")
    print(f"Input SAM: {input_sam}")
    print(f"Output SAM: {deduped_sam}")
    print(f"Headers: {header_count}")
    print(f"Total reads processed: {total_reads:,}")
    print(f"Unique reads written: {unique_reads:,}")
    print(f"Duplicates removed: {duplicates_removed:,}")
    print(f"Invalid UMIs skipped: {invalid_umi_count:,}")
    print("\nReads per chromosome:")
    for chrom, count in sorted(num_chroms.items(),
            key=lambda x: (x[0].isdigit(), int(x[0]) if x[0].isdigit() else float('inf'))):
        print(f"{chrom}\t{count}")
    print("=================================\n")

   # Optional TSV stats
    if args.stats:
        stats = {
            "headers": header_count,
            "total_reads_processed": total_reads,
            "unique_reads_written": unique_reads,
            "duplicates_removed": duplicates_removed,
            "invalid_umis_skipped": invalid_umi_count,
        }
        write_stats_tsv(args.stats, stats)
        sys.stderr.write(f"[deduper] Stats TSV written to: {args.stats}\n")

    # Clean up temp sorted file if we created one
    if temp_sorted_path is not None:
        try:
            os.remove(temp_sorted_path)
        except OSError:
            pass
    
if __name__ == "__main__":
    main()


