#!/usr/bin/env python

# Author: <Amanda Dougherty> <amdo@uoregon.edu>

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework. Current version
holds scripts to convert phred scores and get quality scores, calculate mapped and unmapped percentages, 
calculate the median of a list of integers, convert a fasta file to one line from multiple lines, get gc content, 
functions for removing pcr duplicates, and more. 
'''

__version__ = "0.4"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning


DNAbases = set('ATGCNatcgn')
RNAbases = set('AUGCNaucgn')


import re


def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score'''
    return ord(letter)-33 

# Loop below uses convert_phred to print out a converted string
# for i, letter in enumerate(phred_score):
#     """Iterates over each letter in the phred_score string and prints results out"""
#     print(f"{i}: {letter} - {convert_phred(letter)}")
    
def qual_score(phred_score: str) -> float:
    '''Fuction that iterates through the phred_score string, calculates quality score using convert_phred function,
    and returns the calculated average quality score of the entire string'''
    average=0
    for i, letter in enumerate(phred_score):
        average += convert_phred(phred_score[i])
    average /= len(phred_score)
    return average

if __name__ == "__main__": 
# tests for qual_score funtion
    assert qual_score("EEE") == 36
    assert qual_score("#I") == 21
    assert qual_score("EJ") == 38.5
    print("You calcluated the correct average phred score")

def validate_base_seq(seq,RNAflag=False):
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    return set(seq)<=(RNAbases if RNAflag else DNAbases)

if __name__ == "__main__": 
    assert validate_base_seq("AATAGAT") == True, "Validate base seq does not work on DNA"
    assert validate_base_seq("AAUAGAU", True) == True, "Validate base seq does not work on RNA"
    assert validate_base_seq("Hi there!") == False, "Validate base seq fails to recognize nonDNA"
    assert validate_base_seq("Hi there!", True) == False, "Validate base seq fails to recognize nonDNA"
    print("Passed DNA and RNA tests")

def init_list(lst: list, value: float=0.0) -> list:
    '''This function takes an empty list and will populate it with
    the value passed in "value". If no value is passed, initializes list
    with 101 values of 0.0; can be changed to accommodate different list sizes obvvvv'''
    if len(lst) > 0:
        raise ValueError("List must be empty.")
    lst.extend([value]*101)
    return lst


def populate_list(file: str) -> tuple[list, int]:
    '''Sums Phred quality scores for each base pair position across all reads in a FASTQ file. 
    Iterates through the file, converts quality scores from letters to numbers, and sums the scores for each base position. 
    Returns the list of summed scores and the total number of lines processed'''
    mylist = init_list([])
    # my_list: list = []
# my_list = init_list(my_list)
    with open (file, "r") as fh:
        for num_lines, line in enumerate(fh):   
            pos = 0
            if num_lines % 4 == 3: 
                line = line.strip()
                for x in line: 
                    mylist[pos]=mylist[pos]+ (convert_phred(x))
                    pos+=1
    
    return (mylist, (num_lines+1))

def get_read_length(file_name): 
    """Returns the length of the first sequencing read from a Fastq file"""
    with open (file_name, "r") as fh:
        for num_lines, line in enumerate(fh):   
            if num_lines % 4 == 1: 
                return len(line.strip())

def calc_percent_mapped(file1, file2) -> tuple[float, float]:
    """Reads into two files and line by line removes new lines, and splits the string by \t, 
        then counts both mapped and unmapped reads for each file. Returns the percentage of mapped 
        reads for each."""
    # initializing our mapped and unmapped read variables for both fw/rv
    fw_mapped = 0
    fw_unmapped = 0
    rv_mapped = 0
    rv_unmapped = 0

    # opens each file and assigns a variable

    with open(file1, "r") as fw, open(file2, "r") as rv: 
        # a loop that runs until the first position of either fw or rv line is empty
        while True:
            fw_line = fw.readline().strip("\n")
            rv_line = rv.readline().strip("\n")
            fw_line = fw_line.split("\t")
            rv_line = rv_line.split("\t")
            
            if fw_line[0] == "" or rv_line[0] == "": # if index 0 of the fw or rv line is empty, break 
                break
            if fw_line[0].startswith("ENSMUS"):
                fw_mapped += int(fw_line[1])
            else:
                fw_unmapped += int(fw_line[1])
                
            if rv_line[0].startswith("ENSMUS"):
                rv_mapped += int(rv_line[1])
            else:
                rv_unmapped += int(rv_line[1])
            # returning our percentage of mapped reads for each file 
        return fw_mapped/(fw_mapped + fw_unmapped) * 100, rv_mapped/(rv_mapped + rv_unmapped) * 100


def get_nt_length(f):
    """Calculates and prints the total number of nucleotides in a FASTQ file."""
    with open (f, "r") as fh:
        nt_tot = 0
        for i, line in enumerate(fh):   
            if i % 4 == 1: 
                nt_tot += len(line.strip())
        print(f"{nt_tot}")


def gc_content(DNA):
    '''Returns GC content of a DNA sequence as a decimal between 0 and 1.'''
    assert validate_base_seq(DNA), "String contains invalid characters - are you sure you used a DNA sequence?"
    
    DNA = DNA.upper()
    Gs = DNA.count("G")
    Cs = DNA.count("C")
    return (Gs+Cs)/len(DNA)



def calc_median(lst: list) -> float:
    '''Given a sorted list, returns the median value of the list'''
    if not lst: 
            raise ValueError("How on earth to calc a median of an empty list?")
    sorted_list = sorted(lst)
    length = len(sorted_list)
    med = length // 2
    # making sure it returns None if presented an empty list
    if length == 0:
        raise ValueError("How on earth to calc a median of an empty list?")
    # making sure it works on an even numbered list and takes the avg of the two middlest values
    if length % 2 == 0:
        return (sorted_list[med-1] + sorted_list[med])/2.0
    # if odd, just finds the median
    else:
        return sorted_list[med]

if __name__ == "__main__":
    # tests for calc_median
    assert calc_median([1,3,2]) == 2, "Failed to calculate median of odd numbered list"
    assert calc_median([4,1,2,3]) == 2.5, "Failed on an even numbered list"
    assert calc_median([8,8,8,8]) == 8, "Failed on duplicated number list"
    assert calc_median([1,2,3,4,5]) == 3, "Failed on a sorted list"


def oneline_fasta(f, out):
    """Converts a FASTA file to a one line per sequence format."""
    new_seq:str = ""

    with open (f, "r") as fh, open (f"{out}.fa", "w") as newFile:

        for line in fh:

            if line.startswith(">"):
                if new_seq != "":
                    newFile.write(f"{new_seq}\n")
                    new_seq = ""

                newFile.write(line)
            else:
                line = line.strip('\n')
                new_seq += line
                    
        newFile.write(f"{new_seq}\n")

def parse_line(line):
    cols = line.strip().split('\t')
    umi = cols[0].split(":")[-1]
    flag = int(cols[1])
    chrom = cols[2]
    pos = int(cols[3])
    cig = cols[5]
    strand = True if (flag & 16) else False
    return(umi, flag, chrom, pos, cig, strand, line)

def stranded(flag: int) -> bool: 
    """Returns True if the read is mapped to the reverse strand, False otherwise.""" 
    return (flag & 16) == 16

_CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')

def adjust_pos(pos:int, cig:str, flag:int) -> int: 
    """
    Adjusts the 5' start POS from a SAM/BAM file line based on the CIGAR string and strand.
    LOGIC: 
    pos : int
        1 based leftmost mapping POSition
    cig : str
        CIGAR string from SAM/BAM file (e.g. "10S90M" or "90M10S")
    strand : str
        Strand information, either "+" or "-"
    Returns
    -------
    adj_pos : int
        Adjusted 5' start POS based on CIGAR string and strand
    Notes:
    -----
    Reference-consuming operations in CIGAR strings include:
    M (alignment match), D (deletion), N (skipped region), = (sequence match), X (sequence mismatch)
    Non-reference-consuming operations include:
    I (insertion), S (soft clipping), H (hard clipping), P (padding)
    5' adjustment logic:
    - For "+" strand: subtract soft clipping at the start (if any) from POS
    - For "-" strand: add soft clipping at the end (if any) to POS
    """
    if not cig:
        return pos #case if no cigar string present
    
    cases = [((int)(n), case) for n, case in _CIGAR_RE.findall(cig)]
    if not cases:
        return pos #case if cigar string is malformed
    
    if not stranded(flag): # positive strand, leading soft clip applies
        leading_clip = cases[0][0] if cases and cases[0][1] == ("S") else 0
        return pos - leading_clip
    else: # negative strand, trailing soft clip and ref consuming ops apply
        length = sum(n for n, op in cases if op in ("M", "D", "N", "=", "X"))
        end = pos + length - 1  # -1 to convert to 0-based
        trailing_clip = cases[-1][0] if cases and cases[-1][1] in ("S") else 0
        return end + trailing_clip
    


if __name__ == "__main__":
    # write tests for functions above, Leslie has already populated some tests for convert_phred
    # These tests are run when you execute this file directly (instead of importing it)
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job")

    # tests for adjust_pos function
    # + strand, leading soft clip:
    # POS=100, CIGAR=5S95M → 5' start = 100 - 5 = 95
    assert adjust_pos(100, "5S95M", 0) == 95
    print("Your adjust_pos function is working! Nice job")

    # + strand, no clipping:
    # POS=100, CIGAR=100M → start = 100
    assert adjust_pos(100, "100M", 0) == 100
    print("Your adjust_pos function is working! Nice job")


    # - strand, trailing soft clip:
    # POS=100, CIGAR=95M5S → end = 100 + 95 - 1 = 194; start = 194 + 5 = 199
    assert adjust_pos(100, "95M5S", 16) == 199
    print("Your adjust_pos function is working! Nice job")

    # - strand with D and I:
    # POS=100, CIGAR=10M1I10M2D70M3S
    # ref_len = 10 + 10 + 2 + 70 = 92; end = 100 + 92 - 1 = 191; start = 191 + 3 = 194
    assert adjust_pos(100, "10M1I10M2D70M3S", 16) == 194
    print("Your adjust_pos function is working! Nice job")


