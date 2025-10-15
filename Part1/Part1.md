# Deduper
--------------------------
## Author: Amanda Dougherty
## 10.15.2025


**Goal:** 
Create a strategy for coding a Reference Based PCR Duplicate Removal tool that, when supplied a sorted sam file of uniquely mapped reads, will remove all PCR duplicates and retain only one copy per read. 

**Problem:**   
In a single-end library, we want to remove any PCR duplicates from our data set. PCR duplicates can falsely estimate expression levels or abundance of certain genes or DNA fragments. This leads to inaccurate representations in the data which further complicates downstream analyses and reduces the reliability and robustness of results. 

Our deduper pseudocode will consider two reads duplicates if they share: 
    - same chromosome (column 3: RNAME)
    - same 5' start position (column 4: POS)
    - same strand (column 2: bitwise FLAG)
    - same UMI (parsed from end of QNAME)
    - deal with soft clipping by determining true 5' start (CIGAR string)

--> Develop your algorithm using pseudocode

**Pseudocode:**
``` 
init key_dict{}
READ in and OPEN a sorted SAM file 
WHILE not EOF:
    OPEN files for writing: <output file>
    Loop through each line of INPUT file
    IF line starts with @, write to output file
    ELSE (not a header line)
        GET column 1: QNAME --> UMI = umi
        GET column 2: FLAG --> determine strand = flag
            IF (flag & 16 is 16) flag = "-"
            ELSE flag = "+"
        GET column 3: RNAME --> determine chromosome = chrom
        GET column 4: POS --> determine 5' start position = pos
        GET column 6: CIGAR --> determine if soft-clipped = cig
            IF soft clipped, adjust POS = adj_pos (pos) 
            ELSE do nothing
        put all these into a dictionary 
            key = umi, flag, chrom, pos or adj_pos
        check to see if key already exists
        IF key not in key_dict
            key_dict.add(key)
            write line to output file
        ELSE 
            skip-a-roo

```
--> Determine high level functions  
        Description
        Function headers
        Test examples (for individual functions)
        Return statement


1. **parse_line()** - extracts all fields of interest and returns a dictionary
2. **adjust_pos()** - handles soft clipping and adjusts the 5' start position *if necessary*
3. **is_dupe()** - checks key_dict to see if a "strand" has already been added or not
4. **deduper()** - main function, ties it all together with all three helper functions 

**1**
```
def parse_line(line):
    cols = line.strip().split('\t')
    umi = cols[0].split(":")[-1]
    flag = int(cols[1])
    chrom = cols[2]
    pos = int(cols[3])
    cig = cols[5]
    strand = "-" if (flag & 16) else "+"
    return(umi, strand, chrom, pos, cig, line)
```
**2**
```
def adjust_pos(pos, cig, strand):
    if strand == "+":
        if cig.startswith(str(int)):
            adj_pos = *subtract (int) number of bases from pos*
    else:
        if cig.endswith("S"):
            adj_pos = *add bases for rev strand*
    return adj_pos
```
**3**
```
def is_dupe(umi, strand, chrom, pos, key_dict):
    key = (umi, strand, chrom, pos)
    if key in key_dict:
        return True
    else:
        key_dict[key] = True
        return False
```

**4**
```
def deduper(<input sam>, <output sam>):
    key_dict = {} # initialize dictionary of line keys
    with open(<input sam> as input), open(<output sam>, "w") as output: 
        for line in input: 
            if line.startswith("@"):
                ouput.write(line)
                continue
            umi, strand, chrom, pos, cig, og_line = parse_line(line)
            pos = adjust_pos(pos, cig, strand)
            if not is_dupe(umi, strand, chrom, pos, key_dict)
                output.write(og_line)

```



Write examples:
Include a properly formated sorted input sam file
Include a properly formated expected output sam file


