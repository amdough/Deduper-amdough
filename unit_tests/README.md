# Deduplication Test Files

### File [test1.sam](test1.sam)

This test SAM file contains a small subset of reads designed to test a deduplication script can accurately identify PCR duplicates based on repeated genomic positions (chromosome and start position), strand, and UMI. Each line represents one aligned read. 

### Cases: 

1. Lines 1 & 2: duplicates
    * both map to the same chromosome (1) at same position
    * both have same UMI
    --> one is removed

2. Line 3: unique
3. Line 4: unique
4. Lines 5 & 6: duplicates
    * both map to same chromosome (2) at same position
    * same UMI
    --> one is removed
5. Lines 7 & 8: duplicates
    * both map to same chromosome (X) at same position
    * same UMI
    --> one is removed
6. Lines 8 & 9: duplicates
    * both mape to same chromosome (3) and have different start positions in POS field but because of soft clipping have same start position in actuality
    * same UMI
    --> one is removed 

After deduplication, only 6 lines should remain. [deduped](test1_deduped.sam)