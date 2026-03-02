# Deduper 
A script to remove PCR duplicates from SAM files.

This script removes PCR duplicates from **coordinate-sorted** SAM alignments using:
**UMI + strand + chromosome + adjusted 5' start position**.

## Requirements
- Python 3
- `samtools` available in PATH (only needed if you do NOT use `--no-sort`)
- `bioinfo.py` providing:
  - `parse_line(line)`
  - `adjust_pos(pos, cigar, flag)`
  - (optionally) `stranded(flag)` if used inside `parse_line`
[bioinfo.py](bioinfo.py)
- a file with a list of known UMIs [example UMI file](STL96.txt)

## Usage
[Deduplication python script](dougherty_deduper.py)         
[Slurm .sh script](dougherty_deduper.sh)

### If your SAM is NOT sorted
This will run `samtools sort` internally:

```bash
python dougherty_deduper.py \
  -f input.sam \
  -o deduped.sam \
  -u <umi.file> \
  -s sorted.sam
  ```

### If your SAM is sorted: 

This will skip samtools sorting step:

```bash
python ./dougherty_deduper.py \
    -f sorted.sam \
    -o deduped.sam \
    -u <umi.file> \
    --no-sort
```

# To generate the stats file, add the -- stats flag and specify the output path for the stats file:
``` --stats /gpfs/projects/bgmp/amdo/bioinfo/Bi624/Deduper/output/deduper_stats.tsv ```

[Example --stats output](deduper_stats.tsv)

### Test files:
[unit_tests](/Users/amandadougherty/bioinfo/Bi624/deduper/Deduper-amdough/unit_tests)





