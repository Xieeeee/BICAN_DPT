
## 250203
This directory contains in-house scripts and functions to process and analyze Droplet Paired-Tag data.

### upstoools
To handle the 1st round barcode (pre-index) for cellranger-atac pipeline, we extracted the barcode from the beginning of R2 fastq sequence and append it to read name, which can be kept during the process of alignment. After the alignment, the 1st round barcode will be added back to the cellular barcode field in bam files (we created a new fields, named BB). We used the modified barcode to identify explict collision and downstream analysis

**08.proc_10X_preIndex_DNA.sh**: Shell script for end-to-end process of preidexed library.
**upstoools**: Custome C++ script for pre-index extraction, and fastq trimming.

To setup upstools, run:
```
cd upstoools
sh make.sh
```

