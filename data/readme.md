This folder contains all the data necessary to reproduce the analzyses from the Rolling-circle RNAseq SARS-CoV-2 mutation rate project.

# Description of the files in this folder:

**all-frequencies-combined-and-lethal.rds** : The raw results from running the rolling-circle pipeline.  

  
This file can be loaded into a data.frame in R with the readRDS function.  
The data.frame organization is as follows:  
Each line corresponds to one of the 89,709 base-substitutions possible in the SARS-CoV-2 genome (3 base-substitution possible per position, 29,903 positions in the genome).  
The columns are :  

- position = Genomic position (1-based) in the reference genome
- uid = A concatenation of the position and the mutation considered which can be used as a unique identifier for each mutation.
- ref = The nucleotide in the reference genome.
- bTo = The mutated nucleotide considered.

Followed by two columns for each sample, with the naming convention: \[VARIANT]\_\[REPLICATE]\_\[PASSAGE]\_muts for the number of times this mutation was observed and \[VARIANT]\_\[REPLICATE]\_\[PASSAGE]_cov for the number of calls made at this position.

For example, at row 8,000 (uid = "8000_C_T") the column DELTA_B_4_muts (value = 2) indicates that 2 calls were made with a T instead of a C at this position in the passage #4 of the replicate B of the variant DELTA B.  

The column DELTA_B_4_cov (value = 186044) indicates that 186,044 calls were made at this position in this sample. This would translate into a frequency of C-to-T variants at position 8,000 in this sample (Delta B, passage #4) of 2/186044.
                                                                                                                                                                                                               
                                                                                                                                                                                                               
