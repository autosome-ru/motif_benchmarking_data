# Motif benchmarking data
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3702150.svg)](https://doi.org/10.5281/zenodo.3702150)

This data repository contains performance measures and derived results
related to the following benchmarking study:

> Giovanna Ambrosini, Ilya Vorontsov, Dmitry Penzar, Romain Groux, Oriol Fornes,
> Daria D. Nikolaeva, Benoit Ballester, Jan Grau, Ivo Grosse, Vsevolod Makeev,
> Ivan Kulakovskiy, Philipp Bucher. Insights gained from a comprehensive
> all-against-all transcription factor binding motif benchmarking study.
> Genome Biology 2020.  

## 1. Source annotation files for experiments and matrix collections

 - remap.txt      
 - jolma_yang.txt
 - uniprobe.txt  
 - uniprobe2.txt  
 - jaspar.txt
 - hocomoco.txt
 - cisbp.txt

Format: tab-delimited table
```
Id, descriptive name, gene symbol(s)
```

The experiements or matrices are listed in the same order as 
in the results files (see below). The Id field contains the names
of the experiments and matrices exactly as they are presented in 
the results file. These Ids are typcally derived from the
original sources but may have been modified such as to become
syntactically valid names in R (replacement of illegal 
characters by ".").

Descriptive names are meant to be used in Figures. In `remap.txt`,
`jolma_yang.txt` and `hocomoco.txt`, they are identical to the Ids.
In `jaspar.txt`, Ids consist only of the accession numbers whereas the
descriptive names include additional information about the
target TF(s), e.g.
```
MA0002.2 -> MA0002.2_RUNX1
```

In `cisbp.txt` and `uniprobe.txt`, the descriptive names represent
unmodified or minimally edited versions of Ids, e.g. 
```
M00301_2.00_NKX2.8_Homo_sapiens -> M00301_2.00_NKX2-8_Homo_sapiens
```

The 3rd field contains a comma-separated list of gene symbols. If the 
target TF of the matrix or the experiment is from another
species (usually mouse) the gene symbol corresponds to the human
orthologs. All gene symbols given in the 3rd field are verified
symbols approved by HGNC. In the extremely rare case, where an
experiment or a matrix could not be mapped to an approved gene
symbol, this field is filled with a dash. 

In most cases there is only one gene symbol in this field. The
exceptions are JASPAR matrices for multimeric TFs, e.g.

```
MA0099.3 -> MA0099.3_FOS::JUN -> FOS,JUN
MA0513.1 -> MA0513.1_SMAD2::SMAD3::SMAD4 -> SMAD2,SMAD3,SMAD4
```

There are also 4 exceptions in remap.txt, all relating to
a subfamily-specific antibody recognizing RXRA, RXRB and RXRG. 

Benchmarking with UniPROBE PBMs was done twice, first with a manually
curated compilation of PBMs for human wild-type DNA binding domains
(uniprobe.txt) and a second time with a larger collection including
experiments done with mouse TFs as well (uniprobe2.txt). Note that a 
few mouse TFs could not be mapped to human orthologs. In this case,
the third field contains a dash ("-").

## 2. Gene index tables 

 - genes2remap.txt     
 - genes2jolma_yang.txt 
 - genes2uniprobe.txt  
 - genes2uniprobe2.txt  
 - genes2jaspar.txt
 - genes2hocomoco.txt
 - genes2cisbp.txt


Format: tab-delimited 
```
gene symbol, Id, descriptive name
```

These files are information-wise equivalent to the annotation
files described above. They are given for convenience to simplify
or avoid certain operations in R like parsing comma-seperated
lists of gene symbols or dealing with lists of variable-length
vectors. The main differences to the basic annotation files are: 

 - The gene symbol appears in the first field.
 - The files are sorted by gene symbol. 
 - There are multiple lines for matrices or experiments
   associated with multiple genes.
 - Experiments done with mouse TFs that could not be mapped 
   to humman orthologs are missing.

## 3. Additonal annotation files (all tab-delimited tables)

 - hgnc_gene_symbols.txt 

Format: `gene symbol, gene description`

This file has been downloaded from HGNC on Sep 26, 2019. It is restricted 
to protein-coding genes. 

 - genes2TFclass.txt

Format: `gene symbol, TFclass family number, TFclass family name`

This file has been manually edited based in information provided in
[TFclass](http://www.edgar-wingender.de/huTF_classification.html).

## 4. Benchmark files

 - remap_jaspar_roc.txt    
 - remap_hocomoco_roc.txt  
 - remap_cisbp_roc.txt     

 - jolma_yang_jaspar_roc10.txt
 - jolma_yang_hocomoco_roc10.txt
 - jolma_yang_cisbp_roc10.txt

 - jolma_yang_jaspar_roc50.txt
 - jolma_yang_hocomoco_roc50.txt
 - jolma_yang_cisbp_roc50.txt

 - jolma_yang_shuf_jaspar_roc10.txt 
 - jolma_yang_shuf_hocomoco_roc10.txt
 - jolma_yang_shuf_cisbp_roc10.txt 

 - jolma_yang_shuf_jaspar_roc50.txt 
 - jolma_yang_shuf_hocomoco_100k_roc50.txt
 - jolma_yang_shuf_cisbp_roc50.txt 

 - uniprobe_jaspar_cor.txt
 - uniprobe_hocomoco_cor.txt
 - uniprobe_cisbp_cor.txt     

 - uniprobe2_jaspar_cor.txt
 - uniprobe2_hocomoco_cor.txt
 - uniprobe2_cisbp_cor.txt     

Format: tab-delimited numerical matrix with row and column names
   
Rows correspond to experiments, columns to matrices. Experiment and 
matrix names are identical to the Ids given in the corresponding 
annotations files. The same holds for the order in which they appear. 

## 5. Filtered benchmark files:

The subdirectory `./filtered contains` filtered benchmark files:

 - ./filtered/remap_jaspar_roc.txt
 - ./filtered/remap_hocomoco_roc.txt
 - ./filtered/remap_cisbp_roc.txt

 - ./filtered/jolma_yang_jaspar_roc10.txt
 - ./filtered/jolma_yang_hocomoco_roc10.txt
 - ./filtered/jolma_yang_cisbp_roc10.txt

 - ./filtered/uniprobe2_jaspar_cor.txt
 - ./filtered/uniprobe2_hocomoco_cor.txt
 - ./filtered/uniprobe2_cisbp_cor.txt

 - ./filtered/all_jaspar.txt
 - ./filtered/all_hocomoco.txt
 - ./filtered/all_cisbp.txt

These files provide benchmarking results for subsets of experiments, 
which reached a high performance value (AUC ROC > 0.75 or Pearson 
correlation > 0.35) with at least one matrix.

The files named `all_*.txt` contain peformance values for all types of 
experiment (remap, jolma_yang, uniprobe2) which reached a high 
performance value with at least one of the matrix collections (jaspar,
hocomoco, cisbp). Resuts for different matrix collections are
nevertheless kept separate to keep file sizes reasonably small.

Corresponding annotation files are provided in the root directory:

 - remap_all_filtered.txt
 - jolma_yang_all_filtered.txt
 - uniprobe2_all_filtered.txt

 - all_jaspar_filtered.txt
 - all_hocomoco_filtered.txt
 - all_cisbp_filtered.txt
 - all_filtered.txt

 - genes2all_jaspar_filtered.txt
 - genes2all_hocomoco_filtered.txt
 - genes2all_cisbp_filtered.txt
 
## 6. Subdirectory `./scripts`

 - best_mat.R 

   R script that reproduces the lists of best-performing matrices
   provided in the `./analysis` subdirectory from the benchmarking and
   annotation files provided here.  

## 7. Subdirectory `./analysis`

Contains lists of best matrices for experiments and genes:
### ReMap
best matrices for experiments | best matrices for genes
:-----------------------------|:-----------------------
best4exp_remap_jaspar.txt   | best4gene_remap_jaspar.txt
best4exp_remap_hocomoco.txt | best4gene_remap_hocomoco.txt
best4exp_remap_cisbp.txt    | best4gene_remap_cisbp.txt
best4exp_remap_all.txt      | best4gene_remap_all.txt

### HT-SELEX 10%
best matrices for experiments      | best matrices for genes             
:----------------------------------|:-----------------------
best4exp_jolma_yang10_jaspar.txt   | best4gene_jolma_yang10_jaspar.txt
best4exp_jolma_yang10_hocomoco.txt | best4gene_jolma_yang10_hocomoco.txt
best4exp_jolma_yang10_cisbp.txt    | best4gene_jolma_yang10_cisbp.txt
best4exp_jolma_yang10_all.txt      | best4gene_jolma_yang10_all.txt

### HT-SELEX 50%
best matrices for experiments      | best matrices for genes
:----------------------------------|:-----------------------
best4exp_jolma_yang50_jaspar.txt   | best4gene_jolma_yang50_jaspar.txt
best4exp_jolma_yang50_hocomoco.txt | best4gene_jolma_yang50_hocomoco.txt
best4exp_jolma_yang50_cisbp.txt    | best4gene_jolma_yang50_cisbp.txt
best4exp_jolma_yang50_all.txt      | best4gene_jolma_yang50_all.txt

### Uniprobe
best matrices for experiments   | best matrices for genes
:-------------------------------|:-----------------------
best4exp_uniprobe2_jaspar.txt   | best4gene_uniprobe2_jaspar.txt
best4exp_uniprobe2_hocomoco.txt | best4gene_uniprobe2_hocomoco.txt
best4exp_uniprobe2_cisbp.txt    | best4gene_uniprobe2_cisbp.txt
best4exp_uniprobe2_all.txt      | best4gene_uniprobe2_all.txt

### All benchmark types
best matrices for experiments | best matrices for genes
:-----------------------------|:-----------------------
best4exp_all_jaspar.txt       | best4gene_all_jaspar.txt
best4exp_all_hocomoco.txt     | best4gene_all_hocomoco.txt
best4exp_all_cisbp.txt        | best4gene_all_cisbp.txt
best4exp_all_all.txt          | best4gene_all_all.txt

### All benchmark types (filtered)
best matrices for experiments      | best matrices for genes
:----------------------------------|:-----------------------
best4exp_all_jaspar_filtered.txt   | best4gene_all_jaspar_filtered.txt
best4exp_all_hocomoco_filtered.txt | best4gene_all_hocomoco_filtered.txt
best4exp_all_cisbp_filtered.txt    | best4gene_all_cisbp_filtered.txt
best4exp_all_all_filtered.txt      | best4gene_all_all_filtered.txt

An aggregate performance score is used to choose the best 
performing matrix for a gene. This score is computed as follows. 
Performance measures (AUC ROC or Pearson correlation) for
individual experiments are first converted into ranks. The
aggregate performance score for a given matrix is then defined
as the geometric mean over its rank score for all experiments. 

Format of `best4exp_*.txt`: tab-delimited file with 5 fields:

 - experiment descriptive name
 - matrix descriptive name
 - gene(s) assigned to experiment (possibly comma-seprated list)
 - gene(s) assigned to matrix (possibly comma-seperated list)
 - Performance value (ROC-AUC or correlation coefficient) 

Format of `best4gene_*.txt`: tab-delimited file with 4 fields:

 - gene
 - matrix name
 - aggregate performance score
 - gene(s) assigned to matrix (possibly comma-seperated list)

Note that a few experiments from ReMap as well as some matrices 
from JASPAR are linked to multiple genes.  

## 8. Supplementary scripts
Supplementary scripts used to analyse these data are provided in the repository:
https://github.com/autosome-ru/motif_benchmarking_paper
