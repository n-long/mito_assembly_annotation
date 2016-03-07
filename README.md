####Pipeline for mitochondrial sequencing reads. Performs quality trimming and filtering, de novo assembly, and gene (protein-coding, tRNA, rRNA) annotations in GFF3 and FASTA format.

######Prerequisites:  
[Python 2.7](https://www.python.org/download/releases/2.7/)  
[CAP3](http://seq.cs.iastate.edu/cap3.html)  
[ABySS](http://www.bcgsc.ca/platform/bioinfo/software/abyss)   
[HMMER](http://hmmer.org/)  
[FASTX_toolkit](http://hannonlab.cshl.edu/fastx_toolkit/)  
 

####Add programs to path  
Check to see if programs are in path.  

For unix systems:  
`whereis cap3 abyss-pe fastq_quality_filter hmmbuild`

If any are missing, edit environment file with `nano ~/.bashrc` and add an entry to the directory containing the program  

`export PATH=/path/to/program:$PATH`  

####Clone repository

`git clone https://github.com/n-long/mito_assembly_annotation.git`

####Create index for mitochondrial gene profiles (default is formatted for species closely related to Anolis carolinensis

`hmmpress mito_bank.hmm`

####Alternatively, create your own profiles from fasta sequences, one gene per file (can be single-gene sequence or multi-fasta formatted alignment of single genes from multiple species)

`find . -maxdepth 1 -name "*.fa*" -exec hmmbuild {}_profile.hmm {} \; && cat *profile.hmm > mito_bank.hmm`

tRNA/rRNA models are precompiled* in the mitfi/ subdirectory along with the [Infernal](http://eddylab.org/infernal/) executable (no PATH adding necessary), and will work for any species without calibration.

####Run mito_anno.py in directory containing fastq sequence reads (assumes sequences have been de-multiplexed and are free of barcodes/adaptors)

`python mito_anno.py`

Output filenames will match input files. Assembled mitochondrial genomes will end in `-scaffolds.fa`, genes (including tRNA/rRNA) in `_genes.fasta`, and coordinates in `.gff`.


####*tRNA and rRNA databases come from the published MITOS (web server only) datasets

Bernt, M., Donath, A., Jühling, F., Externbrink, F., Florentz, C., Fritzsch, G., ... & Stadler, P. F. (2013). MITOS: Improved de novo metazoan mitochondrial genome annotation. Molecular phylogenetics and evolution, 69(2), 313-319.

####Parallelized version -- coming soon!

Script contains most of the code to get working, however priority is low until thesis is written.
