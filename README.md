####Pipeline for mitochondrial sequencing reads. Performs quality trimming and filtering, de novo assembly, and gene annotations in GFF3 and FASTA format.

######Prerequisites:  
[Python 2.7](https://www.python.org/download/releases/2.7/)  
[CAP3](http://seq.cs.iastate.edu/cap3.html)  
[ABySS](http://www.bcgsc.ca/platform/bioinfo/software/abyss)   
[HMMER](http://hmmer.org/)  
[FASTX_toolkit](http://hannonlab.cshl.edu/fastx_toolkit/)  
 

####Add programs to path  
On unix systems, check to see if programs are in path:  

`whereis cap3 abyss-pe fastq_quality_filter hmmbuild`

If any are missing, edit environment file with `nano ~/.bashrc` and add an entry to the directory containing the program  

`export PATH=/path/to/program:$PATH`  

####Clone repository

`git clone https://github.com/n-long/mito_assembly_annotation.git`

####Create index for mitochondrial gene profiles (default is formatted for species closely related to Anolis carolinensis

`hmmpress mito_bank.hmm`

####Alternatively, create your own profile from fasta sequences, one gene per file (can be single-gene sequence or multi-fasta formatted alignment of single gene from multiple species)

`find . -maxdepth 1 -name "*.fa*" -exec hmmbuild {}_profile.hmm {} \; && cat *profile.hmm > mito_bank.hmm`

####Run mito_anno.py in directory containing de-multiplexed, barcode-free fastq sequence reads 

`python mito_anno.py`

####Check directory for gene annotation in FASTA and GFF format (output name matches input)

####Parallelized version

Coming soon!
