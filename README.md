####Assemble mitochondrial genomes and retrieve gene annotations from fastq reads.  

######Prerequisites:  
[Python 2.7](https://www.python.org/download/releases/2.7/)  
[CAP3](http://seq.cs.iastate.edu/cap3.html)  
[ABySS](http://www.bcgsc.ca/platform/bioinfo/software/abyss) 
[HMMER](http://hmmer.org/)  
[FASTX_toolkit](http://hannonlab.cshl.edu/fastx_toolkit/)  
 

####Add programs to path

####Clone repository

`git clone https://github.com/n-long/mito_assembly_annotation.git`

####Create index for mitochondrial gene profiles (this one is formatted for species closely related to Anolis carolinensis

`hmmpress mito_bank.hmm`

