#!/usr/bin/python
import sys,os,subprocess,glob,shutil,getopt

def main(argv):
	num_cores = ''
	try:
		opts, args = getopt.getopt(argv,"h:c:",["cores="])
	except getopt.GetoptError:
		print 'mito_anno.py -c [number of cores]'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'mito_anno.py -c [number of cores]'
			sys.exit()
		elif opt in ("-c", "--cores"):
			num_cores = arg
if __name__ == "__main__":
	main(sys.argv[1:])
raw_input("All fastq files in current directory will be processed (*.fq or *.fastq). Press Enter to continue, or Ctrl+C / Ctrl+D to exit.")
fastq_files = sorted(glob.glob("*.f*q"), key=os.path.abspath)

def annotation():
	for element in fastq_files:
        	subprocess.call(['fastq_quality_filter', '-q', '20', '-p', '90', '-i', element, '-o', element+'.processed', '-Q33'])
        	subprocess.call(['fastq_quality_trimmer', '-t', '20', '-l', '50', '-i', element+'.processed', '-o', element+'.processed.filtered', '-Q33'])
        	filtered_fastq=sorted(glob.glob("*.filtered"), key=os.path.abspath)
	for i in range(0,len(filtered_fastq),2):
        	paired_names=filtered_fastq[i:i+2]
        	name=paired_names[0].split(".")
        	outname=paired_names[0].split("_read")
        	clean_names="../"+paired_names[0]+"\t"+"../"+paired_names[1]
        	current_path=os.getcwd()
        	better_path=os.path.join('./',outname[0])
        	if not os.path.exists(better_path):
                	os.mkdir(better_path)
        	os.chdir(better_path)
        	subprocess.call(['abyss-pe', 'name='+name[0], 'k=64', 'in='+clean_names])
        	subprocess.call(['cap3', name[0]+'-scaffolds.fa'])
        	copy=name[0]+'-scaffolds.fa.cap.contigs'
        	shutil.copy2(copy, '../')
        	os.chdir(current_path)
        	dfamtbl=outname[0]+'.tbl'
        	subprocess.call(['nhmmscan', '--dfamtblout', dfamtbl, '--noali', '--max', 'mitochondrial.hmm', copy])
        	infernal_call_tRNA = subprocess.Popen(['java', '-jar', 'mitfi/mitfi.jar', '-cm', 'mitfi/t_rna.cm', '-top', copy], stdout=subprocess.PIPE)
		infernal_out_tRNA = infernal_call_tRNA.communicate()[0]
		infernal_call_rRNA = subprocess.Popen(['java', '-jar', 'mitfi/mitfi.jar', '-cm', 'mitfi/r_rna.cm', '-top', copy], stdout=subprocess.PIPE)
		infernal_out_rRNA = infernal_call_rRNA.communicate()[0]
		with open('tRNA_output_sanity_check', 'w') as out:
			out.write(infernal_out_tRNA)
		with open('rRNA_output_sanity_check', 'w') as out:
			out.write(infernal_out_rRNA)
        	shutil.move(copy, better_path)
        	shutil.move(paired_names[0], better_path)
       		shutil.move(paired_names[1], better_path)
        	processed_names=sorted(glob.glob("*.processed"), key=os.path.abspath)
        	for i in range(0,len(processed_names),2):
                	newpair=processed_names[i:i+2]
                	shutil.move(newpair[0],better_path)
                	shutil.move(newpair[1],better_path)
        	f=open(outname[0]+'.gff','w+r')
        	f.write('##gff-version 3\n')
        	GFF_out = [line.split() for line in open(dfamtbl, 'r') if not line.startswith('#')]
        	tRNA_out = [line.split() for line in infernal_out_tRNA if not line.startswith('#')]
		rRNA_out = [line.split() for line in infernal_out_rRNA if not line.startswith('#')]
		for tRNA in tRNA_out:
			if int(float(tRNA[4])) <= 0.001:
				f.write(tRNA[0] + '\t' + 'MiTFi' + '\t' + 'CDS' + '\t' + tRNA[1] + '\t' + tRNA[2] + '\t' + tRNA[3] + '\t' + tRNA[8] + \
				'\t' + '.' + '\t' + 'trn' + tRNA[6] + '(' + tRNA[5].replace('U','T') + ')' + '\n')
		for rRNA in rRNA_out:
                	f.write(rRNA[0] + '\t' + 'MITOS' + '\t' + 'CDS' + '\t' + rRNA[1] + '\t' + rRNA[2] + '\t' + rRNA[3] + '\t' + rRNA[8] + \
			'\t' + '.' + '\t' + 'rrn' + rRNA[6] + '(' + rRNA[5].replace('U','T') + ')' + '\n')
        		for elem in GFF_out:
                		if abs(int(elem[10])-int(elem[9])) >= (0.90 * int(elem[13])):
                        		f.write(elem[2] + '\t' + 'GenBank' + '\t' + 'CDS' + '\t' + elem[9] + '\t' + elem[10] + '\t' + elem[5] + '\t' + elem[8] + '\t.' + '\t' + elem[0] + '\n')
                        		with open(outname[0]+'_genes.fasta', 'a') as fastaout:
                                		with open(name[0]+'-scaffolds.fa.cap.contigs') as fastain:
                                        		lines=fastain.readlines()[2:]
                                        		index = ''.join([item.replace('\n','') for item in lines])
                                        		index1=int(elem[9])-1
                                        		index2=int(elem[10])-1
                                        		if elem[8] == '-':
                                                		revcom = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A'}[B] for B in x][::-1])
                                                		fastaout.write(">"+elem[0]+'\n' + revcom(index[index2:index1]) + '\n')
                                       			else:
								fastaout.write(">"+elem[0]+'\n' + index[index1:index2] + '\n')
if __name__ == "__main__":
	annotation()	
def parallel(annotation):
	def join(annotation, fastq_files):
		from multiprocessing import Pool
		pool = Pool(processes=int(num_cores))
		pool.close()
		pool.join()
	from functools import partial
	return partial(join, annotation)
if __name__ == "__main__":
	parallel()

run.once = parallel(annotation)
