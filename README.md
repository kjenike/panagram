# Panagram: Interactive, alignment-free pan-genome browser  

Add description and use for panagram

#Useage:
python plot_pangene.py data.file.txt SAMPLENAME
python plot_pangene.py test_data/result_SQUI2.chr1.txt SQUI2

#To view the plots, go to this link in a web browser: http://127.0.0.1:8050/

#Necessary input file for plotting:
	1) phylogenetic tree file in ".treefile" format 
	2) gff files for gene annotations for each sample
	3) gff files for repeat annotations for each sample

########################
If you don't have the above input files, you can run the pre_processing.sh script.
This will require the following files:
	1) fasta files in single line format for each genome in the pangenome



