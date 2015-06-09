# post-process the data from the dumped file: out_###.gmv or prefix_###.gmv
# READ: 	the keyline(s), keyNode(s), and store them in matrix
# WRITE:	output matrix to file; a typical form: lines:
# VISUALize: 	output some figures

# ex. keyLine =''

import PostAnalysis
import sys



if __name__ == "__main__":
    # command-line argument
    print '+++++++ running: ', sys.argv
    my_argv=sys.argv[1:]
    alen=len(my_argv)
    file_input="PostProcessInput"
    if (alen<1):
	print '+++++++ please given file info:'
	print '+++++++ otherwise, using default file', file_input
    else:
	file_input = str(my_argv[0])
    # begin to process the input_file
    PostAnalysis.postAnalysis(file_input);
    print "+++++++ finish post-analysis of the file", file_input
    

