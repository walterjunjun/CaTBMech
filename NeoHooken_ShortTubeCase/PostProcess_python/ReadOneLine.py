#only read one line, used in IBFE 
# ex. keyLine =''
import numpy as np
import string
def ReadOneLine(keyLine, fread):
	PreStr=keyLine;
	IsStart = -1;
	
	with open(fread,'r') as f:		
		while (IsStart < 0): # not find the starting line
			 read_line = f.readline();
			 IsStart = read_line.find(PreStr);
		# then we find, begin to get the 1st line
		read_line = f.readline();
		str_array = string.split(read_line);
		my_array = np.array(map(float, str_array));		

	return my_array
