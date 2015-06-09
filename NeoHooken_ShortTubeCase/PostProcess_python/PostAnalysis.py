# post-process the data from the dumped file: out_###.gmv or prefix_###.gmv
# READ: 	the keyline(s), keyNode(s), and store them in matrix
# WRITE:	output matrix to file; a typical form: lines:
# VISUALize: 	output some figures

# we record: var_info: all variables name
# 	     var_array_plot: to plot the variables
# ex. keyLine =''
import numpy as np
import string
from CPostVariable import PostVariable
import matplotlib.pyplot as plt

def postAnalysis(file_input):
	# deal with postFile_input
	line_no = 0;
	with open(file_input,'r') as f:	
		# 1st line:	
		read_line = f.readline();
		line_no += 1;
		print " 	line #", line_no, ":", read_line;
		
		str_array = string.split(read_line);
		file_pre = str_array[0];
		file_id_bgn = int(str_array[1]);
		file_id_end = int(str_array[2]);
		file_num = file_id_end + 1 -file_id_bgn;	
		# 2nd line: number of variables
		read_line = f.readline();
		line_no += 1;
		print " 	line #", line_no, ":", read_line;		
		str_array = string.split(read_line);
		var_num = int(str_array[0]);		
		
		# read all the lines:
		is_set_common = 0; # not set the common info;
		###########
		# 1st) st varibale: set the common info
		k_var =0; # var_id
		read_line = f.readline();
		line_no += 1;
		print " 	line #", line_no, ":", read_line;
		## deal with the variable		
		str_array = string.split(read_line);
		var_name = str_array[0];

		id_array = getIdArray(str_array);		
		
		print " 	var #", k_var, ":", id_array;	
		## set the class:
		var_info = [var_name];	
		var_array_plot=np.zeros((file_num,var_num));	
		var_a = PostVariable(var_name,id_array);
		var_a.setCommonInfo(file_pre, file_id_bgn, file_id_end);
		var_a.readVarInfo();
		var_a.writeVarInfo();
		var_plot = var_a.ObtainArray(0);
		var_array_plot[:,0] = var_plot; 
		###########
		k_var +=1;
		while (k_var < var_num):
			## read the file
			read_line = f.readline();
			line_no += 1;	
			print " 	line #", line_no, " var #", k_var,":", read_line;

			## deal with the variable
			str_array = string.split(read_line);
			var_name = str_array[0];
			id_array = getIdArray(str_array);
			## set the class:
			var_info.append(var_name);
			var_a = PostVariable(var_name,id_array);
			var_a.setCommonInfo(file_pre, file_id_bgn, file_id_end);
			var_a.readVarInfo();
			var_a.writeVarInfo();
			var_plot = var_a.ObtainArray(0);
			var_array_plot[:,k_var] = var_plot; 
			## move to the next
			k_var +=1;
					
		##########
		print "#####++++ finish the postProcess"
		dealWithPlot(var_info,var_array_plot);
## analysis for the id_info from a string array
# fill up the id_array:: find proper integer before #;	
def getIdArray(str_array):	
	k_pos = 1;
	id_array = [int(str_array[k_pos])];
	total_pos = len(str_array);
	Is_end =0;
	while (Is_end<1 and k_pos < total_pos-1):
		k_pos +=1;
		if (str_array[k_pos][0] == '#'): # end of the id
			Is_end =1;
		else: # not end: fill in the id_array
			id_array.append(int(str_array[k_pos]));
	print "******** check id_array:", id_array	
	return id_array

## we plot the data
def dealWithPlot(var_info, var_array_plot):
	plotname ="pressure-dilation.png";
	fig=plt.figure()
	ax=fig.add_axes([0.1,0.1,0.6,0.75]);
	str_plot = ['r','b','g','m'];
	
	vec_x = np.arange(len(var_array_plot[:,0]));
	k_var =0;
	print "len(var_info):", len(var_info), ",len var_array[0,:]",len(var_array_plot[0,:]),",len var_array[:,0]",len(var_array_plot[:,0])
	while (k_var < 4 and k_var < len(var_info)):
		print "k_var:", k_var
		ax.plot(vec_x,var_array_plot[:,k_var],str_plot[k_var]+"-s",lw=2,label=var_info[k_var]);
		k_var +=1;

	ax.set_ylabel('var', size=12)
	ax.set_xlabel('pressure/step', size=12)
	ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	ax.set_title("pressure-dilation",size=18)
	plt.grid(True)	
	plt.savefig(plotname)
	plt.show();



