#!/usr/bin/python
#Python class: PostProcessedVariable
#"common info": fileInfo, file_num; set_common_info;
#"private info": data: var_name; var_data;
#		 operation: read_var_info; write_var_data
# we store the var_data as: row vector: nrow = gvmfile_num;
# we must use scope to access memeber: (such as className or self)
import numpy as np
import string
class PostVariable:
	#'common base info on file
	file_pre="",
	file_id_bgn=0;
	file_id_end=0;
	set_common_info = 0;

# initialization: 
	def __init__(self, name, id_array):
		self.name = name
		self.id_array = id_array
		self.data_array=np.empty([2,1]);

	def setCommonInfo(self,file_pre, file_bgn, file_end):
		PostVariable.file_pre=file_pre;
		PostVariable.file_id_bgn = file_bgn;
		PostVariable.file_id_end = file_end;
		PostVariable.set_common_info =1;
		print "+++++++++ set once: common info:"
		print "		file_pre: ", PostVariable.file_pre, "; id_begin:", PostVariable.file_id_bgn, "; id_end:", PostVariable.file_id_end
	def readVarInfo(self): # generate the data
		file_num = PostVariable.file_id_end - PostVariable.file_id_bgn +1;
		id_num =len(self.id_array);
		self.data_array=np.zeros((file_num,id_num));
		print"	******** check file_num, id_num:", file_num, id_num
		# read from file
		k_file = PostVariable.file_id_bgn;
		file_count = 0;
		while(k_file<PostVariable.file_id_end +1):
			fread=  PostVariable.file_pre+"_%03d" % k_file + ".gmv";
			keyLine = self.name+ " 1";
			print" 	******* check:read ", keyLine, " in file:", fread;
			my_array = ReadOneLine(keyLine, fread);
			print "	******* check: read data ", my_array[self.id_array];
			self.data_array[file_count,0:id_num] = my_array[self.id_array];	
			#move to next
			file_count +=1;		
			k_file +=1;		


		# print finish reading
		print "+++++++++ finishing reading data"
	def writeVarInfo(self): # dump the data
		file_num = PostVariable.file_id_end - PostVariable.file_id_bgn +1;
		id_num =len(self.id_array);		
		DumpName = self.name+ ".txt";
		#headName = str(file_num) + "," + str(id_num);
		np.savetxt(DumpName, self.data_array); #, header=headName);	
		# print finish outputing
		print "+++++++++ finishing dumping data to file", DumpName

	def ObtainArray(self, col_num): # read the data
		return self.data_array[:,col_num];
	
# read data from the *.gmv file
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
