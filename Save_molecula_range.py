import pandas as pd
import numpy as np
import time
import sys
import getopt
from Data_Processing import *
def printUsage():
	print ('''usage: main.py -i <input> -o <output>
       main.py --in=<input> --out=<output>''')
       

def input_info():
	inputarg=""
	outputarg=""
	try:
		opts, args = getopt.getopt(sys.argv[1:],"hi:o:",["in=","out="])
	except getopt.GetoptError:
		printUsage()
		sys.exit(-1)
	for opt,arg in opts:
		if opt == '-h':
			printUsage()
			sys.exit()
		elif opt in ("-i", "--in"):
			inputarg = arg
		elif opt in ("-o","--out"):
			outputarg = arg
	print ('输入:' + inputarg)
	print ('输出:' + outputarg)
	return inputarg, outputarg


 



if __name__ == '__main__':
    inputarg, outputarg = input_info()
    input_param = Input_Param_save_molecula_range(inputarg)
    input_data = Input_data(input_param)
    element_dict, styl_dict = Init_data()
    element_dict_one_list = Counter(element_dict, input_param, input_data)
    atoms_weight_list =  Atoms_weight(element_dict_one_list, input_data, input_param)
    data_all = Read_file_chosed(input_param)
    data_all_dataframe = Read_file_all_frame_to_dataframe(data_all)
    Save_molecules_in_range_all_frame(input_param, data_all_dataframe, element_dict_one_list, outputarg, atoms_weight_list)

