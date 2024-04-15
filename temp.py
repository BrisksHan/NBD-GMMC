import utils
import numpy as np
import read_file

#all_data = utils.load_any_obj_pkl('input/BRCA_all_data_0.pkl')

#gene_names = all_data[0][0]
#sample_names = all_data[0][1]
#training_expressions = all_data[0][2]
#training_sample_name = all_data[0][3]
#test_expressions = all_data[0][4]
#test_sample_name = all_data[0][5]

#print(np.shape(training_expressions))
#print(len(training_expressions[0]))

#print(len(training_sample_name))

#utils.save_any_obj_pkl(training_expressions, 'input/expression.pkl')
#utils.save_any_obj_pkl(training_sample_name, 'input/names.pkl')

read_file.read_file()

