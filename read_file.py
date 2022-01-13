import copy

def read_gene_file(file_name = 'input/TCGA_LIHC.txt'): 
    f = open(file_name, "r")
    lines = f.readlines()
    all_data = []    
    for i in range(len(lines)):
        single_gene = lines[i].split()
        all_data.append(single_gene)
    f.close()
    return all_data


def read_file(file_name = 'input/TCGA_LIHC.txt'):
    f = open(file_name)
    lines = f.readlines()
    raw_data = []    
    for i in range(len(lines)):
        single_gene = lines[i].split()
        raw_data.append(single_gene)
    f.close()
    del raw_data[0][1]
    sample_names = []
    for i in range(1, len(raw_data[0])):
        sample_names.append(raw_data[0][i])
    gene_names = []
    for i in range(1, len(raw_data)):
        gene_names.append(raw_data[i][0])
    gene_expression_list = [[] for i in range(len(raw_data)-1)]
    for i in range(1, len(raw_data)):
        for j in range(1, len(raw_data[i])):
            gene_expression_list[i-1].append(float(raw_data[i][j]))
    return copy.deepcopy(gene_names), copy.deepcopy(gene_expression_list), copy.deepcopy(sample_names)


