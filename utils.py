import pickle
import platform
import csv

print(platform.architecture())
#-----------------------------------------IO-----------------------------

import platform
platform.architecture()


def load_any_obj_pkl(path):
    ''' load any object from pickle file
    '''
    with open(path, 'rb') as f:
        any_obj = pickle.load(f)
    return any_obj

def save_any_obj_pkl(obj, path):
    ''' save any object to pickle file
    '''
    with open(path, 'wb') as f:
        pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)

def write_gene_names_clusters_csv(gene_names, gene_cluster, name):
    all_infos = []
    for i in range(len(gene_names)):
        single_item = [gene_names[i], gene_cluster[i]]
        all_infos.append(single_item)
    with open(name, "w",newline="") as f:
        writer = csv.writer(f)
        writer.writerows(all_infos)

def write_final_csv(lines, name):
    with open(name, "w",newline="") as f:
        writer = csv.writer(f)
        writer.writerows(lines)


def get_sample_implied_labels(IDs):
    matched_stages = []
    for i in range(len(IDs)):
        sample_type = int(IDs[i][13:15])
        if sample_type >= 10:
            matched_stages.append(1)
        else:
            matched_stages.append(0)
    return matched_stages


