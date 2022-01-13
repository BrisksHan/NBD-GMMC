import utils
import SN_construction
import Network_embedding
#import visualization
import numpy as np
import csv
import copy
import Module_detection
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
import Eigengene_significance
import read_file

#normal_worker = 10, walks_worker = 1, project_file_name = 'THCA/THCA_all_data_'
def parse_args():
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter, conflict_handler='resolve')
    # -----------------------------------------------general settings--------------------------------------------------
    parser.add_argument('--input_folder', default='input/')
    parser.add_argument('--output_folder', default='output/')
    parser.add_argument('--nworker', default=4,
                        help='normalworkernumber')
    parser.add_argument('--wworker', default=4,
                        help='walksworkernumber')
    parser.add_argument('--projectname', default='TCGA_LIHC.txt',
                        help='walksworkernumber')
    parser.add_argument('--window', default=5,
                        help='window size of deepwalk')
    parser.add_argument('--walknumber', default=10,
                        help='the number of indenpedent random walk from each node')
    parser.add_argument('--walklength', default=80,
                        help='the size of length for each random walk')
    parser.add_argument('--ht', default=0.01,
                        help='the parameter of hard threshold')
    parser.add_argument('--st', default=7,
                        help='the parameter of soft threshold')
    parser.add_argument('--dim', default=5,
                        help='the dimension of node embedding')
    parser.add_argument('--save_embedding', default = False,
                        help='save the node embedding')
    parser.add_argument('--translate_implied_name', default= True,
                        help='translate the implied name from TCGA or the sample label must be')
    parser.add_argument('--top_variance_number', default = 2000,
                        help = 'select the top n genes that have the highest variance between sample expressions')

    args = parser.parse_args()
    return args

def data_with_validation(projectname, foldername, top_variance_number = 1000, ht = 0.01, st = 7,  workers = 3, emb_dim = 5, lg_dim = 4, walks_worker = 3, window = 5, walklength = 80, walknumber = 10, save_embedding = False):

    all_gene_names, all_gene_expression, sample_names = read_file.read_file(foldername+projectname)
    gene_names, gene_expression = SN_construction.select_top_variance_gene(all_gene_names, all_gene_expression, top_variance_number)
    graph = SN_construction.create_soft_threshold_graph_abs_mp(gene_names, gene_expression, st = st, worker = workers, ht = ht)
    print('start conducting random walks')
    
    DW = Network_embedding.DeepWalk([graph], workers = workers, emb_dim = emb_dim, num_walks=walknumber, walk_length=walklength, window=window, walks_workers= walks_worker)#577
    print('start skip-gram')
    output = DW.sampling_traning()[0]
    print("gene embedding finished")

    trained_embedding_name = 'output/'+projectname+'_top'+str(top_variance_number)+'_'+"emb"+str(emb_dim)+"and"+str(st)+'_'+str(ht)+'_'+str(window)+".pkl"
    if save_embedding == True:
        utils.save_any_obj_pkl(output, trained_embedding_name)
    
    embedded_gene_names, gene_cluster, gmm_probility_list = Module_detection.BGM_soft(output)

    gene_cluster_info_name = 'output/'+projectname+'_'+"gene_cluster_es"+'_top'+str(top_variance_number)+'_'+str(st)+'_'+str(ht)+'_'+str(window)+'_'+str(emb_dim)+".csv"
    #write_gene_names_clusters_csv(embedded_gene_names, gene_cluster, gene_cluster_info_name)

    np_cluster = np.array(gene_cluster)
    np_cluster_media = []
    for i in range(len(np_cluster)):#labels -1 would be ignored
        if np_cluster[i] >= 0:
            np_cluster_media.append(np_cluster[i])
    np_cluster = np.array(np_cluster_media)
    count = 0
    all_gene_info_string = ""
    np_unique = np.unique(np_cluster)
    print("unique cluster:",np_unique)

    module_es = {}

    for g in np.unique(np_cluster):

        print(count," gene cluster out of ",len(np_unique))
        count+=1
        cluster_indexes = [i for i, x in enumerate(np_cluster) if x == g]
        #cluster_indexes = list(cluster_indexes)
        cluster_gene_names = []
        for i in range(len(cluster_indexes)):
            cluster_gene_names.append(gene_names[int(cluster_indexes[i])])

        new_cluster_expressions = get_sample_expressions(gene_names, cluster_gene_names, gene_expression, sample_names)
        
        labels = utils.get_sample_implied_labels(sample_names)
        print("sample number:",len(labels))
        transposed_result = np.transpose(new_cluster_expressions)
        eigengene = Eigengene_significance.PCA_decomposition(transposed_result)#astype(float)
        #module_es.append()
        #print('train sample gene expressions length:',len(new_cluster_expressions[0]))
        cur_module_es = Eigengene_significance.Eigengene_significance(eigengene, labels)
        print("current module es:", cur_module_es)
        module_es[g] = cur_module_es

    write_gene_names_clusters_es_csv(gene_names, gene_cluster, module_es, gene_cluster_info_name)

        #print('test sample gene expressions length:',len(new_cluster_expressions[0]))
    #return [highest_train, highest_test]#[highest_train_Gmeans, HT_test_Gmeans, highest_train_mcc, HT_test_mcc, highest_train_BM, HT_test_BM]

def write_gene_names_clusters_es_csv(gene_names, gene_cluster, module_es, name):
    all_infos = []
    for i in range(len(gene_names)):
        cur_cluster = module_es[gene_cluster[i]]
        single_item = [gene_names[i], gene_cluster[i], cur_cluster]
        all_infos.append(single_item)
    with open(name, "w",newline="") as f:
        writer = csv.writer(f)
        writer.writerows(all_infos)

def get_sample_expressions(gene_names, cluster_gene_names, training_expressions, training_sample_name):#training_expressions
    print("selected_gene_numbers:",len(cluster_gene_names))
    #sample_names, sample_expressions = get_sample_names_and_expressions(with_healthy_samples = with_healthy_samples)
    #gene_names, gene_expressions = get_gene_names_and_expressions(with_heathy_samples= with_healthy_samples)
    selected_sample_expressions = [[] for i in range(len(training_sample_name))]
    for i in range(len(gene_names)):
        if (gene_names[i] in cluster_gene_names) == True:
            #print("to be completed")
            for j in range(len(selected_sample_expressions)):
                selected_sample_expressions[j].append(training_expressions[i][j])
    return copy.deepcopy(selected_sample_expressions)


def write_final_csv(lines, name):
    print("write all result")
    with open(name, "w",newline="") as f:
        writer = csv.writer(f)
        writer.writerows(lines)

def NBD_GMMC(args):#KIRC/KIRC_all_data_
    
    normal_worker =  int(args.nworker)
    walks_worker = int(args.wworker)

    top_variance_number = int(args.top_variance_number)

    project_file_name = args.projectname

    save_embedding = bool(args.save_embedding)

    input_folder = args.input_folder
    output_folder = args.output_folder

    window = int(args.window)
    walknumber = int(args.walknumber)
    walklength = int(args.walklength)

    ht = float(args.ht)
    st = float(args.st)

    emb_dim = int(args.dim)

    print('project_name:', project_file_name)
    all_evaluation = []


    current_evaluation = data_with_validation(project_file_name, input_folder, top_variance_number = top_variance_number, ht = ht, st = st, workers= normal_worker, emb_dim = emb_dim, walks_worker = walks_worker, window = window, walklength = walklength, walknumber = walknumber, save_embedding = save_embedding)
    #all_evaluation.append(current_evaluation)

    #output_name = output_folder + project_file_name+'_'+str(window)+'_'+str(emb_dim)+'_'+str(st)+'_'+str(ht)+'.csv'
    #write_final_csv(all_evaluation, output_name)


if __name__ == "__main__":
    NBD_GMMC(parse_args())
    