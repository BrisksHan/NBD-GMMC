import networkx as nx
import copy
import numpy as np
import scipy.stats
import multiprocessing
import time


class create_correlation_dict:
    def __init__(self, names, expressions, power = 7, worker = 3): # recommend: parallel_walks = number_walks; if memory error, plz reduce it
        self.names = names
        self.expressions = expressions
        self.power = power
        self.worker = worker
        self.all_walks = []
        self.length = len(names)

    def parallel_calculate_all_correlation(self):
        #t1 = time.time()
        print('start parallel computationing')
        pool = multiprocessing.Pool(processes=self.worker)
        range_number = len(self.names)-1
        all_correlations = pool.map(self.calculate_correlation, range(range_number))
        pool.close()  # Waiting for all subprocesses done..
        pool.join()
        return all_correlations

    def calculate_correlation(self, gene_index):
        correlation_info = []
        for i in range(gene_index+1, self.length):
            coefficient = scipy.stats.pearsonr(np.asarray(self.expressions[gene_index]), np.asarray(self.expressions[i]))[0]
            correlation_info.append(coefficient)
        return correlation_info


def create_soft_threshold_graph_abs_mp(names, expressions, st = 7, worker = 3, ht = 0.01):#this is the one used
    t1 = time.time()
    Graph = nx.Graph()
    print('start calculating correlation')
    ccd = create_correlation_dict(names, expressions, st, worker)
    all_corr = ccd.parallel_calculate_all_correlation()

    for i in range(len(names)):
        Graph.add_node(names[i])

    added_edge_counter = 0
    discarded_edge_counter = 0
    for i in range(len(all_corr)):
        for j in range(len(all_corr[i])):
            #print(np.asarray(expression_values[i]))
            #print(np.asarray(expression_values[j]))
            efficient = all_corr[i][j]
            new_efficient = abs(efficient)
            power_efficent = new_efficient**st
            if power_efficent >= ht:
                Graph.add_edge(names[i], names[i+j+1], weight = power_efficent)
                #Graph.add_edge(names[i+j+1], names[i], weight = power_efficent)
                added_edge_counter += 1
            else:
                discarded_edge_counter +=1
    print('added edges:', added_edge_counter)
    print('discarded edges:',discarded_edge_counter)
    #return Graph
    t2 = time.time()
    print(f'Time for all all_correlations: {(t2-t1):.2f}')  # use multiple cores, total time < sum(time@itr)
    return Graph



def create_soft_threshold_graph_mapped_mp(names, expressions, power = 7, worker = 3, minimum_threshold = 0.01):#this is used as well
    t1 = time.time()
    Graph = nx.Graph()
    print('start calculating correlation')
    ccd = create_correlation_dict(names, expressions, power, worker)
    all_corr = ccd.parallel_calculate_all_correlation()

    for i in range(len(names)):
        Graph.add_node(names[i])

    added_edge_counter = 0
    discarded_edge_counter = 0
    for i in range(len(all_corr)):
        for j in range(len(all_corr[i])):
            #print(np.asarray(expression_values[i]))
            #print(np.asarray(expression_values[j]))
            efficient = all_corr[i][j]
            new_efficient = abs(efficient+1)/2
            power_efficent = new_efficient**power
            if power_efficent >= minimum_threshold:
                Graph.add_edge(names[i], names[i+j+1], weight = power_efficent)
                #Graph.add_edge(names[i+j+1], names[i], weight = power_efficent)
                added_edge_counter += 1
            else:
                discarded_edge_counter += 1
    print('added edges:', added_edge_counter)
    print('discarded edges:',discarded_edge_counter)
    #return Graph
    t2 = time.time()
    print(f'Time for all all_correlations: {(t2-t1):.2f}')  # use multiple cores, total time < sum(time@itr)
    return Graph




def select_top_variance_gene(gene_names, expressions, top_variance_gene_number):
    alls = []
    for i in range(len(expressions)):
        std = np.std(expressions[i])
        alls.append(std)

    expression_orders = np.argsort(np.array(alls))[-top_variance_gene_number:]
    new_gene_names = []
    new_expressions = []
    for i in range(len(expression_orders)):
        new_gene_names.append(copy.deepcopy(gene_names[expression_orders[i]]))
        new_expressions.append(copy.deepcopy(expressions[expression_orders[i]]))

    return copy.deepcopy(new_gene_names), copy.deepcopy(new_expressions)





