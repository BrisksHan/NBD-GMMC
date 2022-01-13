import time
import utils
import numpy as np
from sklearn import mixture
import math

def read_embedding_result(embeddings):
    #embeddings = utils.load_any_obj_pkl(path)
    IDs = list(embeddings)
    final_embeddings = []
    for i in range(len(IDs)):
        #print(embeddings[IDs[i]])
        final_embeddings.append(embeddings[IDs[i]])
    return IDs, final_embeddings



def BGM_soft(embedding_result, k_max = None, concentration_prior = 0.1, minimum_community = 4, covariance_type = 'diag'):
    print('start DP soft')
    t1 = time.time()
    IDs, embeddings = read_embedding_result(embedding_result)
    if k_max == None:
        k_max = int(math.sqrt(len(embeddings)))*2
    embeddings = np.array(embeddings)

    bic = []
    models = []
    lowest_aic = 1
    
    model = mixture.BayesianGaussianMixture(n_components=k_max, covariance_type = covariance_type, weight_concentration_prior = concentration_prior/(k_max), max_iter = 1000)

    model.fit(embeddings)

    initial_labels = model.predict(embeddings)
    initial_labels_probas = model.predict_proba(embeddings)

    refined_labels = soft_get_refined_labels(initial_labels, initial_labels_probas, minimum_community)

    t2 = time.time()
    print('time for clustering:',t2-t1)
    return IDs, refined_labels, initial_labels_probas


def soft_get_refined_labels(initial_labels, initial_labels_probas, minimum_community):
    np_initial_labels = np.array([initial_labels], dtype = int)
    initial_comm_IDs, counts = np.unique(np_initial_labels, return_counts=True)

    ignored_communities = []
    for i in range(len(counts)):
        if counts[i] < minimum_community:
            ignored_communities.append(initial_comm_IDs[i])
    
    new_labels = []

    for i in range(len(initial_labels)):
        current_label = initial_labels[i]
        #print(current_label)
        #print(ignored_communities)
        if (current_label in ignored_communities) == True:
            #print(initial_labels_probas[i])
            new_label = get_desired_label(initial_labels_probas[i], ignored_communities)
            new_labels.append(new_label)
        else:
            new_labels.append(current_label)
    return new_labels

def get_desired_label(a_prob, ignored_communities):
    sorted_prob = sorted(a_prob)
    count = -1
    while True:
        #print(list(a_prob).index(sorted_prob[count]))
        new_index = list(a_prob).index(sorted_prob[count])
        if (new_index in ignored_communities) == False:
            return new_index
        count -= 1
