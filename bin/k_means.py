import sys
import numpy as np
import pandas as pd
import json
from sklearn.preprocessing import LabelEncoder
import pickle
from sklearn.cluster import KMeans
from aggloclustering import get_clustering_metrics


algoirithm = ["lloyd", "elkan"]


def main(json_object, meta, embedding, training_set, validation_set):
   
    #Subset embedding and metadata by training and validation_set
    #Remember emdbedding currenlty has ['ID'] as the last column
    train = embedding[embedding['ID'].isin(training_set['ID'])]
    valid = embedding[embedding['ID'].isin(validation_set['ID'])]
    del train[train.columns[-1]]
    del valid[valid.columns[-1]]
    del embedding[embedding.columns[-1]]

    #Subset MetaData
    train_meta = meta[meta['ID'].isin(training_set['ID'])]
    valid_meta = meta[meta['ID'].isin(validation_set['ID'])]

    #Run Clustering Algorithm
    metrics_rows = []
    label_encoder = LabelEncoder()
    models = []

    true_col = 'ontology_term'
    true_labels = label_encoder.fit_transform(meta[true_col])
    true_train_labels = label_encoder.fit_transform(train_meta[true_col])
    #Do I just use the true_label encodings or the encodings from validation set separately?? Same for training
    true_valid_labels = label_encoder.fit_transform(valid_meta[true_col]) 

    for index, algorithm in enumerate(algorithm):
        params = {**json_object, 'algorithm': algorithm}
        try:
	        #Full Embedding
            clustering_model = KMeans(**params)
            clustering_model.fit(embedding)

	        #Train and Validation
            test_model = KMeans(**params) 
            valid_labels = test_model.predict(valid)
        except:
            continue
        clustered_labels =  clustering_model.labels_
        training_labels = test_model.labels_
       
        params_str = json.dumps(params)
        
	#Full Set
    ami, ari, fm, _, _, jaccard, _ = get_clustering_metrics(clustered_labels, true_labels)
    metrics_rows.append([index, params_str, ami, ari, fm, jaccard])
    models.append(
        (index, pd.DataFrame(clustered_labels), clustering_model)
    )

    # Training Set
    ami, ari, fm, _, _, jaccard, _ = get_clustering_metrics(training_labels, true_train_labels) 
    metrics_rows.append([index, params_str, ami, ari, fm, jaccard])
    models.append(
            (index, pd.DataFrame(training_labels), test_model)
        )

    # Validation Set
    ami, ari, fm, _, _, jaccard, _ = get_clustering_metrics(valid_labels, true_valid_labels)
    metrics_rows.append([index, params_str, ami, ari, fm, jaccard]) 
    models.append(
	    (index, pd.DataFrame(valid_labels), test_model)  # Will be the same model as training
	)

    metrics_df = pd.DataFrame.from_records(metrics_rows,
        columns=['id', 'params', 'ami', 'ari', 'fm', 'jaccard'])
    return metrics_df, models


if __name__ == '__main__':
    with open(sys.argv[1], 'r') as openfile:
        params = json.load(openfile)
    embedding = np.load(sys.argv[2])
    samples = pd.read_table("/net/seq/data2/projects/sabramov/SuperIndex/dnase-0108/output/index/indivs_order.txt", header=None)
    samples = samples.T
    samples.columns = ['ID']
    embedding['ID'] = samples['ID']
    
    #Load Metadata
    meta_columns = ['tc_number', 'ds_number', 'ln_number', 'id', 
                'passed_agg_qc', 'hotspot1_spot', 'hotspot2_spot', 'nuclear_percent_duplication', 'taxonomy_name', 
                'ontology_id', 'ontology_term', 'system', 'subsystem', 'organ', 'organ_region', 'side_position', 
                'cell_type', 'cell_markers', 'disease_pathology', 'sampling_origin', 'state', 'germ_layer', 
                'growth_stage', 'description', 'donor_id', 'donor_age', 'donor_age_units', 'donor_sex', 
                'donor_ethnicity', 'donor_health_condition', 'vendor_name', 'perturbation', 'treatment_name', 
                'treatment_details', 'dose', 'dose_units', 'filtered_alignments_bam', 'bam_index', 'hotspot_peaks', 
                'hotspot_peaks_1per', 'hotspot_peaks_point1per', 'adapter', 'total_alignments', 'all_aligned', 
                'all_mapq_filter', 'autosomal_align', 'chr_m', 'hotspot1_num_bases', 'hotspot1_num_spots', 
                'hotspot2_num_bases', 'hotspot2_num_spots', 'insert_ft_eleven', 'insert_ls_ratio', 'paired_aligned', 
                'paired_autosomal_align', 'paired_nuclear_align', 'picard_median_insert_size', 'picard_percent_duplication',
                'raj_dataset', 'tag', 'footprintdir', 'num_footprints']
    
    meta = pd.read_table(sys.argv[3], header=None, names=meta_columns)
    prefix=sys.argv[4]
    training_set = pd.read_table(sys.argv[5], header=None, names=['ID'])
    validation_set = pd.read_table(sys.argv[6], header=None, names=['ID'])

    metrics_df, models = main(params, meta, embedding, training_set, validation_set)
    for index, labels, clustering in models:
        labels.to_csv(f"{prefix}.{index}.annotations.txt", header=None, index=False)
        labels.to_csv(f"{prefix}.{index}.annotations.txt", header=None, index=False)
        with open(f"{prefix}.{index}.model.pkl", 'wb') as out:
            pickle.dump(clustering, out)

    metrics_df.to_csv(f'{prefix}.metrics.tsv', sep="\t", header=None, index=False)
