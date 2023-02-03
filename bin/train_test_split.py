train_ags = []
test_ags = []

for cell_type in cell_list:
    subset = id_taxonomy[id_taxonomy['taxonomy_name'] == cell_type]
    num_samples = len(subset)
    count = .8*num_samples
    
    train = subset.sample(round(count), random_state=1010)
    train_ags.append(train['ID'])
    test = subset[~subset['ID'].isin(train['ID'])]
    test_ags.append(test['ID'])