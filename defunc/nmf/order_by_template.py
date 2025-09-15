from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
import pandas as pd
import numpy as np
import matplotlib.colors as mcolors
import os


# W_old_path = '/net/seq/data2/projects/sabramov/SuperIndex/dnase-index0415/embeddings/NMF/nndsvda_29c_050124//nmf_results/nndsvda.top35_max100_no_cancer_seed2674.peak_mask@cf_0.05@mv_0.7.29.W.npy'
W_old_path = '/net/seq/data2/projects/sabramov/SuperIndex/hotspot3/w_babachi_new.v17/embeddings/NMF/full_matrix/output/nmf/hotspot3_index.ontology_extended_spot_20_weights_norm.29/hotspot3_index.ontology_extended_spot_20_weights_norm.29.W.npy'


def reorder_components(W1, W2):
    print(W1.shape, W2.shape)
    cosine_similarity = 1 - cdist(W1, W2, 'cosine')
    _, col_ind = linear_sum_assignment(-cosine_similarity)
    return col_ind


def get_component_data(W, W_old=None):
    component_data = pd.DataFrame(
        [
        { 'color': '#BB2DD4',
            'index': 18,
            'name': 'T-cell',
            'short_name': 'TCL'},
        { 'color': '#E6009B',
            'index': 12,
            'name': 'B-lymphoblast',
            'short_name': 'BLY'},
        { 'color': '#FF7600',
            'index': 4,
            'name': 'Macrophage/Monocyte',
            'short_name': 'MAC'},
        { 'color': '#ED2024',
            'index': 14,
            'name': 'Erythroid',
            'short_name': 'ERY'},
        { 'color': '#AE1618',
            'index': 25,
            'name': 'CMP',
            'short_name': 'CMP'},
        { 'color': '#32A63A',
            'index': 8,
            'name': 'Adult tissues',
            'short_name': 'ADT'},
        { 'color': '#00AD95',
            'index': 28,
            'name': 'Adrenal',
            'short_name': 'ADR'},
        { 'color': '#3BE796',
            'index': 27,
            'name': 'Thyroid',
            'short_name': 'THY'},
        { 'color': '#204784',
            'index': 2,
            'name': 'Basal ganglia',
            'short_name': 'BGL'},
        { 'color': '#0479FD',
            'index': 24,
            'name': 'Caudate',
            'short_name': 'CAU'},
        { 'color': '#37C1F1',
            'index': 7,
            'name': 'Cerebellar',
            'short_name': 'CER'},
        { 'color': '#9ADBF9',
            'index': 23,
            'name': 'Brain (fetal)',
            'short_name': 'FBR'},
        { 'color': '#7278DA',
            'index': 21,
            'name': 'Eye',
            'short_name': 'VIS'},
        { 'color': '#D0D2F8',
            'index': 17,
            'name': 'Fetal tissues',
            'short_name': 'FTL'},
        { 'color': '#C2FFED',
            'index': 11,
            'name': 'Lung (fetal)',
            'short_name': 'LNG'},
        { 'color': '#C2FF99',
            'index': 5,
            'name': 'Muscle (fetal)',
            'short_name': 'MSC'},
        { 'color': '#8CC63F',
            'index': 16,
            'name': 'Heart (fetal)',
            'short_name': 'HRT'},
        { 'color': '#4A6876',
            'index': 15,
            'name': 'Intestine (fetal)',
            'short_name': 'INT'},
        { 'color': '#FFDAB9',
            'index': 1,
            'name': 'Kidney (fetal)',
            'short_name': 'KDN'},
        { 'color': '#FFC26A',
            'index': 26,
            'name': 'Epithelial (kidney)',
            'short_name': 'EPK'},
        { 'color': '#FFDE00',
            'index': 13,
            'name': 'Epithelial',
            'short_name': 'EPI'},
        { 'color': '#ECED9D',
            'index': 6,
            'name': 'Myoblast',
            'short_name': 'MYB'},
        { 'color': '#759CD5',
            'index': 22,
            'name': 'Endothelial',
            'short_name': 'END'},
        { 'color': '#6B54A3',
            'index': 3,
            'name': 'Stromal A',
            'short_name': 'STR'},
        { 'color': '#45316B',
            'index': 20,
            'name': 'Stromal B',
            'short_name': 'STR2'},
        { 'color': '#D48EB8',
            'index': 19,
            'name': 'Placental',
            'short_name': 'PLC'},
        { 'color': '#FFABF1',
            'index': 10,
            'name': 'Embryonic A',
            'short_name': 'EMB'},
        { 'color': '#FFD7F8',
            'index': 9,
            'name': 'Embryonic B',
            'short_name': 'EMB2'},
        { 'color': '#E7E9E9',
            'index': 0,
            'name': 'TSS',
            'short_name': 'TSS'}
        ]
    )

    component_colors = component_data['color'].values

    if W_old is None:
        if os.path.exists(W_old_path):
            W_old = np.load(W_old_path).T
            if W_old.shape[1] != W.shape[1] or W.shape[0] <= 29:
                W_old = None

        if W_old is None:
            if W.shape[0] <= 29:
                colors = component_colors[:W.shape[0]]
            else:
                colors = define_colors(W.shape[0])
            return pd.DataFrame({
                'index': np.arange(W.shape[0]),
                'color': colors,
                'name': ['Component {}'.format(i) for i in range(W.shape[0])],
            })
    reorder = reorder_components(W, W_old)

    # #FIXME
    component_data = component_data.reset_index(
        names='color_order'
    ).set_index('index').loc[
        reorder
    ].reset_index(
        drop=True
    ).reset_index(
        names='index'
    ).sort_values(
        'color_order'
    ).reset_index(
        drop=True
    ).drop(
        columns='color_order'
    )

    return component_data


def define_colors(n_components):
    # Function to define colors for components. Borrowed from Meuleman et al. 2020
    # Based on 2020 component colors
    comp_colors = ["#ffe500",
                    "#fe8102",
                    "#ff0000",
                    "#07af00",
                    "#4c7d14",
                    "#414613",
                    "#05c1d9",
                    "#0467fd",
                    "#009588",
                    "#bb2dd4",
                    "#7a00ff",
                    "#4a6876",
                    "#08245b",
                    "#b9461d",
                    "#692108",
                    "#c3c3c3",
                    "#6630a6",
                    "#ffc26a",
                    "#fc197e",
                    "#759cd5",
                    "#a6da57",
                    "#343331",
                    "#d04299",
                    "#a1efff",
                    "#ffadf1",
                    "#fef3bb",
                    "#61567b",
                    "#ffc9c6",
                    "#d0d2f8"]
    neworder = np.array([16,10,7,11,2,12,1,8,4,15,14,5,9,6,3,13]).astype(int) - 1
    
    component_colors = list(np.array(comp_colors)[neworder]) + comp_colors[16:]

    maxassigned = len(component_colors)
    
    if (n_components > maxassigned):
        # somewhat defunct but whatever. Adds extra "random" colors if you use more than 16 
        colornames = np.sort(list(mcolors.CSS4_COLORS.values()))
        colornames = list(set(colornames) - set(comp_colors))
        count = maxassigned
        np.random.seed(100)
        while (count < n_components):
            new_color = colornames[np.random.randint(len(colornames))]
            trialcount = 0
            while ((new_color in component_colors) and (trialcount < 100)):
                trialcount += 1
            #print('new color', count, new_color)
            component_colors.append(new_color)
            count += 1
    return component_colors[:n_components]