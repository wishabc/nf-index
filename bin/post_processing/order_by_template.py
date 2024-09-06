from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
import pandas as pd
import numpy as np
import matplotlib.colors as mcolors
import os


W_old_path = '/net/seq/data2/projects/sabramov/SuperIndex/dnase-index0415/embeddings/NMF/nndsvda_29c_050124//nmf_results/nndsvda.top35_max100_no_cancer_seed2674.peak_mask@cf_0.05@mv_0.7.29.W.npy'


def reorder_components(W1, W2):
    print(W1.shape, W2.shape)
    cosine_similarity = 1 - cdist(W1, W2, 'cosine')
    _, col_ind = linear_sum_assignment(-cosine_similarity)
    return col_ind


def get_component_data(W, W_old=None):
    component_order = [
        27,
        6,
        2,
        25,
        11,
        14,
        5,
        15,
        28,
        24,
        10,
        18,
        3,
        13,
        26,
        23,
        1,
        7,
        12,
        17,
        4,
        20,
        21,
        19,
        16,
        8,
        9,
        22,
        0,
    ]

    component_colors = [
        '#BB2DD4',
        '#E6009B',
        '#FFDE00',
        '#FF7600',
        '#ED2024',
        '#AE1618',
        '#414613',
        '#32A63A',
        '#3BE796',
        '#00AD95',
        '#4A6876',
        '#08245B',
        '#0479FD',
        '#37C1F1',
        '#A1EFFF',
        '#45316B',
        '#6B54A3',
        '#759CD5',
        '#FEF3BB',
        '#FFDAB9',
        '#FFC26A',
        '#D0D2F8',
        '#C2FFED',
        '#C2FF99',
        '#8CC63F',
        '#D48EB8',
        '#FF9BC7',
        '#FFABF1',
        '#E7E9E9',
    ]

    component_names = [
        'Lymphoid (B-cell)',
        'Lymphoid (T-cell/NK-cell)',
        'Myeloid (Monocyte)',
        'Myeloid (Macrophage)',
        'Erythroid',
        'Myeloid progenitor',
        'Adult tissues / Cardiac',
        'Adult tissues / Endocrine',
        'Thyroid',
        'Adrenal',
        'Digestive',
        'Adult tissues / Muscle',
        'Neural (Basal ganglia)',
        'Neural (Cereballar)',
        'Neural (Fetal)',
        'Stromal A',
        'Stromal B',
        'Endothelial',
        'Epithelial A',
        'Epithelial B',
        'Renal (Epithelial)',
        'Renal (Fetal)',
        'Respiratory (Fetal)',
        'Musculoskeletal (Fetal)',
        'Cardiac (Fetal)',
        'Placental A',
        'Placental B',
        'Embryonic',
        'Tissue invariant',
    ]

    if W_old is None:
        if os.path.exists(W_old_path) and W.shape[0] <= 29:
            W_old = np.load(W_old_path).T
            if W_old.shape[1] != W.shape[1]:
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

    component_data = pd.DataFrame({
        'index': component_order,
        'color': component_colors,
        'name': component_names,
    })

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
        myrandint = np.random.randint(len(colornames))
        while (count < n_components):
            new_color = colornames[np.random.randint(len(colornames))]
            trialcount = 0
            while ((new_color in component_colors) and (trialcount < 100)):
                #print('what am i doing here')
                newcolor = colornames[np.random.randint(0, len(colornames))]
                trialcount += 1
            #print('new color', count, new_color)
            component_colors.append(new_color)
            count += 1
    return component_colors[:n_components]