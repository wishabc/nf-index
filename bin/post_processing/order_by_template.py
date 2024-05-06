from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
import pandas as pd
import numpy as np


W_old_path = '/net/seq/data2/projects/sabramov/SuperIndex/dnase-index0415/embeddings/NMF/nndsvda_29c_050124//nmf_results/nndsvda.top35_max100_no_cancer_seed2674.peak_mask@cf_0.05@mv_0.7.29.W.npy'


def reorder_components(W1, W2):
    cosine_similarity = 1 - cdist(W1, W2, 'cosine')
    _, col_ind = linear_sum_assignment(-cosine_similarity)
    return col_ind


def get_component_data(W, W_old=None):
    if W_old is None:
        W_old = np.load(W_old_path)
    reorder = reorder_components(W, W_old)
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