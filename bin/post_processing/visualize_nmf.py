import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
from tqdm import tqdm

from order_by_template import get_component_data

def barplot_at_scale(matrix, metadata, order=None, component_data=None, tolerance=0.95, label_colors=None, agst=None):
    assert len(metadata) == matrix.shape[1]
    if component_data is None:
        raise NotImplementedError
    else:
        color_list = component_data['color']
        matrix = matrix[component_data['index'], :]
    matrix = matrix / matrix.sum(axis=0)
    if agst is None:
        agst = np.argsort(matrix, axis=0)[::-1, :]
    sorted_W = np.take_along_axis(matrix, agst, axis=0)
    
    sep = np.max(agst) + np.max(sorted_W) + 1
    if order is None:
        # order = np.argsort(agst[0, :] * sep**2 +  agst[1, :] * sep + sorted_W[0])[::-1]
        order = np.argsort(agst[0, :] * sep + sorted_W[0] / sorted_W.sum(axis=0))[::-1]
    
    colors = np.array(color_list)[agst[:, order]]
    
    heights = sorted_W[:, order]
    
    per_bar = 100
    chunks = matrix.shape[1] // per_bar + 1
    
    fig, axes = plt.subplots(chunks, 1, figsize=(20, 4*chunks))
    if chunks == 1:
        axes = [axes]
    fig.subplots_adjust(hspace=1.5)
    
    maxv = 0
    for k in tqdm(np.arange(chunks)):
        ax = axes[k]
        sl = slice(per_bar*k, per_bar*(k+1), 1)
        for j in np.arange(matrix.shape[1])[sl]:
            cumul_h = 0
            for i in np.arange(matrix.shape[0]):
                ax.bar(x=j, height=heights[i, j], bottom=cumul_h, color=colors[i, j])
                cumul_h += heights[i, j]
                maxv = max(maxv, cumul_h)
                if cumul_h >= tolerance:
                    break
        ax.set_xticks(np.arange(matrix.shape[1])[sl])
        ax.set_xticklabels(metadata.iloc[order, :][sl].apply(lambda row:
                                                             f"{row['core_ontology_term']} {row['SPOT1_score']:.1f}",
                                                             axis=1), rotation=90, )
        if label_colors is not None:
            assert len(label_colors) == matrix.shape[1]
            for xtick, col in zip(ax.get_xticklabels(), label_colors[sl]):
                xtick.set_color(col)
        ax.set_xlim(np.arange(matrix.shape[1])[sl].min(), np.arange(matrix.shape[1])[sl].min() + per_bar)
        # break
    
    for ax in axes:
        ax.set_ylim(0, maxv*1.05)

    return order, fig

def get_barsortorder(matrix):
    # assumes rows are the data, columns are NMF components
    max_component_idx = np.argmax(matrix, axis=1)
    barsortorder = []
    for i in range(matrix.shape[1]):
        desired_order = np.argsort(-matrix[:,i])
        relevant_cut = max_component_idx[desired_order] == i
        barsortorder.append(desired_order[relevant_cut])
    matrix = matrix.astype(int)
    return np.concatenate(barsortorder)


def plot_component(data, labels, color, ax=None, top_count=15):
    n_samples = data.shape[0]
    top_count = min(top_count, n_samples)  # Adjust to plot up to 15 samples, or fewer if less are available

    sorted_indices = np.argsort(data)[-top_count:]
    sorted_data = data[sorted_indices]
    sorted_names = labels[sorted_indices]

    ax.barh(np.arange(top_count), sorted_data, color=color)
    ax.set_yticks(np.arange(top_count))
    ax.set_yticklabels(sorted_names, ha='right', va='center', fontsize=20, color='k')
    ax.set_xlim(0, 1.1 * sorted_data.max())
    return ax
    

def plot_nmf(basis, annotations, top_count=15, component_data=None, common_scale=False):
    n_components = basis.shape[0]

    #bar_graph_sort_order = get_barsortorder(basis)

    ncols = int(np.ceil(np.sqrt(n_components)))
    nrows = int(np.ceil(n_components / ncols))
    
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10 * ncols, 7 * nrows))
    
    axes = axes.flatten()
    assert len(axes) >= n_components
    
    [fig.delaxes(ax) for ax in axes[n_components:]]

    xlims = []
    for ax, (_, row) in zip(axes[:n_components], component_data.iterrows()):
        component_is_major = np.argmax(basis, axis=0) == row['index']
        
        ax = plot_component(
            basis[row['index'], component_is_major],
            annotations[component_is_major], 
            color=row['color'], 
            ax=ax,
            top_count=top_count
        )
        ax.set_title(f'{row["name"]}')
        xlims.append(ax.get_xlim())
    if common_scale:
        for ax in axes:
            ax.set_xlim(0, max(xlims, key=lambda x: x[1])[1])
    plt.tight_layout()
    return fig, axes


def get_order(agst, sorted_matrix, by='primary'):
    sep = np.max(agst) + np.max(sorted_matrix) + 1
    if by == 'secondary':
        order = np.argsort(agst[0, :] * sep**2 +  agst[1, :] * sep + sorted_matrix[0])[::-1]
    elif by == 'primary':
        order = np.argsort(agst[0, :] * sep + sorted_matrix[0])[::-1]
    return order


def get_tops_and_bottoms(agst, heights):
    tops = heights.cumsum(axis=0)
    bottoms = tops - heights
    
    idxs = np.argsort(agst, axis=0)
    return np.take_along_axis(tops, idxs, axis=0), np.take_along_axis(bottoms, idxs, axis=0)


def plot_stacked(matrix, colors, ax=None, lims=None, order_by='primary', order=None, agst=None):
    if lims is None:
        lims = 0, matrix.shape[1]
    matrix = matrix[:, lims[0]:lims[1]]
    n_components, n_objects = matrix.shape

    if agst is None:
        agst = np.argsort(matrix, axis=0)[::-1, :]
    heights = np.take_along_axis(matrix, agst, axis=0)

    if order is None:
        order = get_order(agst, heights, by=order_by)

    tops, bottoms = get_tops_and_bottoms(agst[:, order], heights[:, order])

    fb_tops = np.repeat(tops, 2, axis=1)
    fb_bottoms = np.repeat(bottoms, 2, axis=1)
    xvals = np.concatenate([[0], np.repeat(np.arange(1, matrix.shape[1]), 2), [matrix.shape[1]]])

    if ax is None:
        fig, ax = plt.subplots(figsize=(matrix.shape[1]/100, 2))
    for i, color in enumerate(colors):
        ax.fill_between(xvals, fb_bottoms[i], fb_tops[i], lw=0, color=color)
    ax.set_xlim(-0.5, matrix.shape[1] + 0.5)

    return ax, agst, order


def plot_barplots(matrix, component_data=None, n=10_000, order_by='primary', order=None, agst=None, normalize=True, ax=None):
    if matrix.shape[1] > n:
        np.random.seed(0)
        H_dsp = matrix[:, np.random.choice(np.arange(matrix.shape[1]), n)]
    else:
        H_dsp = matrix

    if component_data is None:
        raise ValueError
    else:
        H_dsp = H_dsp[component_data['index'], :]
        colors = component_data['color']

    if normalize:
        H_dsp = H_dsp / H_dsp.sum(axis=0)

    if ax is None:
        fig, ax = plt.subplots(figsize=(20, 2))

    return plot_stacked(H_dsp, colors, ax=ax, order_by=order_by, order=order, agst=agst)


def adjust_agst(agst, i):
    n_rows, n_cols = agst.shape
    # Create a new array that will hold the adjusted indices
    new_agst = np.empty_like(agst)
    
    # For each column, find the position of the zero index
    zero_positions = (agst == i)
    
    # Iterate through each column to set the zero at the top and reorder the remaining elements
    for col in range(n_cols):
        col_indices = agst[:, col]
        zero_pos = np.where(col_indices == i)[0][0]
        # Place the zero at the top of the new column
        new_agst[0, col] = col_indices[zero_pos]
        # Fill the rest of the column with the remaining indices, excluding the original zero position
        new_agst[1:, col] = np.concatenate([col_indices[:zero_pos], col_indices[zero_pos + 1:]])
        
    return new_agst


def main(binary_matrix, W, H, metadata, samples_mask, peaks_mask, dhs_annotations, vis_path):
    component_data = get_component_data(W)

    print('Order samples by component contribution')
    relative_W = W / W.sum(axis=0)
    for i, row in component_data.iterrows():
        weights = np.ones(W.shape[0])
        weights[i] *= W.shape[0] * 2

        agst = np.argsort(relative_W * weights[:, None], axis=0)[::-1, :]
    
        _, fig = barplot_at_scale(
            W,
            metadata,
            component_data=component_data,
            order=np.argsort(-relative_W[row['index'], :]),
            agst=agst
        )
        comp_name = row["name"].replace("/", "_")
        fig.savefig(f'{vis_path}/detailed_barplot_all_normal_samples.{comp_name}.pdf', transparent=True, bbox_inches='tight')
        plt.close(fig)
    
       # Plot samples
    print('All samples')
    ax, _, _ = plot_barplots(W, component_data)
    plt.savefig(f'{vis_path}/Barplot_all_normal_samples.pdf', transparent=True, bbox_inches='tight')
    plt.close(ax.get_figure())

    print('Train stratified samples set')
    ax, _, _ = plot_barplots(W[:, samples_mask], component_data)
    plt.savefig(f'{vis_path}/Barplot_stratified_train_samples.pdf', transparent=True, bbox_inches='tight')
    plt.close(ax.get_figure())

    # Plot peaks
    print('Reference peaks')
    ax, _, _ = plot_barplots(H[:, peaks_mask], component_data, normalize=True)
    plt.savefig(f'{vis_path}/Barplot_reference_train_DHSs.pdf', transparent=True, bbox_inches='tight')
    plt.close(ax.get_figure())

    print('All peaks')
    ax, _, _ = plot_barplots(H, component_data)
    plt.savefig(f'{vis_path}/Barplot_all_DHSs.pdf', transparent=True, bbox_inches='tight')
    plt.close(ax.get_figure())

    print('All peaks not normalized')
    ax, _, _ = plot_barplots(H, component_data, normalize=False)
    plt.savefig(f'{vis_path}/Barplot_all_DHSs.not_norm.pdf', transparent=True, bbox_inches='tight')
    plt.close(ax.get_figure())

    #Only reproduced DHSs
    print('>3 peaks supproting a DHS')
    reproduced_peaks = binary_matrix.sum(axis=1) > 3
    ax, _, _ = plot_barplots(H[:, reproduced_peaks], component_data, normalize=True)
    plt.savefig(f'{vis_path}/Barplot_DHS_supported_by_4+samples.pdf', transparent=True, bbox_inches='tight')
    plt.close(ax.get_figure())

    print('Detailed barplot all samples')
    s_order, fig = barplot_at_scale(W, metadata, component_data=component_data)
    if samples_mask.sum() < samples_mask.shape[0]:
        plt.close(fig)
        s_mask = samples_mask[s_order]
        barplot_at_scale(
            W,
            metadata,
            s_order,
            component_data=component_data,
            label_colors=[
                'r' if s else 'k' for s in s_mask
            ]
        )
    plt.savefig(f'{vis_path}/detailed_barplot_all_normal_samples.pdf', transparent=True, bbox_inches='tight')
    plt.close(fig)

    print('Top 20 samples per component')
    annotations = metadata["taxonomy_name"].values
    fig, axes = plot_nmf(W, annotations, top_count=20, component_data=component_data)
    plt.savefig(f'{vis_path}/Top20_all_samples_barplot.pdf', bbox_inches='tight', transparent=True)
    plt.close(fig)

    fig, axes = plot_nmf(W, annotations, top_count=20, component_data=component_data, common_scale=True)
    plt.savefig(f'{vis_path}/Top20_all_samples_barplot.common_scale.pdf', bbox_inches='tight', transparent=True)
    plt.close(fig)


    if dhs_annotations is not None:
        ax = plot_dist_tss(H, dhs_annotations, component_data)
        plt.savefig(f'{vis_path}/Distance_to_tss.pdf', bbox_inches='tight', transparent=True)
        plt.close(fig)



def plot_dist_tss(H, index, component_data, ax=None):
    max_component = np.argmax(H, axis=0)
    if ax is None:
        fig, ax = plt.subplots(figsize=(2, 2))
    for i, row in component_data.iterrows():
        data = np.abs(index['dist_tss'][max_component == row['index']])
        ax.plot(np.sort(data), np.linspace(0, 1, len(data)), color=row['color'])
    ax.set_xlim(-50, 5000)
    ax.set_xlabel('Distance to TSS')
    ax.set_ylabel('Cumulative proportion of DHSs')
    return ax


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Matrix normalization using lowess')
    parser.add_argument('matrix', help='Path to matrix to run NMF on')
    parser.add_argument('W', help='W matrix of perform NMF decomposition')
    parser.add_argument('H', help='H matrix of perform NMF decomposition')
    parser.add_argument('metadata', help='Path to metadata file')
    parser.add_argument('dhs_meta', help='Path to DHS index')
    parser.add_argument('n_components', help='Number of components to use in NMF', type=int)
    parser.add_argument('--samples_mask', help='Mask of used samples, numpy array')
    parser.add_argument('--peaks_mask', help='Mask of used samples, numpy array')
    #parser.add_argument('--samples_weights', help='Path to samples weights (for weighted NMF)', default=None)
    parser.add_argument('--outpath', help='Path to save visualizations', default='./')
    parser.add_argument('--dhs_annotations', help='Path to DHS annotations', default=None)
    args = parser.parse_args()

    mat = np.load(args.matrix).astype(float)
    samples_m = np.loadtxt(args.samples_mask).astype(bool)
    peaks_m = np.loadtxt(args.peaks_mask).astype(bool)

    metadata = pd.read_table(args.metadata)
    W = np.load(args.W).T
    H = np.load(args.H).T
    dhs_meta = pd.read_table(args.dhs_meta, header=None)[3].to_numpy()
    if args.dhs_annotations is not None:
        dhs_annotations = pd.read_table(args.dhs_annotations)
        dhs_annotations = dhs_annotations[dhs_annotations['dhs_id'].isin(dhs_meta)]
    else:
        dhs_annotations = None

    main(mat, W, H, metadata, samples_m, peaks_m, dhs_annotations, args.outpath)

