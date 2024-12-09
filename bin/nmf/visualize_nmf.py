import numpy as np
import matplotlib.pyplot as plt
import argparse
from tqdm import tqdm

from order_by_template import get_component_data, define_colors
from perform_NMF import NMFInputData, parse_nmf_args, setup_parser


def barplot_at_scale(matrix, metadata, colors, order=None, agst=None, label_colors=None):
    assert len(metadata) == matrix.shape[1], f'Metadata ({len(metadata)}) should have the same number of samples as the matrix ({matrix.shape[1]})'
    
    if agst is None:
        agst = np.argsort(matrix, axis=0)[::-1, :]
    if order is None:
        sep = np.max(agst) + np.max(matrix) + 1
        max_r = matrix.max(axis=0)
        order = np.argsort(agst[0, :] * sep + max_r / matrix.sum(axis=0))[::-1]

    ordered_matrix = matrix[:, order]
    
    per_bar = 100
    chunks = np.ceil(matrix.shape[1] / per_bar).astype(int)
    
    fig, axes = plt.subplots(chunks, 1, figsize=(20, 4*chunks))
    if chunks == 1:
        axes = [axes]
    fig.subplots_adjust(hspace=1.5)
    
    maxv = np.max(matrix.sum(axis=0))
    for k in tqdm(np.arange(chunks)):
        ax = axes[k]
        sl = slice(per_bar*k, per_bar*(k+1), 1)
        num_elements = order[sl].shape[0]
        plot_stacked(ordered_matrix[:, sl], colors, ax=ax,
                     order=np.arange(num_elements),
                     agst=agst[:, order[sl]])
        ax.set_xticks(np.arange(num_elements) + 0.5)
        ax.set_xticklabels(
            metadata.iloc[order, :]['sample_label'][sl],
            rotation=90
        )
        if label_colors is not None:
            assert len(label_colors) == matrix.shape[1]
            for xtick, col in zip(ax.get_xticklabels(), label_colors[sl]):
                xtick.set_color(col)
        ax.set_xlim(0, per_bar)
    
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
    top_count_actual = min(top_count, n_samples)

    sorted_indices = np.argsort(data)[-top_count_actual:]
    sorted_data = data[sorted_indices]
    sorted_names = labels[sorted_indices]

    ax.barh(np.arange(top_count_actual), sorted_data, color=color)
    ax.set_yticks(np.arange(top_count_actual))
    ax.set_ylim(top_count_actual - top_count - 0.5, top_count_actual - 0.5)
    ax.set_yticklabels(sorted_names, ha='right', va='center', fontsize=20, color='k')
    if len(sorted_data) > 0:
        ax.set_xlim(0, 1.1 * sorted_data.max())
    return ax
    

def plot_nmf(basis, annotations, top_count=15, component_data=None, common_scale=False):
    n_components = basis.shape[0]

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


def plot_stacked(matrix, colors, ax=None, lims=None, order_by='primary', order=None, agst=None, orient='horizontal'):
    if lims is None:
        lims = 0, matrix.shape[1]
    matrix = matrix[:, lims[0]:lims[1]]

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
        if orient == 'horizontal':
            ax.fill_between(xvals, fb_bottoms[i], fb_tops[i], lw=0, color=color)
            ax.set_xlim(0, matrix.shape[1])
            ax.set_ylim(0, None)
        elif orient == 'vertical':
            ax.fill_betweenx(xvals, fb_bottoms[i], fb_tops[i], lw=0, color=color)
            ax.set_ylim(0, matrix.shape[1])
            ax.set_xlim(0, None)

    return ax, agst, order


def plot_barplots(matrix, component_data=None, n=10_000, normalize=True, ax=None, **kwargs):
    if matrix.shape[1] > n:
        np.random.seed(0)
        H_dsp = matrix[:, np.random.choice(np.arange(matrix.shape[1]), n)]
    else:
        H_dsp = matrix

    if component_data is None:
        colors = define_colors(H_dsp.shape[0])
    else:
        H_dsp = H_dsp[component_data['index'], :]
        colors = component_data['color']

    if normalize:
        H_dsp = H_dsp / H_dsp.sum(axis=0)

    if ax is None:
        fig, ax = plt.subplots(figsize=(20, 2))

    return plot_stacked(H_dsp, colors, ax=ax, **kwargs)


def main(
        nmf_data: NMFInputData,
        W: np.ndarray,
        H: np.ndarray,
        vis_path,
        component_data=None
    ):
    if component_data is None:
        component_data = get_component_data(W)

    print('Order samples by component contribution')
    relative_W = W / W.sum(axis=0)
    relative_W = relative_W[component_data['index'], :]

    metadata = nmf_data.samples_metadata
    samples_mask = nmf_data.samples_mask
    peaks_mask = nmf_data.peaks_mask

    binary_matrix = nmf_data.matrix[peaks_mask, :]
    dhs_meta = nmf_data.dhs_metadata

    # for i, row in component_data.iterrows():
    #     weights = np.ones(W.shape[0])
    #     weights[i] = W.shape[0]

    #     agst = np.argsort(relative_W + weights[:, None], axis=0)[::-1, :]
    #     print((agst[0] == i).all())

    #     _, fig = barplot_at_scale(
    #         relative_W,
    #         metadata.iloc[:, :],
    #         colors=component_data['color'],
    #         order=np.argsort(relative_W[i, :])[::-1],
    #         agst=agst
    #     )
    #     comp_name = row["name"].replace("/", "_").replace(' ', '_')
    #     fig.savefig(f'{vis_path}/Detailed_barplot.{comp_name}.pdf', transparent=True, bbox_inches='tight')
    #     plt.close(fig)

    ######### Plot samples #########
    if samples_mask.sum() < samples_mask.shape[0]:
        print('Reference samples set')
        ax, _, _ = plot_barplots(W[:, samples_mask], component_data)
        plt.savefig(f'{vis_path}.Barplot_reference_train_samples.pdf', transparent=True, bbox_inches='tight')
        plt.close(ax.get_figure())
    

    print('All samples')
    ax, _, _ = plot_barplots(W, component_data)
    plt.savefig(f'{vis_path}.Barplot_all_samples.pdf', transparent=True, bbox_inches='tight')
    plt.close(ax.get_figure())


    ######## Plot peaks #########
    print('All peaks')
    ax, _, _ = plot_barplots(H, component_data)
    plt.savefig(f'{vis_path}.Barplot_all_DHSs.pdf', transparent=True, bbox_inches='tight')
    plt.close(ax.get_figure())

    print('All peaks not normalized')
    ax, _, _ = plot_barplots(H, component_data, normalize=False)
    plt.savefig(f'{vis_path}.Barplot_all_DHSs.not_norm.pdf', transparent=True, bbox_inches='tight')
    plt.close(ax.get_figure())

    #Only reproduced DHSs
    print('>3 peaks supproting a DHS')
    reproduced_peaks = binary_matrix.sum(axis=1) > 3
    ax, _, _ = plot_barplots(H[:, reproduced_peaks], component_data, normalize=True)
    plt.savefig(f'{vis_path}.Barplot_DHS_supported_by_4+samples.pdf', transparent=True, bbox_inches='tight')
    plt.close(ax.get_figure())

    print('Detailed barplot all samples')
    s_order, fig = barplot_at_scale(relative_W, metadata, colors=component_data['color'])
    if samples_mask.sum() < samples_mask.shape[0]:
        plt.close(fig)
        s_mask = samples_mask[s_order]
        barplot_at_scale(
            relative_W,
            metadata,
            colors=component_data['color'],
            order=s_order,
            label_colors=[
                'r' if s else 'k' for s in s_mask
            ]
        )
    plt.savefig(f'{vis_path}.Detailed_barplot_all_samples.pdf', transparent=True, bbox_inches='tight')
    plt.close(fig)

    print('Top 20 samples per component')
    annotations = metadata["sample_label"].values
    fig, axes = plot_nmf(W, annotations, top_count=20, component_data=component_data)
    plt.savefig(f'{vis_path}.Top20_all_samples_barplot.pdf', bbox_inches='tight', transparent=True)
    plt.close(fig)

    fig, axes = plot_nmf(W, annotations, top_count=20, component_data=component_data, common_scale=True)
    plt.savefig(f'{vis_path}.Top20_all_samples_barplot.common_scale.pdf', bbox_inches='tight', transparent=True)
    plt.close(fig)


    if 'dist_tss' in dhs_meta.columns:
        ax = plot_dist_tss(H, dhs_meta['dist_tss'], component_data)
        plt.savefig(f'{vis_path}.Distance_to_tss.pdf', bbox_inches='tight', transparent=True)
        plt.close(fig)



def plot_dist_tss(H, dist_tss, component_data, ax=None):
    max_component = np.argmax(H, axis=0)
    if ax is None:
        fig, ax = plt.subplots(figsize=(2, 2))
    for i, row in component_data.iterrows():
        data = np.abs(dist_tss[max_component == row['index']])
        ax.plot(np.sort(data), np.linspace(0, 1, len(data)), color=row['color'])
    ax.set_xlim(-50, 5000)
    ax.set_xlabel('Distance to TSS')
    ax.set_ylabel('Cumulative proportion of DHSs')
    return ax


if __name__ == '__main__':
    print('Visualizing NMF results')
    parser = argparse.ArgumentParser('Matrix normalization using lowess')
    parser = setup_parser()
    print('Adding options to parser')
    parser.add_argument('W', help='W matrix of perform NMF decomposition')
    parser.add_argument('H', help='H matrix of perform NMF decomposition')
    parser.add_argument('--outpath', help='Path to save visualizations', default='./')
    args = parser.parse_args()

    nmf_data = parse_nmf_args(args)
    
    if 'sample_label' not in nmf_data.samples_metadata.columns:
        if 'SPOT1_score' not in nmf_data.samples_metadata.columns and 'core_ontology_term' not in nmf_data.samples_metadata.columns:
            raise ValueError('sample_label column is missing in samples metadata')
        print('Assuming sample_label is core_ontology_term + SPOT1_score')
        nmf_data.samples_metadata.loc[:, 'sample_label'] = (
            nmf_data.samples_metadata['core_ontology_term'].astype(str) + " " + 
            nmf_data.samples_metadata['SPOT1_score'].astype(str).str[:3]
        )

    W = np.load(args.W).T # NMF components x samples
    H = np.load(args.H).T # NMF components x peaks
    outprefix = f"{args.outpath}/{args.prefix}"
    main(nmf_data, W, H, outprefix)

