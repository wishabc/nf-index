import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import umap

import matplotlib.colors as mcolors
import seaborn as sns

sys.path.append('/home/sboytsov/Vocabulary')
import OONMF
import OONMFhelpers
from OONMFhelpers import get_barsortorder

from common import get_matrix_path

def make_stacked_bar_plot_sorted(self, Nrelevant, BarMatrix, bargraph_out, names=[], barsortorder=[], plot_title='', figdim1=150, figdim2=40):
    if len(barsortorder) < 1:
        barsortorder = np.arange(Nrelevant)

    if len(names) < 1:
        names = [str(i) for i in range(Nrelevant)]
        names = np.array(names)

    start = 0
    end = Nrelevant
    Xpositions = np.arange(Nrelevant)

    self.define_colors()

    plt.clf()
    plt.figure(figsize=(figdim1, figdim2))

    for column in range(Nrelevant):
        ground_pSample = 0
        bar_heights = BarMatrix[:, start:end].T[barsortorder][column]
        indices = np.flip(np.argsort(bar_heights))
        # indices = np.arange(self.Ncomps)
        sorted_heights = bar_heights[indices]
        colours = [self.Comp_colors[i] for i in indices]

        for component in range(self.Ncomps):
            plt.bar(Xpositions[column], sorted_heights[component], bottom=ground_pSample, color=colours[component], alpha=1)
            ground_pSample = np.sum(sorted_heights[0:component + 1])

    OONMFhelpers.increase_axis_fontsize()

    plt.ylabel('sum of signal in matrix', fontsize=70)
    if (len(plot_title) > 0):
        plt.title(plot_title)

    samplenamesize = (1 / Nrelevant)**0.5 * 300
    thebottom = min([(1 / Nrelevant)**0.3 * 1.2, 0.3])

    plt.xticks(Xpositions, names[barsortorder], rotation='vertical', fontsize=samplenamesize)
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=thebottom)

    plt.savefig(bargraph_out)
    plt.show()
    plt.close()


OONMF.NMFobject.make_stacked_bar_plot_sorted = make_stacked_bar_plot_sorted


def get_colors_order(n_components):
    # code from Wouter's library, to have the similar colors as in barplot
    colors_order = ['#FFE500', '#FE8102', '#FF0000', '#07AF00', '#4C7D14', '#414613', '#05C1D9', '#0467FD', '#009588', '#BB2DD4', '#7A00FF', '#4A6876', '#08245B', '#B9461D', '#692108', '#C3C3C3']
    neworder = np.array([16,10,7,11,2,12,1,8,4,15,14,5,9,6,3,13]).astype(int) - 1
    colors_order = list(np.array(colors_order)[neworder])

    if (n_components>len(colors_order)):
        colornames = np.sort(list(mcolors.CSS4_COLORS.keys()))
        count = len(colors_order)
        np.random.seed(10)
        myrandint = np.random.randint(len(colornames))
        while (count < n_components):
            myrandint =    np.random.randint(len(colornames))
            newcolor = colornames[myrandint]
            trialcount = 0
            while ((newcolor in colors_order) and (trialcount < 100)):
                newcolor = colornames[np.random.randint(0,len(colornames))]
                trialcount+=1
            #print('new color ',count,newcolor)
            colors_order.append(newcolor)
            count+=1
    return colors_order


def load_meta(sample_order_path, meta_path, cluster_meta_path, gen_meta_path):
    # Sample order matrix to metadata 
    with open(sample_order_path) as f:
        indivs_order = np.array(f.readline().strip().split())

    # metadata with cluster names
    clust = pd.read_table(meta_path)
    clust = clust.set_index('id')

    # sorting metada in matrix's samples order
    metadata = pd.read_table(gen_meta_path)
    metadata = metadata.set_index('ag_id').loc[indivs_order]
    metadata = metadata.join(clust.loc[:, ['cluster']], how='left')

    # Here in column 'VA_cluster' I load Sasha's clustering labels
    z = pd.read_table(cluster_meta_path)
    metadata['VA_cluster'] = metadata['ln_number'].map(z.set_index('ln_number')['cluster'])

    metadata.reset_index(drop=True, inplace=True)

    return metadata


def vis_nmf(dir_path_mat, dir_path_pic, n_components, method, metadata, data_slice_size='full'):
    colors_order = get_colors_order(n_components)

    # Data loading
    allnames = metadata['taxonomy_name'].to_list()
    va_tags = metadata['VA_cluster'].to_list()

    decomp = OONMF.NMFobject(n_components)
    decomp.matrix_input_name(
        get_matrix_path(dir_path_mat, method, data_slice_size, n_components, 'W'),
        get_matrix_path(dir_path_mat, method, data_slice_size, n_components, 'H')
    )
    decomp.read_matrix_input(compressed=False)
    print(decomp.Basis.shape, decomp.Mixture.shape, "matrices loaded for n_components", n_components)
    bar_graph_sort_order = get_barsortorder(decomp.Basis)
    decomp.normalize_matrices()
    
    # heatmap
    decomp.make_standard_heatmap_plot(decomp.Basis.shape[0], 
                                decomp.Basis, 
                                dir_path_pic+str(method)+'.'+str(data_slice_size)+'.'+str(n_components)+'.heatmap.svg', 
                                names=np.array(allnames), 
                                barsortorder= bar_graph_sort_order)
    
    # Stacked barplot
    decomp.make_stacked_bar_plot_sorted(decomp.Basis.shape[0], 
                        decomp.Basis.T, 
                        dir_path_pic+str(method)+'.'+str(data_slice_size)+'.'+str(n_components)+'.stacked_bar_plot.svg', 
                        names=np.array(allnames), 
                        barsortorder=bar_graph_sort_order)

    # Normolized stcked barplot
    decomp.make_stacked_bar_plot_sorted(decomp.Basis.shape[0], 
                            decomp.NormedBasis.T, 
                            os.path.join(dir_path_pic, f'{method}.{data_slice_size}.{n_components}.stacked_bar_plot.normed.sorted.pdf'), 
                            names=np.array(allnames), 
                            barsortorder=bar_graph_sort_order)
    
    
    # UMAP 
    majcomp = np.argmax(decomp.Basis, axis=1)
    
    reducer = umap.UMAP(min_dist=0.5, n_neighbors=200, random_state=20)
    embedding = reducer.fit_transform(decomp.Basis)

    plt.clf()
    plt.figure(figsize=(10,10))
    decomp.define_colors()
    plt.scatter(embedding[:,0], embedding[:,1], color=np.array(decomp.Comp_colors)[majcomp], alpha=1, marker='.')
    plt.ylabel('UMAP axis 2')
    plt.xlabel('UMAP axis 1')
    plt.title(f'VST data, {method} Frobenius norm, {decomp.Mixture.shape[1]} DHSs')
    plt.savefig(dir_path_pic+str(method)+'.'+str(data_slice_size)+'.'+str(n_components)+'.umap.png', )
    plt.close(fig)
    
    # UMAP, comparing with VA and meta
    color_columns = ['VA_cluster', 'system']
    fig, axes = plt.subplots(1, 3, figsize=(30, 10))
    for idx, color_column in enumerate(color_columns):
        cluster_colors = sns.color_palette("tab10", len(metadata[color_column].unique()))
        color_dict = dict(zip(metadata[color_column].unique(), cluster_colors))
        metadata['cluster_color'] = metadata[color_column].map(color_dict)
        ax = axes[idx]
        for cluster, color in color_dict.items():
            mask = metadata[color_column] == cluster
            ax.scatter(embedding[mask, 0], embedding[mask, 1], color=color, label=cluster)
        ax.set_xlabel('UMAP axis 1')
        ax.set_ylabel('UMAP axis 2')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
    ax = axes[2]
    decomp.define_colors()
    ax.scatter(embedding[:, 0], embedding[:, 1], color=np.array(decomp.Comp_colors)[majcomp], alpha=1, marker='.')
    ax.set_ylabel('UMAP axis 2')
    ax.set_xlabel('UMAP axis 1')
    ax.set_title(f'VST data, {method} Frobenius norm, {decomp.Mixture.shape[1]} DHSs')

    plt.savefig(dir_path_pic+str(method)+'.'+str(data_slice_size)+'.'+str(n_components)+'.umap.comparing.png', bbox_inches='tight')
    plt.close(fig)

    ### Iterpretation part
    # Determine the number of rows and columns for the subplots
    nrows = int(np.ceil(np.sqrt(n_components)))
    ncols =int(np.ceil(n_components / nrows))

    # Create a figure with subplots
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(20, 10 * nrows))
    axes = axes.flatten()  # Flatten the axes array to easily loop through

    for i_component in range(n_components):
        # Plot i_component major
        green_cut = np.argmax(decomp.Basis, axis=1) == i_component
        print("Number of samples in this category ", len(green_cut[green_cut]))

        if len(green_cut[green_cut]) >= 1: # if not, there are too few samples to count 10% top
            # Plot top 10% of this major as described:
            colors = colors_order
            decomp_NormedBasis = decomp.NormedBasis[green_cut] 

            # Get the top 10% of samples
            #top_percentage = 0.1
            top_count = 10 #int(decomp_NormedBasis.shape[0] * top_percentage)

            # Sort the samples by the value of the i_component
            sorted_indices = np.argsort(decomp_NormedBasis[:, i_component])[-top_count:]
            sorted_matrix = decomp_NormedBasis[sorted_indices, :]
            sorted_names = np.array(allnames)[green_cut][sorted_indices]

            # Create a horizontal stacked bar plot on the subplot
            ax = axes[i_component]
            ax.barh(np.arange(sorted_indices.shape[0]), sorted_matrix[:, i_component], color=colors[i_component % len(colors)], label=f'Component {i_component}')

            # Add labels to the plot
            for i, name in enumerate(sorted_names):
                ax.text(sorted_matrix[i, i_component] / 2, i, name, ha='center', va='center', color='white', fontsize=20)

            # Customize the plot
            ax.set_yticks(np.arange(top_count))
            ax.set_yticklabels([])
            ax.legend(loc='upper right')
            ax.set_xlabel('Component values')
            #ax.set_ylabel('Top 10% samples')
            #ax.set_title(f'Horizontal bar plot for top 10% samples (Component {i_component})')

        else:
            continue

    # Remove unused subplots
    for i in range(i_component + 1, len(axes)):
        fig.delaxes(axes[i])

    # Adjust layout and show the figure
    #fig.tight_layout()
    plt.savefig(dir_path_pic+str(method)+'.'+str(data_slice_size)+'.'+str(n_components)+'.iterpretation.10barplot.png', bbox_inches='tight')
    plt.close(fig)


def main(sample_order_path, meta_path, cluster_meta_path, gen_meta_path, dir_path_mat, dir_path_pic, n_components, method):
    metadata = load_meta(sample_order_path, meta_path, cluster_meta_path, gen_meta_path)
    vis_nmf(dir_path_mat, dir_path_pic, int(n_components), method, metadata)


if __name__ == '__main__':
    main(sys.argv[1:])