#!/usr/bin/env python
# coding: utf-8

# # scVelo Human Cattaneo Data

# Data has been integrated using Seurat. Next, Velocities will be estimated using scVelo package and an AnnData Object will be created.
# - Read TSV file
# - Compute velocities
# - Visualization 2D UMAP plots

# ## Dependancies and functions

# In[1]:


# velocity version
vlo_ver="v1"


# In[2]:


import matplotlib.pyplot as plt
import scvelo as scv # 0.2.5
import scipy.stats as st
import tqdm
import loompy as lmp
import scanpy as sc
import pandas as pd
import numpy as np
import ipywidgets
import warnings
import anndata
import igraph
import tqdm
import sys
import os
import re

warnings.filterwarnings('ignore')


# In[3]:


np. __version__ 


# In[4]:


pd.__version__


# In[5]:


scv.__version__


# In[6]:


scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.settings.set_figure_params('scvelo')  # for beautified visualization


# In[7]:


def read_in(targProj,targSamp,seuTab):

    # define base dir (si cattaneo o cardif)
    bsedirProj="/gpfs/projects/bsc83/MN4/bsc83/Projects/Creatio/" + targProj + "/"
    
    # get sample table
    smpTab = seuTab.loc[seuTab["sample"]==targSamp]
    smpTab['cell'] = smpTab['cell'].apply(lambda x: re.sub(r"-.*$", "-1", x))
    smpTab.set_index('cell', inplace=True)
    
    # recover original cell name (needs to be generalized)
    
    # Define Path to cellranger output, targSamp being the name of the sample
    Path10x= bsedirProj + 'cellranger_new/' + targSamp + '/outs/filtered_feature_bc_matrix/'
    
    # Read filtered feature bc matrix output from cellranger count
    scData = sc.read_10x_mtx(Path10x,var_names='gene_symbols',cache=True)
    
    # Filter out existing row index labels
    # print(scData.obs.index.tolist())
    # print(smpTab.index)
    valid_cells = [row for row in scData.obs.index.tolist() if row in smpTab.index]
    
    # keep only processed cells
    scData = scData[valid_cells]
    
    # get cell annotacion for cells
    scData.obs["clusters"] = smpTab.loc[valid_cells,'cell_annot'].astype('category')

    # get stage
    scData.obs["phase"] = smpTab.loc[valid_cells,'phase'].astype('category')

    # get sample
    scData.obs["sample"] = smpTab.loc[valid_cells,'sample'].astype('category')
    
    # get pca for cells
    scData.obsm["X_pca"] = smpTab.loc[valid_cells,["PC_" + str(i) for i in range(1, 10)]].to_numpy()

    # get umap for cells
    scData.obsm["X_umap"] = smpTab.loc[valid_cells,["UMAP_1","UMAP_2"]].to_numpy()
    
    # Replace cell names in the AnnData object
    scData.obs_names = [targSamp + ":" + cell_name.replace('-1','') + "x" for cell_name in scData.obs_names]

    # return empty anndata if no cells
    if scData.n_obs==0:
        print("Sample " + targSamp + " had no cells")
        return anndata.AnnData()
    
    # Read velocyto output
    smpVelo = scv.read(bsedirProj + 'cellranger_new/' + targSamp + '/velocyto/' + targSamp + '.loom', cache=True)
    
    # keep only processed cells
    smpVelo = smpVelo[scData.obs_names]
    
    # Merge velocyto and cellranger outputs
    scData = scv.utils.merge(scData, smpVelo)

    # return object
    return scData

def run_scvelo(x,adataFile):
    
    if os.path.exists(adataFile):
        print("Previous scVelo results were found: "+ adataFile)
        return scv.read(adataFile)
    else:
        print("Running scvelo and storing output here: "+ adataFile)
        
        # filter, normalize and log transform
        scv.pp.filter_genes(x)
        scv.pp.normalize_per_cell(x)
        scv.pp.filter_genes_dispersion(x)
        scv.pp.log1p(x)

        # estimate moments
        scv.pp.moments(x, n_pcs=10, n_neighbors=15)
        
        # recover dynamics
        scv.tl.recover_dynamics(x, n_jobs=120)
        
        #  to project the velocities into a lower-dimensional embedding, transition 
        # probabilities of cell-to-cell transitions are estimated
        scv.tl.velocity(x, n_jobs=120, mode='dynamical', filter_genes=True)
        scv.tl.velocity_graph(x, n_jobs=120)
        
        # estimate pseudotime
        scv.tl.recover_latent_time(x)

        # get root and terminal states
        scv.tl.terminal_states(x)

        # psuedotime
        scv.tl.velocity_pseudotime(x)

        # velocity confidence
        scv.tl.velocity_confidence(x)

        # store
        x.write(adataFile)

        # return updated object
        return x


# # Read in

# In[8]:


# base directories
bsedirInt="/gpfs/projects/bsc83/MN4/bsc83/Projects/Creatio/integration_human/v6/scvelo/"
outdir=f"/gpfs/projects/bsc83/MN4/bsc83/Projects/Creatio/integration_human/v6/scvelo/"

# mkdir outdir
if not os.path.exists(outdir):
    os.makedirs(outdir)
    
# Read in seurat table
seuTab = bsedirInt + "seurat_export_for_scvelo_v6.tsv"
seuTab = pd.read_csv(seuTab, sep='\t')

# read in metadata from cattaneo
bsedirProjCat="/gpfs/projects/bsc83/MN4/bsc83/Projects/Creatio/CattaneoPaper/cellranger_new/"
metaCat = pd.read_csv(bsedirProjCat+"sample_names.txt", header=None)
excCat = ["ERR3988191","ERR3988192","ERR3988193","ERR3988194","ERR3988195","ERR3988196", "ERR3988197", "ERR3988198", "ERR3988147", "ERR3988148", "ERR3988149", "ERR3988150", "ERR3988181", "ERR3988190"]
#excCat=["ERR3988181", "ERR3988190"] no hi ha la folder de outs
metaCat = metaCat[~metaCat[0].isin(excCat)] 

# Define the file name
fn = f"{outdir}/adata_unnorm_{vlo_ver}.h5ad"

# Check if the file exists
if os.path.exists(fn):
    # Load the file if it exists
    adata = anndata.read_h5ad(fn)
    print(f"Loaded existing file: {fn}")
else:
    ## read in data
    scDataLstCat = []
    for targSamp in tqdm.tqdm(metaCat[0].tolist(), desc="Processing samples"):
        scDataLstCat.append(read_in("CattaneoPaper", targSamp, seuTab))
    
    ## concatenate samples
    scDataLstCat = [x for x in scDataLstCat if x.n_obs > 0]
    adataCat = anndata.concat(scDataLstCat)
    adata = adataCat
    
    # make clusters categorical
    adata.obs["clusters"] = adata.obs["clusters"].astype('category')
    
    # store anndata file with unnormalized matrix
    adata.write(fn)


# Plot UMAP just to make sure everything has been imported correctly.

# In[9]:


scv.pl.scatter(adata, color='clusters', cmap='gnuplot', show=True, alpha=0.5, legend_loc='right margin')


# # RNA velocity

# In[10]:


fn=f"/gpfs/projects/bsc83/MN4/bsc83/Projects/Creatio/integration_human/v6/scvelo/adata_{vlo_ver}_scvelo.h5ad"

if os.path.exists(fn):
    adata = anndata.read_h5ad(fn)
else:
    run_scvelo(adata,fn)


# # Velocity embeddings 

# We can embed the velocity on the UMAP

# In[11]:


# Create a figure with subplots
fig, ax = plt.subplots(figsize=(7, 5))

# plot velocities
scv.pl.velocity_embedding_stream(adata, basis='umap', show=True, ax=ax)
ax.set_title('Embedding')

# Adjust the layout of the subplots
plt.tight_layout()

# add title
fig.suptitle("Velocity embeddings")

# Save the combined plot
#plt.savefig(f"{adata_dir}/umap_velo_emb_{vlo_ver}_filtered.png")
plt.show()


# # Root cells and end points

# In scvelo, every cell gets assigned a root cell probability and an end point probability. We can plot these and compare them to identify difference in the location of the both types of cells.

# In[12]:


# Create a figure with subplots
adata_dir= "/gpfs/projects/bsc83/MN4/bsc83/Projects/Creatio/integration_human/v6/scvelo/"
fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(20, 10))

# plot root and end point probability
scv.pl.scatter(adata, color=['root_cells'], cmap='viridis', show=False, ax=ax1)
scv.pl.scatter(adata, color=['end_points'], cmap='viridis', show=False, ax=ax2)

# Set titles for the subplots
ax1.set_title('Root')
ax2.set_title('Ends')

# Adjust the layout of the subplots
plt.tight_layout()

# Save the combined plot
#plt.savefig(f"{adata_dir}figures/umap_root_ends{vlo_ver}_new.png")
plt.show()
plt.close()


# # Pseudotime

# Velocity pseudotime is a random-walk based distance measures on the velocity graph. After computing a distribution over root cells obtained from the velocity-inferred transition matrix, it measures the average number of steps it takes to reach a cell after start walking from one of the root cells. Contrarily to diffusion pseudotime, it implicitly infers the root cells and is based on the directed velocity graph instead of the similarity-based diffusion kernel.

# In[13]:


scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot', show=True)  
#plt.savefig(f"{adata_dir}figures/umap_pseudotime_{vlo_ver}_new.svg")


# In[16]:


adata


# In[17]:


scv.pl.scatter(adata, color='latent_time', cmap='gnuplot', show=True)  


# # Top genes

# We can check which genes are driving the differentiation process by finding genes with high correlation with pseudotime.

# In[19]:


# genes driving differentiation
top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:50]
scv.pl.heatmap(adata, var_names=top_genes, sortby='latent_time', col_color='latent_time', n_convolve=100)


# # Cluster markers (differential velocity)

# We can test which genes have **cluster-specific differential velocity expression**, being significantly higher/lower compared to the remaining population. The module scv.tl.rank_velocity_genes runs a differential velocity t-test and outputs a gene ranking for each cluster. Thresholds can be set (e.g. min_corr) to restrict the test on a selection of gene candidates.

# In[20]:


# Identify important genes
scv.tl.rank_velocity_genes(adata, groupby='clusters', min_corr=.3)


# # Velocity length and confidence

# We can also compare the velocity length. However, we have to be careful here since we might not be comparing the same scale...

# In[21]:


# Create a figure with subplots
fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(20, 10))

# plot velocity length and confidence
scv.pl.scatter(adata, c='velocity_length', cmap='coolwarm', perc=[5, 95], show=False, ax=ax1)
scv.pl.scatter(adata, c='velocity_confidence', cmap='coolwarm', perc=[5, 95], show=False, ax=ax2)

# Set titles for the subplots
ax1.set_title('Velocity')
ax2.set_title('Confidence')

# Adjust the layout of the subplots
plt.tight_layout()

# add title
fig.suptitle("Velocity and confidence")

# Save the combined plot
#plt.savefig(f"{adata_dir}figures/umap_velo_len_conf_{vlo_ver}.png")
plt.show()
plt.close()

