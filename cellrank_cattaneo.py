#!/usr/bin/env python
# coding: utf-8

# # CellRank Human Cattaneo Data
# 
# Data has been integrated using Seurat package and RNA velocities have been computed using Velocyto and scVelo. 
# We have AnnData as input with the RNA velocities results.
# 
# - 2D UMAp plot for RNA velocity Visualization
# - CellRank V2 Velocity Kernel and connectivity kernel are used to compute transition matrixes
# - Combined kernels
# - Visualization of transition matrix using low dimensional embedding velocity streamline plot and Random walk
# - Initial and terminal states computation using advanced shur decomposition and GPCCA Estimator

# ## Import packages

# In[ ]:


import os
import sys

# import standard packages
import numpy as np
import pandas as pd 

# Plotting packages
import matplotlib.pyplot as pl
import seaborn as sns
import plotly.graph_objects as go  # for 3D plotting of UMAP and velocities
import plotly.express as px

#import single-cell packages
import scanpy as sc
import scanpy.external as sce
import scvelo as scv
import cellrank as cr

# other Packages
# import scipy.stats as st
from scipy.stats import ranksums 
from scipy.sparse import find # use for removing outliers in UMAP
import loompy as lmp
import ipywidgets as widgets
from IPython.display import display
import warnings
import anndata as ad
import igraph
import tqdm
import leidenalg # for leidian algorithm 


# Optimization Packages
import petsc
import slepc
import petsc4py
import slepc4py
import jax
import rpy2

warnings.filterwarnings('ignore')
print(sys. version)


# ## Set global parameters

# In[4]:


# show errors(0), warnings(1), info(2), hints(3)
# set verbosity levels
sc.settings.verbosity = 3
cr.settings.verbosity = 3
scv.settings.verbosity = 3 

scv.settings.presenter_view = True  # set max width size for presenter view
scv.settings.set_figure_params('scvelo', dpi_save=100, dpi=80, transparent=False, fontsize=8, color_map='viridis')
scv.settings.plot_prefix = ""


# In[5]:


# should figures just be displayed or also saved?
save_figure = True


# In[6]:


from typing import Optional, Iterable, TypeVar, Tuple
AnnData = TypeVar("AnnData")


# ## Set directories

# In[30]:


out_dir="/gpfs/projects/bsc83/MN4/bsc83/Projects/Creatio/integration_human/v6/cell_fate"


# ## Read Anndata Object

# In[10]:


#Read anndata object 
adata = sc.read("/gpfs/projects/bsc83/MN4/bsc83/Projects/Creatio/integration_human/v6/scvelo/adata_v1_scvelo.h5ad")


# ## Explore adata object

# In[11]:


adata


# In[12]:


adata.layers['velocity']


# In[14]:


# Plot velocities avaliable in adata object UMAP
scv.pl.velocity_embedding_stream(adata, 
                                 basis='umap',
                                 density=3,
                                 title='RNA Velocity',
                                 fontsize=15, 
                                 legend_fontsize=12, 
                                 figsize=(14,10),
                                 save=False)


# In[15]:


# proportions of spliced vs unspliced
fig_kwargs = {'figsize': (10, 10)}
if save_figure: fig_kwargs['save'] = "proportions_of_spliced_vs_unspliced.pdf"
scv.pl.proportions(adata, **fig_kwargs)


# ## CellRank Meets RNA Velocity

# ### Compute Velocity Kernel and transition matrix

# In[16]:


# Compute Velocity Kernel
vk = cr.kernels.VelocityKernel(adata, vkey='velocity')
# Compute Transition Matrix
vk.compute_transition_matrix()
# Underlying base kernels.
print(vk.kernels)
# shape of matrix
print(vk.shape)


# In[17]:


vk.plot_projection(basis='umap', recompute=True)


# In[18]:


adata


# ### Compute Connectivity Kernel and transition matrix

# In[19]:


# Compute similarity-based Connectivity Kernel
ck = cr.kernels.ConnectivityKernel(adata)
ck.compute_transition_matrix()
# Underlying base kernels.
print(ck.kernels)
# shape of matrix
print(ck.shape)


# In[20]:


ck.plot_projection(basis='umap', recompute=True)


# ### Combined Kernel for further Analysis (Estimator)
# 
# RNA velocity can be a very noise quantity; to make our computations more robust, we combine the VelocityKernel with the similarity-based ConnectivityKernel.

# In[21]:


# Combine Velocity Kernel with gene expression similarity kernel
combined_kernel = 0.9 * vk + 0.1 * ck
print(combined_kernel)


# In[22]:


combined_kernel.plot_projection(basis='umap', recompute=True)


# In[24]:


combined_kernel.write_to_adata(key=None, copy=False)


# ### Transition matrix 

# In[25]:


# Get the transition matrix
T = adata.obsp['T_fwd'].A


# ## Computing Initial and Terminal States (GPCCA Estimator)

# ## Advanced Estimator usage 

# In[26]:


# Generalized Perron Cluster Cluster Analysis (GPCCA) Estimator
g = cr.estimators.GPCCA(combined_kernel)
print(g)


# In[27]:


# Compute the Schur decomposition
g.compute_schur()
# Plot the top eigenvalues in a real or a complex plane.
g.plot_spectrum(real_only=True)
# real_only (Optional[bool]) – Whether to plot only the real part of the spectrum. If None, plot real spectrum if no complex eigenvalues are present.
g.plot_spectrum(real_only=False)


# In[28]:


# Plot the Schur matrix
g.plot_schur_matrix(title='schur matrix', cmap='viridis', figsize=None, dpi=80, save=None)


# In[29]:


g.schur_vectors


# ### Compute fixed Macrostates

# Using the graph of eigenvalues, we decide to pick 4 macrostates.

# In[32]:


# compute macro states specify n_states to compute 
g.compute_macrostates(n_states=4, n_cells=400, cluster_key="clusters")


# In[33]:


# plot macro states
# GPCCA.plot_macrostates(which, states=None, color=None, discrete=True,mode=PlotMode.EMBEDDING, time_key='latent_time', same_plot=True, title=None, cmap='viridis', **kwargs)

g.plot_macrostates(which="all",discrete=True, legend_loc="right", s=100)
g.plot_macrostates(which="all",discrete=False, legend_loc="right", save=f'{out_dir}/4_macrostates_new.png')


# #### Macrostates Composition
# Now we have macrostates. We assign one label per macrostate based on the underlying 'clusters' annotation. However, that does not imply that all cells within one macrostate are from the same underlying cluster as we use majority voting. We can visualize the relationship between clusters and macrostates. We show below the distribution over cluster membership for each of the cells most confidently assigned to a certain macrostate.
# 
# With some exceptions, most macrostates recruit their top-cells from a single underlying cluster. This plotting function works for any observation-level covariate

# In[34]:


# Plot macrostates composition frequency depends on n_cells
g.plot_macrostate_composition(key="clusters", figsize=(10, 5))


# #### Coarse-grained transition matrix
# 
# To get an idea of how these macrostates are related, we plot the coarse-grained transition matrix. This transition matrix aggregates the single-cell Markov chain to a macrostate-level Markov chain. Entries in this matrix tell us how likely macrostates are to transition into one-another. We identify initial and terminal states based on the following criteria: 
# - *terminal_states* are very stable (large value on the diagonal). They can have incoming probability mass, but almost no outgoing probability mass.
# - *initial_states* are less stable (smaller values on the diagonal and in the coarse-grained stationary distribution). They can have outgoing, but almost no incoming probability mass.
# - *intermediate_states* are just all remaining macrostates which we classified neither as terminal nor as initial.
# 
# By default, macrostates are ordered according to their diagonal value, increasing from left to right. The diagonal value is a proxy for a states’ metastability, i.e. the probability of leaving the state once entered. (check the corresponding matrix element).ctly. 

# In[35]:


# Plot Coarse-grained transition matrix
g.plot_coarse_T(show_stationary_dist=True, 
                show_initial_dist=True, 
                order='stat_dist', 
                cmap='viridis',
                xtick_rotation=45, 
                annotate=True, 
                show_cbar=True, 
                title=None, 
                figsize=(8, 8), 
                dpi=80, 
                save=None)


# In[36]:


g.coarse_stationary_distribution


# ### Compute terminal states using stability index threshold

# In[37]:


g.predict_terminal_states(allow_overlap=True)


# In[40]:


g.plot_macrostates(which="terminal", discrete=True, legend_loc="right", s=100) 
g.plot_macrostates(which="terminal", discrete=False, legend_loc="right")


# In[50]:


#variable definition
terminal_states=['0', '3', '7', '8']


# ### Compute Initial macrostate

# Force initial state 8, according to higher initial dist in coarse-grained transition matrix

# In[42]:


initial_states=["8"]
g.set_initial_states(states=initial_states,
                       n_cells=200, 
                       allow_overlap=True
                      )


# In[43]:


# Plot initial macrostate
g.plot_macrostates(which="initial", discrete=True, legend_loc="right", s=100)
g.plot_macrostates(which="initial", discrete=False, legend_loc="right")


# ### EDA on Macrostates

# In[44]:


g.plot_macrostates(which="all",
                   discrete=False,
                   color=None,
                   same_plot=False,
                   time_key='velocity_pseudotime',
                   mode='time',
                   legend_loc=None)


# In[45]:


g.plot_macrostates(which="all",
                   discrete=False,
                   color=None,
                   same_plot=False,
                   time_key='latent_time',
                   mode='time',
                   legend_loc="right")


# In[46]:


# macrostates_memberships is Soft assignment of microstates (cells) to macrostates
g.macrostates_memberships


# In[47]:


# plot cell count and probability for those cells belonging to each macrostate
df = pd.DataFrame(g.macrostates_memberships)
df.columns = g.macrostates.cat.categories.tolist()

# Create a histogram using Plotly Express
fig = px.histogram(df)
fig.update_layout(title='Proability of number of cells belonging to each macrostate',
                  yaxis_title='Cell Count', 
                  xaxis_title='Probability', 
                  legend_title_text='Macrostates', 
                  height=600, width=800)

# Show the figure
fig.show()


# In[51]:


# Plot histogram distribution of Initial, terminal or macro states
macrostate_names = terminal_states
# Calculate the number of macrostates and the number of rows and columns for subplots
num_macrostates = len(macrostate_names)
num_cols = 2  # Number of columns for subplots
num_rows = (num_macrostates + 1) // num_cols  # Calculate the number of rows

# Create a figure with a grid of subplots
fig, axes = pl.subplots(num_rows, num_cols, figsize=(10, 8))
fig.suptitle('Distribution of Macrostates')

# Flatten the axes array for easy iteration
axes = axes.flatten()

# Create histograms for each macrostate in a loop
for i, macrostate_name in enumerate(macrostate_names):
    result = pd.DataFrame(g.terminal_states_memberships[macrostate_name])
    
    # Plot the histogram in the appropriate subplot
    axes[i].hist(result, bins=100)
    axes[i].set_xlabel('probability')
    axes[i].set_ylabel('cell count')
    axes[i].set_title(f'Distribution of {macrostate_name} memberships probabilities')

# Hide any empty subplots
for i in range(num_macrostates, num_cols * num_rows):
    fig.delaxes(axes[i])

# Adjust subplot spacing
pl.tight_layout()
pl.show()


# ## Compute fate probabilities

# We compute fate_probabilities by aggregating over all random walks that start in a given cell and end in some terminal population. 
# [compute_fate_probabilities](https://cellrank.readthedocs.io/en/latest/api/_autosummary/estimators/cellrank.estimators.GPCCA.html#cellrank.estimators.GPCCA.compute_fate_probabilities)
# 
# For each cell, this computes the probability of being absorbed in any of the terminal_states. In particular, this corresponds to the probability that a random walk initialized in transient cell 
#  will reach any cell from a fixed transient state before reaching a cell from any other transient state.

# In[52]:


g.compute_fate_probabilities(solver='gmres', 
                             use_petsc=True, 
                             n_jobs=None, 
                             backend='loky', 
                             show_progress_bar=True, 
                             tol=1e-06, 
                             preconditioner=None)


fig_kwargs = {'figsize': (15, 10)}
if save_figure: fig_kwargs['save'] = f"{out_dir}/umap_terminal_macrostates_fate_probabilities.png"
g.plot_fate_probabilities(states=None, 
                          color=None, 
                          mode='embedding', 
                          time_key=None, 
                          same_plot=False, 
                          title=None,
                          dpi= 200,
                          cmap='viridis',
                          **fig_kwargs)


# We can visualize these probabilities in compact form by coloring cells according to their most-likely fate, with color intensity reflecting the degree of lineage-bias by using plot_fate_probabilities().

# In[53]:


# Fate probabilities in terms of pseudotime
for state in terminal_states:
    g.plot_fate_probabilities(states=[state], 
                              color=None, 
                              mode='time', 
                              time_key="latent_time", 
                              same_plot=False, 
                              title=None, 
                              legend_loc='right', 
                              cmap='viridis',
                              figsize= (15, 10),
                              dpi=100,
                              save = f"{out_dir}/fate_probabilities_velocity_pseudotime_plot_terminal_macrostates_{state}.pdf")


# In[54]:


# plotly histogram for fate probabilities for terminal states/lineages
df = pd.DataFrame(g.fate_probabilities)
df.columns = terminal_states

# Create a histogram using Plotly Express
fig = px.histogram(df)
fig.update_layout(title='histogram for fate probabilities for terminal states/lineages',
                  yaxis_title='Cell Count', 
                  xaxis_title='Probability', 
                  legend_title_text='Macrostates', 
                  height=600, width=800)

# Show the figure
fig.show()


# ## Compute lineage driver Genes 

# To infer putative driver genes for any of these trajectories, we correlate expression values with fate probabilities.
# 
# compute_lineage_drivers:
# Correlates gene expression with lineage probabilities, for a given lineage and set of clusters. Often, it makes sense to restrict this to a set of clusters which are relevant for the specified lineages.

# In[55]:


# Compute lineage driver genes
g.compute_lineage_drivers(['3' , '7'],
                          method='fisher',
                          cluster_key='clusters', 
                          clusters=None,
                          layer=None, 
                          use_raw=False, 
                          confidence_level=0.95, 
                          n_perms=1000, 
                          seed=None,
                          show_progress_bar=True
                          )


# In[60]:


#save result in a csv file 
result = g.compute_lineage_drivers(
    ['3', '7'],
    method='fisher',
    cluster_key='clusters', 
    clusters=None,
    layer=None, 
    use_raw=False, 
    confidence_level=0.95, 
    n_perms=1000, 
    seed=None,
    show_progress_bar=True
)

df = pd.DataFrame(result)

df.to_csv(f'{out_dir}/lineage_driver_genes.csv', index=True)


# #### Plot driver genes on UMAP 

# In[58]:


# Plot top driver genes in terms of correlation value (corr)  
g.plot_lineage_drivers(lineage= '3', 
                           n_genes=10, 
                           use_raw=False, 
                           ascending=False, 
                           ncols=2, 
                           title_fmt='{gene} corr={corr:.4e} qval={qval:.4e}', 
                           figsize=None, 
                           dpi=100,
                           fontsize=14, 
                           save=f'{out_dir}/top_driver_genes_immature_ iMSNs.pdf')


# In[59]:


# Plot top driver genes in terms of correlation value (corr) 
g.plot_lineage_drivers(lineage= '7', 
                           n_genes=10, 
                           use_raw=False, 
                           ascending=False, 
                           ncols=2, 
                           title_fmt='{gene} corr={corr:.4e} qval={qval:.4e}', 
                           figsize=None, 
                           dpi=100, 
                           fontsize=14, 
                           save=f'{out_dir}/top_driver_genes_dMSNs.pdf')


# ## Save Object

# In[ ]:


adata


# In[ ]:


adata.obsm['lineages_fwd']


# In[ ]:


#to save the kernel to the disk for further analysis.h5ad file 
# adata.write_h5ad('/gpfs/projects/bsc83/Projects/Creatio/integration_mouse/cellfate_analysis/v15/adata_cellrank_wt_v4.h5ad')


# In[ ]:


adata_cr_wt= sc.read("/gpfs/projects/bsc83/Projects/Creatio/integration_mouse/cellfate_analysis/v15/adata_cellrank_wt_v4.h5ad")


# In[ ]:


adata_cr_wt


# In[ ]:


adata_cr_wt.obsm['lineages_fwd']


# In[ ]:


type(adata.obsm['lineages_fwd']) # before saving object


# In[ ]:


type(adata_cr_wt.obsm['lineages_fwd']) # in saved object


# In[ ]:


file_size_bytes = os.path.getsize('/gpfs/projects/bsc83/Projects/Creatio/integration_mouse/cellfate_analysis/v15/adata_cellrank_wt_v4.h5ad')

# Convert bytes to megabytes
file_size_mb = file_size_bytes / (1024 ** 2)

print(f"The size of the file is approximately {file_size_mb:.2f} MB")

