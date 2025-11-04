"""Spatial count operator

Performs a spatial count using neighbourhood analysis. It counts the frequency 
of neighbour cluster per cell using either radius or knn method.

Uses the scimap library's spatial_count function.
"""

import numpy as np
import pandas as pd
import anndata as ad
import scimap as sm
import scimap.tl
import scimap.pl
from tercen.client import context as ctx

# For local testing, uncomment and provide your credentials:
# tercenCtx = ctx.TercenContext(
#     workflowId="YOUR_WORKFLOW_ID",
#     stepId="YOUR_STEP_ID",
#     serviceUri="https://tercen.com/api/v1",
#     authToken="YOUR_TOKEN_HERE"
# )

tercenCtx = ctx.TercenContext()

# Get operator properties
method = tercenCtx.operator_property('method', typeFn=str, default='radius')
radius = tercenCtx.operator_property('radius', typeFn=float, default=30.0)
knn = tercenCtx.operator_property('knn', typeFn=int, default=10)
n_clusters = tercenCtx.operator_property('n_clusters', typeFn=int, default=8)

# Get all available data from Tercen
df = tercenCtx.select(df_lib="pandas")

# Map columns - look for cluster/phenotype column and label column
# Common patterns: ds3.cluster_id, cluster_id, phenotype, colors
cluster_col = None
for col in df.columns:
    if 'cluster' in col.lower() or 'phenotype' in col.lower() or col == 'colors':
        cluster_col = col
        break

if cluster_col is None:
    raise RuntimeError(f"Could not find cluster/phenotype column. Available columns: {df.columns.tolist()}")

# Look for label column (cell_id)
# Common patterns: Object, label, cell_id, id
label_col = None
for col in df.columns:
    if col == 'Object' or 'label' in col.lower() or col == 'cell_id' or col == 'id':
        label_col = col
        break

if label_col is None:
    raise RuntimeError(f"Could not find label column. Available columns: {df.columns.tolist()}")

# Keep the actual column names from the data for output relation
# Just create a 'colors' alias for easier processing in the loop
df['colors'] = df[cluster_col]

# Remove duplicate rows - keep only one row per unique cell (Object)
# This handles cases where the same cell appears multiple times in the input
df_unique = df.drop_duplicates(subset=[label_col, '.ci', '.ri'], keep='first')

# Prepare data for scimap: need X_centroid, Y_centroid, and imageid
df_unique = df_unique.rename(columns={'.x': 'X_centroid', '.y': 'Y_centroid'})
df_unique['imageid'] = df_unique['.ci'].astype(str) + '_' + df_unique['.ri'].astype(str)

# Create AnnData object for scimap
# Extract marker data (if any numeric columns exist, otherwise use dummy data)
numeric_cols = df_unique.select_dtypes(include=[np.number]).columns.tolist()
# Remove coordinate and index columns
exclude_cols = ['X_centroid', 'Y_centroid', '.ci', '.ri', label_col]
marker_cols = [col for col in numeric_cols if col not in exclude_cols]

if len(marker_cols) > 0:
    data_matrix = df_unique[marker_cols].values
else:
    # Create dummy data if no markers available
    data_matrix = np.zeros((len(df_unique), 1))

# Create metadata dataframe
obs_data = df_unique[[label_col, 'colors', 'imageid', 'X_centroid', 'Y_centroid', '.ci', '.ri']].copy()
obs_data.index = obs_data.index.astype(str)

# Create AnnData object
adata = ad.AnnData(X=data_matrix, obs=obs_data)

# Run scimap spatial_count
adata = sm.tl.spatial_count(
    adata, 
    phenotype='colors',
    method=method,
    radius=radius if method == 'radius' else None,
    knn=int(knn) if method == 'knn' else None,
    imageid='imageid',
    x_coordinate='X_centroid',
    y_coordinate='Y_centroid',
    label='spatial_count'
)

# Cluster the spatial_count results using k-means to identify neighbourhood regions
adata = sm.tl.spatial_cluster(
    adata, 
    df_name='spatial_count', 
    method='kmeans', 
    k=n_clusters, 
    label='neighbourhood_cluster'
)

# Extract the neighbourhood cluster assignments from adata.obs
result_df = adata.obs[[label_col, '.ci', '.ri', 'neighbourhood_cluster']].copy()

# Ensure .ci and .ri are integers as required by Tercen
result_df['.ci'] = result_df['.ci'].astype(np.int32)
result_df['.ri'] = result_df['.ri'].astype(np.int32)

# Add namespace and save
result_df = tercenCtx.add_namespace(result_df)
tercenCtx.save(result_df)

