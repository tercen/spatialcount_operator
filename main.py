"""Spatial count operator

Performs a spatial count using neighbourhood analysis. It counts the frequency 
of neighbour cluster per cell using either radius or knn method.

Uses the scimap library's spatial_count function.
"""

import numpy as np
import pandas as pd
import anndata as ad
import scimap as sm
from tercen.client import context as ctx

# For local testing with live Tercen connection, uncomment and provide credentials:
# tercenCtx = ctx.TercenContext(
#     workflowId="YOUR_WORKFLOW_ID",
#     stepId="YOUR_STEP_ID",
#     serviceUri="https://tercen.com/api/v1",
#     authToken="YOUR_TOKEN_HERE"
# )

tercenCtx = ctx.TercenContext()

# Get operator properties
method = tercenCtx.operator_property('method', typeFn=str, default='radius')
radius = tercenCtx.operator_property('radius', typeFn=float, default=80.0)
knn = tercenCtx.operator_property('knn', typeFn=int, default=10)
n_clusters = tercenCtx.operator_property('n_clusters', typeFn=int, default=8)

# Get all data including colors
df = tercenCtx.select(df_lib="pandas")

# Get column factors in order: [imageid, cellid, y_coordinate, x_coordinate]
# Based on Tercen projection order: filename, Object, centroid-0, centroid-1
col_df = tercenCtx.cselect(df_lib="pandas")
col_df['.ci'] = col_df.index

# Get column names in order (excluding .ci which we added)
col_names = [c for c in col_df.columns if c != '.ci']

if len(col_names) >= 4:
    # Use positional mapping based on order
    # Order: filename, Object, centroid-0 (X), centroid-1 (Y)
    imageid_col = col_names[0]  # filename
    cellid_col = col_names[1]   # Object
    x_col = col_names[2]         # centroid-0 (X coordinate)
    y_col = col_names[3]         # centroid-1 (Y coordinate)
else:
    raise RuntimeError(f"Expected at least 4 column factors. Found: {col_names}")

# Merge column data
df = df.merge(col_df, on='.ci', how='left')

# Get cluster/phenotype info from colors projection
# tercenCtx.colors returns the column name(s) used for colors
color_col_names = tercenCtx.colors
if isinstance(color_col_names, list) and len(color_col_names) > 0:
    color_col = color_col_names[0]
elif isinstance(color_col_names, str):
    color_col = color_col_names
else:
    raise RuntimeError(f"Could not determine color column from tercenCtx.colors: {color_col_names}")

# The color data is already in the dataframe
if color_col in df.columns:
    df['colors'] = df[color_col].astype(str)
else:
    raise RuntimeError(f"Color column '{color_col}' not found in dataframe. Available: {df.columns.tolist()}")

# Prepare data for scimap
df['X_centroid'] = df[x_col]
df['Y_centroid'] = df[y_col]
df['imageid'] = df[imageid_col].astype(str)
df['cell_id'] = df[cellid_col]

# Remove duplicate rows
df_unique = df.drop_duplicates(subset=['.ci', '.ri'], keep='first')

# Create AnnData object for scimap
# Extract marker expression data from the main table
# Get numeric columns that are likely markers (exclude coordinates, indices, labels)
numeric_cols = df_unique.select_dtypes(include=[np.number]).columns.tolist()
exclude_cols = ['X_centroid', 'Y_centroid', '.ci', '.ri', 'cell_id', x_col, y_col]
marker_cols = [col for col in numeric_cols if col not in exclude_cols and not col.startswith('.')]

if len(marker_cols) > 0:
    data_matrix = df_unique[marker_cols].values
else:
    # Create dummy data if no markers available
    data_matrix = np.zeros((len(df_unique), 1))

# Create metadata dataframe
obs_data = df_unique[['colors', 'imageid', 'X_centroid', 'Y_centroid', '.ci', '.ri']].copy()
obs_data.index = obs_data.index.astype(str)
obs_data['cell_id'] = obs_data['.ci']

# Ensure colors is string type
obs_data['colors'] = obs_data['colors'].astype(str)

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

# Check if spatial_count has any data
if adata.uns['spatial_count'].shape[1] > 0:
    # Cluster the spatial_count results using k-means to identify neighbourhood regions
    adata = sm.tl.spatial_cluster(
        adata, 
        df_name='spatial_count', 
        method='kmeans', 
        k=n_clusters, 
        label='neighbourhood_cluster'
    )
    # Extract the neighbourhood cluster assignments from adata.obs
    result_df = adata.obs[['.ci', '.ri', 'neighbourhood_cluster']].copy()
else:
    # Return just the cell identifiers without clustering
    result_df = obs_data[['.ci', '.ri']].copy()
    result_df['neighbourhood_cluster'] = 0

# Ensure .ci and .ri are integers as required by Tercen
result_df['.ci'] = result_df['.ci'].astype(np.int32)
result_df['.ri'] = result_df['.ri'].astype(np.int32)

# Add namespace and save
result_df = tercenCtx.add_namespace(result_df)

# Save to Tercen
tercenCtx.save(result_df)

