"""Spatial count operator

Performs a spatial count using neighbourhood analysis. It counts the frequency 
of neighbour cluster per cell using either radius or knn method.

Uses the scimap library's spatial_count function.
"""

import numpy as np
import pandas as pd
import anndata as ad
import scimap as sm
import gc
from tercen.client import context as ctx

# Initialize Tercen context
try:
    tercenCtx = ctx.TercenContext()
except Exception as e:
    raise RuntimeError(f"Failed to initialize Tercen context: {e}")

# Get operator properties
method = tercenCtx.operator_property('method', typeFn=str, default='radius')
radius = tercenCtx.operator_property('radius', typeFn=float, default=80.0)
knn = tercenCtx.operator_property('knn', typeFn=int, default=10)
n_clusters = tercenCtx.operator_property('n_clusters', typeFn=int, default=8)

# Validate parameters
if method not in ['radius', 'knn']:
    raise ValueError(f"method must be 'radius' or 'knn', got '{method}'")
if radius <= 0:
    raise ValueError(f"radius must be positive, got {radius}")
if knn < 1:
    raise ValueError(f"knn must be >= 1, got {knn}")
if n_clusters < 2:
    raise ValueError(f"n_clusters must be >= 2, got {n_clusters}")

# Get only essential columns to minimize memory
df = tercenCtx.select(['.ci', '.ri'], df_lib="pandas")

# Get column factors
col_df = tercenCtx.cselect(df_lib="pandas")
col_names = col_df.columns.tolist()

if len(col_names) >= 4:
    # Order: imageid, cellid, x_coord, y_coord
    imageid_col, cellid_col, x_col, y_col = col_names[0], col_names[1], col_names[2], col_names[3]
else:
    raise RuntimeError(f"Expected at least 4 column factors. Found: {col_names}")

# Get color column name
color_col_names = tercenCtx.colors
color_col = color_col_names[0] if isinstance(color_col_names, list) else color_col_names

# Select only needed columns from main data
needed_cols = ['.ci', '.ri', color_col]
df = tercenCtx.select(needed_cols, df_lib="pandas")

# Add column factors directly (avoid merge for memory efficiency)
col_df_reset = col_df.reset_index()
col_df_reset['.ci'] = col_df_reset.index

# Merge only needed column data
df = df.merge(col_df_reset[['.ci', imageid_col, cellid_col, x_col, y_col]], on='.ci', how='left')

# Rename for scimap
df.rename(columns={
    color_col: 'colors',
    x_col: 'X_centroid',
    y_col: 'Y_centroid',
    imageid_col: 'imageid'
}, inplace=True)

# Convert types in place
df['colors'] = df['colors'].astype(str)
df['imageid'] = df['imageid'].astype(str)

# Remove duplicates
df.drop_duplicates(subset=['.ci', '.ri'], keep='first', inplace=True)

# Process each image separately to reduce memory usage
results = []
unique_images = df['imageid'].unique()

for img in unique_images:
    # Filter data for this image only
    img_df = df[df['imageid'] == img].copy()
    
    # Create minimal AnnData object for this image
    data_matrix = np.zeros((len(img_df), 1), dtype=np.float32)
    obs_data = img_df[['colors', 'imageid', 'X_centroid', 'Y_centroid', '.ci', '.ri']].copy()
    obs_data.index = obs_data.index.astype(str)
    
    adata = ad.AnnData(X=data_matrix, obs=obs_data)
    
    # Run spatial_count for this image
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
    
    # Check if spatial_count has data
    if adata.uns['spatial_count'].shape[1] > 0:
        # Cluster for this image
        adata = sm.tl.spatial_cluster(
            adata, 
            df_name='spatial_count', 
            method='kmeans', 
            k=n_clusters, 
            label='neighbourhood_cluster'
        )
        img_result = adata.obs[['.ci', '.ri', 'neighbourhood_cluster']].copy()
    else:
        img_result = obs_data[['.ci', '.ri']].copy()
        img_result['neighbourhood_cluster'] = 0
    
    results.append(img_result)
    
    # Clean up to free memory
    del adata, data_matrix, obs_data, img_df, img_result
    gc.collect()

# Combine results from all images
result_df = pd.concat(results, ignore_index=True)
del results
gc.collect()

# Ensure .ci and .ri are integers as required by Tercen
result_df['.ci'] = result_df['.ci'].astype(np.int32)
result_df['.ri'] = result_df['.ri'].astype(np.int32)

# Add namespace and save
result_df = tercenCtx.add_namespace(result_df)

# Save to Tercen
tercenCtx.save(result_df)

