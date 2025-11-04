"""Spatial count operator

Performs a spatial count using neighbourhood analysis. It counts the frequency 
of neighbour cluster per cell using either radius or knn method.

Based on the scimap spatial_count concept but implemented with scikit-learn.
"""

import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors
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

# Get all unique phenotypes across the dataset
all_phenotypes = df_unique['colors'].unique()

# Compute spatial counts per image group (.ci, .ri represent column/row grouping in crosstab)
# Within each group, we analyze spatial neighbors for each cell
results = []

for (ci, ri), group in df_unique.groupby(['.ci', '.ri']):
    if len(group) < 2:
        continue
    
    coords = group[['.x', '.y']].to_numpy(dtype=float)
    
    if method == 'knn':
        n_neighbors = min(knn + 1, len(group))
        nn = NearestNeighbors(n_neighbors=n_neighbors, algorithm='auto').fit(coords)
        distances, indices = nn.kneighbors(coords)
        # Remove self (first neighbor)
        neighbor_idx = [inds[1:] for inds in indices]
    else:  # radius
        nn = NearestNeighbors(radius=radius, algorithm='auto').fit(coords)
        distances, indices = nn.radius_neighbors(coords)
        # Remove self from neighbors
        neighbor_idx = [inds[inds != i] for i, inds in enumerate(indices)]
    
    # Count phenotype frequencies for each cell
    # Reset index to use positional indexing within this group
    group_reset = group.reset_index(drop=True)
    for i, neigh in enumerate(neighbor_idx):
        # Use the actual label value from the input data
        cell_label = group_reset.iloc[i][label_col]
        
        # Get phenotypes of neighbors
        if len(neigh) == 0:
            neigh_phenotypes = pd.Series([], dtype=object)
        else:
            neigh_phenotypes = group_reset.iloc[neigh]['colors']
        
        counts = neigh_phenotypes.value_counts()
        
        # Add a row for each possible phenotype (with count 0 if not found)
        for phenotype in all_phenotypes:
            count = counts.get(phenotype, 0)
            results.append({
                '.ci': ci,
                '.ri': ri,
                label_col: cell_label,
                'name_neigbours': phenotype,
                'count_neighbours': float(count)
            })

# Create result dataframe
result_df = pd.DataFrame(results)

# Ensure .ci and .ri are integers as required by Tercen
result_df['.ci'] = result_df['.ci'].astype(np.int32)
result_df['.ri'] = result_df['.ri'].astype(np.int32)

# Add namespace and save
result_df = tercenCtx.add_namespace(result_df)
tercenCtx.save(result_df)

