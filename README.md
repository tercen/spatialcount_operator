# Spatial count Operator

### Description

Performs a spatial count using neighbourhood analysis. It counts the frequency of neighbour cluster per cell using scikit-learn's NearestNeighbors algorithm.

### Usage

Input|.
---|---
`x-axis`        | number, the x centroid position of each cell  
`y-axis`        | number, the y centroid position of each cell
`row`           | text, an image grouping 
`column`        | text, an image grouping
`colors`        | text, phenotypic cluster 
`label`         | integer, a unique cell id

Settings|.
---|---
`method`   | radius or knn, default radius
`radius`   | number, setting for the size of radius to use, default 30
`knn`      | integer, setting for the number of neighbors, default 10
`n_clusters` | integer, number of neighborhood clusters for k-means, default 8


Output|.
---|---
`neighbourhood_cluster` | integer, the neighborhood cluster assignment (0 to n_clusters-1)

### Details

Performs spatial neighborhood analysis using the scimap library:
1. Counts neighbor phenotypes within a radius or using k-nearest neighbors
2. Performs k-means clustering on the spatial count results to identify neighborhood regions
3. Returns one neighborhood cluster assignment per cell

The operator uses:
- `sm.tl.spatial_count()` - counts phenotype frequencies in neighborhoods
- `sm.tl.spatial_cluster()` - clusters cells based on their neighborhood composition

### Example Output

Each cell gets one neighborhood cluster assignment:

```
Object, neighbourhood_cluster
9.0, 4
48.0, 4
56.0, 4
75.0, 4
110.0, 4
```

Cells with similar neighborhood compositions are assigned to the same cluster (0 to n_clusters-1).
