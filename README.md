# Spatial count Operator

### Description

Performs a spatial count using neighbourhood analysis. It counts the frequency of neighbour cluster per cell using scikit-learn's NearestNeighbors algorithm.

### Usage

Input|.
---|---
`y-axis`        | number, marker expression values (matrix)
`column` (1st)  | text, image identifier (e.g., filename)
`column` (2nd)  | number, unique cell identifier (e.g., Object)
`column` (3rd)  | number, X centroid position (e.g., centroid-0)
`column` (4th)  | number, Y centroid position (e.g., centroid-1)
`colors`        | text, phenotypic cluster assignment
`row`           | text, marker names (optional)

Settings|.
---|---
`method`   | radius or knn, default radius
`radius`   | number, setting for the size of radius to use, default 80
`knn`      | integer, setting for the number of neighbors, default 10
`n_clusters` | integer, number of neighborhood clusters for k-means, default 8


Output|.
---|---
`neighbourhood_cluster` | integer, the neighborhood cluster assignment (0 to n_clusters-1)

### Details

Performs spatial neighborhood analysis using the scimap library:
1. Reads marker expression data and cell metadata from Tercen projection
2. Uses column factors in order: imageid, cellid, x_coordinate, y_coordinate
3. Gets phenotype/cluster information from the colors projection
4. Counts neighbor phenotypes within a radius or using k-nearest neighbors
5. Performs k-means clustering on the spatial count results to identify neighborhood regions
6. Returns one neighborhood cluster assignment per cell

The operator uses:
- `sm.tl.spatial_count()` - counts phenotype frequencies in neighborhoods
- `sm.tl.spatial_cluster()` - clusters cells based on their neighborhood composition

**Note:** The operator uses generic column mapping based on factor order, not hardcoded column names. This makes it flexible for different datasets.

### Example Output

Each cell gets one neighborhood cluster assignment:

```
.ci, .ri, neighbourhood_cluster
0, 0, 0
1, 0, 6
0, 1, 0
1, 1, 6
2, 0, 2
```

Cells with similar neighborhood compositions are assigned to the same cluster (0 to n_clusters-1).

The output maintains the crosstab structure with `.ci` and `.ri` indices for proper integration with Tercen workflows.
