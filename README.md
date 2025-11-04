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
`radius`   | number,setting for the size of radius to use, default 30
`knn`      | integer, setting for the number of centers,  default 10


Output|.
---|---
`Object` (or label column name) | integer, the cell id from input (creates output relation)
`name_neigbours`  | text, name of neighbourhood cluster phenotype
`count_neighbours`| integer, count of that neighbour phenotype for this cell

### Details

Performs spatial count analysis using scikit-learn's NearestNeighbors algorithm:
- Implements the spatial_count concept from scimap but using scikit-learn
- Uses spatial coordinates (X, Y) from single-cell data to perform neighborhood analysis
- Supports two methods: radius-based or k-nearest neighbors (knn)
- Automatically deduplicates input rows to ensure each cell is counted once
- Preserves the original label column name (e.g., `Object`) for proper output relation

### Example Output

For each cell, the operator returns one row per neighbor phenotype found:

```
Object, name_neigbours, count_neighbours
9.0, c00, 4
9.0, c15, 3
9.0, c14, 3
48.0, c15, 4
48.0, c11, 3
```

This shows that cell 9.0 has 4 neighbors of phenotype c00, 3 of c15, etc.
