# Spatial Count Operator Tests

## Test Data

### Input: `spatial_input.csv`
- **12 cells** across **2 images** (Image1, Image2)
- **5 markers**: CD3, CD4, CD8, CD20, CD68
- **3 phenotypes**: T-cell, B-cell, Macrophage
- **Spatial layout**: Cells are grouped by phenotype to create distinct neighborhoods

### Spatial Layout (Image1):
```
Cluster 1 (T-cells):     cells 1,2,6,9  at X~10-15, Y~20-25
Cluster 2 (Macrophages): cells 4,5     at X~50-55, Y~20-25
Cluster 3 (B-cells):     cells 3,7,8   at X~90-95, Y~20-25
```

### Expected Output: `expected_output.csv`
- Each cell gets a **neighbourhood_cluster** assignment (0, 1, or 2)
- Results are expanded to all `.ri` values (one row per cell per marker)
- Total: 12 cells × 5 markers = 60 rows

### Test Configuration: `test.json`
- **Columns**: filename (imageid), Object (cellid), centroid-0 (X), centroid-1 (Y)
- **Rows**: CD3, CD4, CD8, CD20, CD68 (markers)
- **Colors**: phenotype (cluster assignment)
- **Parameters**: 
  - method: radius
  - radius: 10.0 pixels
  - n_clusters: 3

## Running Tests

Tests are run automatically by Tercen's CI/CD pipeline when:
1. Code is pushed to the repository
2. A new release tag is created

The test validates that the operator:
- Correctly reads spatial coordinates from column factors
- Uses phenotype data from colors projection
- Computes neighborhood clusters based on spatial proximity
- Returns results in the correct Tercen format (.ci, .ri, neighbourhood_cluster)

## Notes

- The expected output is approximate since k-means clustering can vary slightly
- The test primarily validates data structure and format, not exact cluster assignments
- For production use, test with larger datasets in Tercen workflows
