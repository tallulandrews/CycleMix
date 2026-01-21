# CycleMix
### Gaussian Clustering for CellCycle Assignment in scRNAseq

This is an R package which assigns single cell RNASeq data to cell-cycle stages using annotated marker genes. Various sets of marker genes used in the literature are provided or users can provide their own, both positive and negative markers are supported. 

Each stage is assigned independently using a gaussian mixture-model and cells can be assigned to one, multiple or no cell-cycle stage.

## Installation:

```r
require("remotes")
install_github('tallulandrews/CycleMix')
```

## Utilization:

There are two core functionality for CycleMix. The first is assigning cells to cell-cycle phases:

```r
# Human data
phases <- classifyCells(seurat_object, HGeneSets$Tirosh)
# Mouse data
phases <- classifyCells(seurat_object, MGeneSets$Cyclone)
````

Second is to regress differences between various cell-cycle phases:
```r
# Remove differences between G1S and G2M cells
corrected_expression <- regressCyclePartial(GetAssayData(seurat_object, assay="RNA", layer="data"), phases, method="phase", phases=c("G1S", "G2M"), allow_negative=TRUE, type="norm")
```

Please see the vignette for step-by-step guide to use CycleMix.
