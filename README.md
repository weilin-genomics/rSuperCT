rSuperCT is the R version of SuperCT, implement the supervised-learning framework for enhanced characterization of single-cell transcriptomic profiles

## Installation
```
devtools::install_github('weilin-genomics/rSuperCT')
library(rSuperCT)
```
## Usage of rSuperCT
Follow below steps to predict cell types for scRNA-seq data. For now, only human and mouse are supported.

1) (Skip this step if myces is a matrix) Take pbmc_small dataset in Seurat package as example, the input format for ImportData() function can be Seurat object, (sparse) matrix or data frame in which the rows represent for features and the columns for cells.
```{r}
library(Seurat)
myces <- ImportData(pbmc_small)
myces # show the number of cells and features
```
2) To predict cell types using pre-trained and well-named models published by https://github.com/weilin-genomics/SuperCT,
specify the project name by model parameter and corresponding compressed model file will be downloaded and uncompressed to local 'models' directory from 'https://github.com/weilin-genomics/rSuperCT_models'. Myces is a matrix with cells as columns and features as rows.
As prediction done, A column pred_types was saved in meta.data slot of your CellESet object. 'generic_38CellTypes' is
```{r}
dir.create('./models', showWarnings = FALSE)
myces <- PredCellTypes(myces, species = 'human', model = 'generic_38CellTypes', results.dir = './models')
```
3) Visualize the predictive results for summarization.
```{r}
plotHist(myces)
```
## References
[1] Xie Peng and Gao Mingxuan (2019). SuperCT: a supervised-learning framework for enhanced characterization of single-cell transcriptomic profiles. https://doi.org/10.1093/nar/gkz116. Nucleic Acids Research.

[2] https://github.com/weilin-genomics/SuperCT
## Contact
If you have any suggestion, questions and bugs report, feel free to contact weilin.baylor@gmail.com.
