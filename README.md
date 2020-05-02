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
2) To predict cell types using pre-trained models, you can find the models from 'https://github.com/weilin-genomics/rSuperCT_models' that can be downloaded and uncompressed to the local 'models' directory. As prediction done, A column pred_types was saved in meta.data slot of your CellESet object. 'generic_38CellTypes' is the model you used for predictions.
```{r}
dir.create('./models', showWarnings = FALSE)
###download and uncompress the model to ./models directory##
myces <- PredCellTypes(myces, species = 'human', model = 'generic_38celltypes', results.dir = './models')
####you can find the predictions under myces@meta.data$pred_types
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
