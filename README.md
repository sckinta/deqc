# deqc

### normalize count matrix
```R
count2tpm(count, geneLength, pseudoCount=0, log=F, pseudoTPM=1)
count2fpkm(count, geneLength, pseudoCount=0, log=F, pseudoFPKM=1)
```


### plot PCA 
plot PCA using normalized count and sample annotation
```R
plot_pca(tpm, col_anno, color_col, shape_col=NULL, text_col=NULL, point_size=5, plot_file="")
```
