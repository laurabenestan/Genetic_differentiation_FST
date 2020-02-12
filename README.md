# Genetic_differentiation_FST

Calculate Weir and Cockram FST using assigner package available in R.
To download and install **assigner** package see [Thierry Gosselin github] (https://github.com/thierrygosselin/assigner)

Then, in your R environment prepare the workspace.
### Remove last features
``` r
rm(list=ls())
ls()
```

### Download libraries
``` r
library(assigner)
library(ggplot2)
library(reshape)
library(ggdendro)
library(gridExtra)
library(radiator)
library(ape)
library(poppr)
```

#### Import your vcf on a tidy format
``` r
tidy_data <- tidy_genomic_data(data = "test.vcf",strata="test_n.txt")
```

#### Calculate fst on this tidy data using Confidenec Intervals of 97.5% and 2.5%
``` r
fst_markers <- fst_WC84(
  data = tidy_data,
  holdout.samples = NULL, 
  pairwise = TRUE,
  ci = TRUE,
  iteration.ci = 100,
  quantiles.ci = c(0.025, 0.975),
  digits = 9,
  parallel.core = 8
)
```

### Check what is present in the fst_marker object
``` r
names(fst_markers)
```

### Save the dataframe recording each pairwise FST and CI intervals 
``` r
write.table(as.data.frame(fst_markers$pairwise.fst), "Fst_matrix.txt",quote=FALSE, row.names=FALSE, sep="\t") 
head(fst_markers$pairwise.fst)
```

### Explore the Fst values
``` r
summary(fst_markers$pairwise.fst$FST)
fst <- fst_markers$pairwise.fst.full.matrix
class(fst)
```

### Save the results into a data frame
``` r
fst_df <- as.data.frame(fst_markers$fst.markers)
```

### Add the markers information in a cloumn named CHROM
``` r
fst_df$CHROM <- seq(1,24603,1)
```

### Do a graph on the distribution of Fst values
``` r
g1 <- ggplot(data=fst_df, aes(fst_df$FST)) + 
  geom_histogram(breaks=seq(0, 1, by=0.01), 
                 col="black", 
                 fill="black", 
                 alpha=.2)+
  xlab("FST values across populations")+
  ylab("Number of SNPs")+
  theme_classic()
```

### Save the Fst graph
``` r
ggsave("Fst_graph.pdf", width=5, height=5)
```

### Save results
``` r
write.table(as.data.frame(fst_df), "Fst_markers.txt",quote=FALSE, row.names=FALSE, sep="\t") 
```

### Check the summary of the FST values 
``` r
summary(fst_df)
```

### Check the quantile
``` r
quantile(fst_df$FST, c(.01,.99))
```

### Subset the loci related to the quantile values
``` r
outliers_99 <- subset(fst_df, subset=fst_df$FST>0.827055)
```

### Save results
``` r
write.table(as.data.frame(outliers_99), "Outliers_99.txt",quote=FALSE, row.names=FALSE, sep="\t") 
save.image(file = "Fst_sebastes.RData")
```

### Create a dendrogram
``` r
hc= hclust(dist(fst))
ggdendrogram(hc, theme_dendro=F, color= "black", size=6,rotate = TRUE)+
  theme(axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
```

### Save the dendrogram
``` r
ggsave("Dendrogram.pdf", width=5, height=5)
```

