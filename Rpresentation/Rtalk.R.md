Group Meeting
========================================================
author: Jens Kleinjung
date: 06.07.2017
autosize: true


Overview
========================================================
1. Semi-supervised RNA-Seq analysis
2. Data Pipelines (with 'GNU make')
3. Single-cell RNA-Seq


1. RNA-Seq analysis of ENS at day E13
========================================================



ggplot2
========================================================
(Taken from the ggplot2 wbesite and Wikipedia)
ggplot2 is a plotting system for R,
which tries to take the good parts of base and lattice graphics
and none of the bad parts.

Created by Hadley Wickham in 2005, ggplot2 is an implementation of
Leland Wilkinson's Grammar of Graphics—a general scheme for data visualization
which breaks up graphs into semantic components such as scales and layers.

The 'R Graphics Cookbook' by Winston Chang provides a set of recipes
to solve common graphics problems.
***

```r
install.packages("ggplot2")
```
![plot of chunk unnamed-chunk-2](Rtalk.R-figure/unnamed-chunk-2-1.png)

Biplot
========================================================
![alt text](biplot_semisuperv.png)
***

```r
library(devtools)
install_github("vqv/ggbiplot")

g = qplot(df.u$xvar, df.u$yvar, colour = expr.gradient, size = expr.mean, xlab = u.axis.labs[1], ylab = u.axis.labs[2]);

g = g + scale_colour_gradient(low = "blue", high = "orange"); 

g = g + scale_size(range = c(1, 7));

g = g + geom_segment(data = df.v[gene.glial.nidx, ], aes(x = 0, y = 0, xend = xvar, yend = yvar, size = 1), ...
```


Heatmap
========================================================
![alt text](heatmap_key.png)
***
R code


Gene Bars
========================================================
![alt text](Gbars.Tubb3.cutoff.png)
***
R code


3D Profile
========================================================
![alt text](bifurk_3d_box.png)
***
R code


RNA-Seq Alignment Pipeline
========================================================
- Dependency resolution: .fastq.gz <- fastq <- bam <- counts <- table <- normalised table
- GNU Make
    * Macros
    * Rules with dependencies + commands
    * Pattern rules

- Demos
    1. C program
    2. Variable compilation




