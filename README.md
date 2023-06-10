## J-score

J-score measures similarities between two clustering structures. It is implemented as an R/CRAN library. https://cran.r-project.org/web/packages/jScore/index.html

## Install 
````
install.packages('jScore')
````
## Usage 
````
library(jScore)

## an example comparing two clustering assignments
truth=c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3)  # A numeric vector of truth class labels.
pred= c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5)  # A numeric vector of predicted class labels.
jscore(truth, pred)  ## 0.5414634
````
## Reference
Please cite our publication in arXiv.
Ahmadinejad N, Liu L (2021). J-Score: A Robust Measure of Clustering Accuracy https://arxiv.org/abs/2109.01306
(Manuscript under review in PeerJ)
