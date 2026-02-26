library(caret)
library(quantregForest)

qrf_caret <- list(
  
  type = "Regression",
  library = "quantregForest",
  
  parameters = data.frame(
    parameter = c("mtry", "nodesize"),
    class = c("numeric", "numeric"),
    label = c("mtry", "Node size")
  ),
  
  grid = function(x, y, len = NULL, search = "grid") {
    expand.grid(
      mtry = unique(pmax(1, floor(seq(1, ncol(x), length.out = len)))),
      nodesize = c(5, 10, 15)
    )
  },
  
  fit = function(x, y, wts, param, lev, last, classProbs, ...) {
    quantregForest(
      x = x,
      y = y,
      mtry = param$mtry,
      nodesize = param$nodesize,
      keep.inbag = TRUE,
      ...
    )
  },
  
  predict = function(modelFit, newdata, submodels = NULL) {
    predict(modelFit, newdata, what = c(0.5, 0.05, 0.95))
  },
  prob = NULL, 
  
  sort = function(x) x[order(x$mtry), ],
  
  levels = function(x) NULL
)
