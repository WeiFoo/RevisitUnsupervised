# Revisit Unsupervised Effort-Aware Just-in-Time Defect prediction


This package is modified basd on the code provoided by Yang et al. [`` Effort-aware just-in-time defect prediction: Simple unsupervised models could be better than supervised models.''](https://cs.nju.edu.cn/yangyibiao/) 



# Dependence

* R
* RWeka, car, usdm, effsize, ScottKnott packages are required. Please run the following commands.

```R
	install.packages("RWeka")
	install.packages("car")
	install.packages("usdm")
	install.packages("effsize")
	install.packages("ScottKnott")
```

# Run
```
 source("exeMain.r")
```
# Folder Structure


```
.
├── LICENSE
├── README.md
└── jit
    ├── input   ## input data
    │   ├── bugzilla.csv
    │   ├── columba.csv
    │   ├── jdt.csv
    │   ├── mozilla.csv
    │   ├── platform.csv
    │   └── postgres.csv
    ├── output  ## reuslts are saved here accordingly
    │   ├── cross-project
    │   ├── cross-validation
    │   └── cross-validation-timewise
    └── script_r
        ├── ReportResults.r
        ├── core.r
        ├── evaluate.r
        ├── exeMain.r  ## the start point of R code
        ├── learner.r
        ├── packages
        │   ├── EMImputation1.0.1.zip
        │   ├── RBFNetwork1.0.8.zip
        │   ├── citationKNN1.0.1.zip
        │   ├── conjunctiveRule1.0.4.zip
        │   ├── gaussianProcesses1.0.1.zip
        │   ├── grading1.0.1.zip
        │   ├── gridSearch1.0.7.zip
        │   ├── gridSearch1.0.8.zip
        │   ├── gridSearch1.0.9.zip
        │   ├── isotonicRegression1.0.1.zip
        │   ├── kernelLogisticRegression1.0.0.zip
        │   ├── lazyBayesianRules1.0.1.zip
        │   ├── leastMedSquared1.0.1.zip
        │   ├── linearForwardSelection1.0.1.zip
        │   ├── multiBoostAB1.0.1.zip
        │   ├── paceRegression1.0.1.zip
        │   ├── partialLeastSquares1.0.4.zip
        │   ├── racedIncrementalLogitBoost1.0.1.zip
        │   ├── ridor1.0.1.zip
        │   ├── rotationForest1.0.2.zip
        │   ├── stackingC1.0.1.zip
        │   └── userClassifier1.0.3.zip
        ├── results
```


