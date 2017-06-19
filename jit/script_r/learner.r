### following RWeka packages can be downloaded in site: http://sourceforge.net/projects/weka/files/weka-packages/
# RWeka::WPM("install-package", "./packages/linearForwardSelection1.0.1.zip")
# RWeka::WPM("install-package", "./packages/ridor1.0.1.zip")
# RWeka::WPM("install-package", "./packages/EMImputation1.0.1.zip")
# RWeka::WPM("install-package", "./packages/citationKNN1.0.1.zip")
# RWeka::WPM("install-package", "./packages/isotonicRegression1.0.1.zip")
# RWeka::WPM("install-package", "./packages/paceRegression1.0.1.zip")
# RWeka::WPM("install-package", "./packages/leastMedSquared1.0.1.zip")
# RWeka::WPM("install-package", "./packages/RBFNetwork1.0.8.zip")
# RWeka::WPM("install-package", "./packages/conjunctiveRule1.0.4.zip")
# RWeka::WPM("install-package", "./packages/rotationForest1.0.2.zip")

getWekaClassifier <- function(name) {
  WC <- NULL
  
  if (name %in% c("LinearRegression", "Logistic", "SMO", "IBk", "LBR", "AdaBoostM1", "Bagging", "LogitBoost", "MultiBoostAB", "Stacking", "CostSensitiveClassifier", "JRip", "M5Rules", "OneR", "PART", "J48", "LMT", "M5P", "DecisionStump", "SimpleKMeans")) {
    WC <- get(name, asNamespace("RWeka"))
  } else {
    if (name=="NaiveBayes") {
      WC <- RWeka::make_Weka_classifier("weka/classifiers/bayes/NaiveBayes")
    } else if (name %in% c("GaussianProcesses", "SimpleLogistic", "SimpleLinearRegression", "IsotonicRegression", "LeastMedSq", "MultilayerPerceptron", "PaceRegression", "PLSClassifier", "RBFNetwork", "SMOreg")) {
      WC <- RWeka::make_Weka_classifier(paste("weka/classifiers/functions/", name, sep=""))
    } else if (name %in% c("KStar", "LWL")) {
      WC <- RWeka::make_Weka_classifier(paste("weka/classifiers/lazy/", name, sep=""))
    } else if (name %in% c("Bagging", "AdaBoostM1", "RotationForest", "RandomSubSpace", "AdditiveRegression", "AttributeSelectedClassifier", "CVParameterSelection", "Grading", "GridSearch", "MultiScheme", "RegressionByDiscretization", "StackingC", "Vote")) {
      WC <- RWeka::make_Weka_classifier(paste("weka/classifiers/meta/", name, sep=""))
    } else if (name %in% c("ConjunctiveRule", "DecisionTable", "ZeroR", "Ridor")) {
      WC <- RWeka::make_Weka_classifier(paste("weka/classifiers/rules/", name, sep=""))
    } else if (name %in% c("REPTree", "UserClassifier", "RandomForest")) {
      WC <- RWeka::make_Weka_classifier(paste("weka/classifiers/trees/", name, sep=""))
    }
  }
  return(WC)
}

