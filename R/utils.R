DataSelection=function(featureDataAll,CellType=NULL,Sample=NULL,featureType=NULL) {
  featureDataSubset=featureDataAll
  SampleDataAllTable=featureDataAll[["SampleData"]]

  if (!is.null(CellType) & !is.null(Sample)) {
    selectedSamples=which((SampleDataAllTable$Sample %in% Sample) & (SampleDataAllTable$CellType %in% CellType))
    print(paste0(paste(Sample,collapse=";")," were selected as Sample"))
    print(paste0(paste(CellType,collapse=";")," were selected as CellType"))
  } else if (!is.null(CellType)) {
    selectedSamples=which(SampleDataAllTable$CellType %in% CellType)
    print(paste0(paste(CellType,collapse=";")," were selected as CellType"))
  } else if (!is.null(Sample)) {
    selectedSamples=which((SampleDataAllTable$Sample %in% Sample))
    print(paste0(paste(Sample,collapse=";")," were selected as Sample"))
  } else {
    selectedSamples=NULL
  }
  if (!is.null(selectedSamples)) {
    featureDataSubset[["FeatureData"]]=featureDataSubset[["FeatureData"]][selectedSamples,]
    featureDataSubset[["SampleData"]]=featureDataSubset[["SampleData"]][selectedSamples,]

    if ("FeatureDataNormlized" %in% names(featureDataSubset)) {
      featureDataSubset[["FeatureDataNormlized"]]=featureDataSubset[["FeatureDataNormlized"]][selectedSamples,]
    }

    print("Number of Sample/CellType after selection:")
    print(table(featureDataSubset[["SampleData"]]$Sample,featureDataSubset[["SampleData"]]$CellType))

  }

  if (!is.null(featureType)) {
    selectedFeatures=which(featureDataSubset[["FeatureType"]] %in% featureType)
    featureDataSubset[["FeatureData"]]=featureDataSubset[["FeatureData"]][,selectedFeatures]
    featureDataSubset[["FeatureType"]]=featureDataSubset[["FeatureType"]][selectedFeatures]

    if ("FeatureDataNormlized" %in% names(featureDataSubset)) {
      featureDataSubset[["FeatureDataNormlized"]]=featureDataSubset[["FeatureDataNormlized"]][,selectedFeatures]
    }
    if ("FeatureDataNormlizedSummary" %in% names(featureDataSubset)) {
      featureDataSubset[["FeatureDataNormlizedSummary"]]=featureDataSubset[["FeatureDataNormlizedSummary"]][selectedFeatures,]
    }

    print(paste0(paste(featureType,collapse=";")," were selected as featureType"))
    print("Number of features after selection:")
    print(table(featureDataSubset[["FeatureType"]]))
  }

  return(featureDataSubset)
}
