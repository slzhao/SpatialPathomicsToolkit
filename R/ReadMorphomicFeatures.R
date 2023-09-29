ReadMorphomicFeatures=function(metaTable) {
  featureTableList=list()
  for (i in 1:nrow(metaTable)) {
    sampleName=metaTable[i,1]
    fileName=metaTable[i,2]
    featureTable=read.csv(fileName)

    #colnames bugfix
    colnames(featureTable)=gsub("_podocyte_nuclei","",colnames(featureTable))

    featureTypes=gsub("_.*","",colnames(featureTable)[-c(1:2)])

    #location and number to meta table
    featuresToMetaInd=which(featureTypes=="Number" |featureTypes=="Location")
    featuresToMetaTable=featureTable[,c(1:2,featuresToMetaInd+2)]
    featuresToMetaTable=data.frame(as.vector(metaTable[i,]),featuresToMetaTable)
    row.names(featuresToMetaTable)=paste0(featuresToMetaTable$File,
                                          "_",
                                          featuresToMetaTable$ImageNumber,
                                          "_",
                                          featuresToMetaTable$ObjectNumber)

    featureTableData=featureTable[,-c(1:2,featuresToMetaInd+2)]
    colnames(featureTableData)=make.names(gsub("^[A-Za-z]+_","",colnames(featureTableData)))
    row.names(featureTableData)=row.names(featuresToMetaTable)

    featureTypes=featureTypes[-c(featuresToMetaInd)]
    names(featureTypes)=colnames(featureTableData)

    featureTableList[["FeatureData"]]=rbind(featureTableList[["FeatureData"]],featureTableData)
    featureTableList[["SampleData"]]=rbind(featureTableList[["SampleData"]],featuresToMetaTable)
    featureTableList[["FeatureType"]]=featureTypes
  }
  return(featureTableList)
}

#Normalization
MorphomicFeaturesNormlization=function(featureDataAll,nUnique=10,CellType=NULL,Sample=NULL,featureType=NULL) {
  #library(tidyverse)
  #source("D:/source/r_cqs/myPkg/R/dataTransfomation.R")

  featureDataSubset=DataSelection(featureDataAll,CellType=CellType,Sample=Sample,featureType=featureType)


  dataTable=featureDataSubset[["FeatureData"]]
  #SampleDataTable=featureDataSubset[["SampleData"]]
  #featureTypeTable=featureDataSubset[["FeatureType"]]

  dataTableFiltered=dataTable

  print("Remove variables less than 10 unique values")
  temp=apply(dataTable,2,function(x) length(unique(x)))
  featuresToRemoveInd=which(temp<=nUnique)
  if (length(featuresToRemoveInd)>0) {
    featureToRemoveNames=colnames(dataTable)[featuresToRemoveInd]
    dataTableFiltered=dataTable[,-featuresToRemoveInd]

    featureDataSubset[["FeatureData"]]=featureDataSubset[["FeatureData"]][,-featuresToRemoveInd]
    featureDataSubset[["FeatureType"]]=featureDataSubset[["FeatureType"]][-featuresToRemoveInd]

    print(paste0(length(featuresToRemoveInd)," features removed"))
    print("Number of features after selection:")
    print(table(featureDataSubset[["FeatureType"]]))
  } else {
    dataTableFiltered=dataTable
  }


  print("Remove variables with NA")
  temp=apply(dataTableFiltered,2,function(x) length(which(is.na(x))))
  featuresToRemoveInd=which(temp>0)
  if (length(featuresToRemoveInd)>0) {
    featureToRemoveNames=colnames(dataTableFiltered)[featuresToRemoveInd]
    dataTableFiltered=dataTableFiltered[,-featuresToRemoveInd]

    featureDataSubset[["FeatureData"]]=featureDataSubset[["FeatureData"]][,-featuresToRemoveInd]
    featureDataSubset[["FeatureType"]]=featureDataSubset[["FeatureType"]][-featuresToRemoveInd]

    print(paste0(length(featuresToRemoveInd)," features removed"))
    print(paste(featureToRemoveNames,collapse=";"))
    print("Number of features after selection:")
    print(table(featureDataSubset[["FeatureType"]]))
  }

  print("Remove variables less than 2% subjects with unique values")
  temp=apply(dataTableFiltered,2,function(x) quantile(x,c(0.01,0.99)))
  featuresToRemoveInd=which(temp[1,]==temp[2,])
  if (length(featuresToRemoveInd)>0) {
    featureToRemoveNames=colnames(dataTableFiltered)[featuresToRemoveInd]
    dataTableFiltered=dataTableFiltered[,-featuresToRemoveInd]

    featureDataSubset[["FeatureData"]]=featureDataSubset[["FeatureData"]][,-featuresToRemoveInd]
    featureDataSubset[["FeatureType"]]=featureDataSubset[["FeatureType"]][-featuresToRemoveInd]

    print(paste0(length(featuresToRemoveInd)," features removed"))
    print(paste(featureToRemoveNames,collapse=";"))
    print("Number of features after selection:")
    print(table(featureDataSubset[["FeatureType"]]))
  }

  rawDataValuesNormlizedResult=dataTransfomation(data.frame(dataTableFiltered))
  rawDataValuesNormlized=rawDataValuesNormlizedResult[[1]]
  rawDataValuesNormlizedSummary=rawDataValuesNormlizedResult[[2]]

  # if (!is.null(featureTypeTable)) {
  #   rawDataValuesNormlizedSummary=rawDataValuesNormlizedSummary %>% left_join(rownames_to_column(data.frame(featureTypeTable)),by=c("feature"="rowname"))
  # }

  print("Normalization done")
  # print(paste0("Number of variables: ",ncol(rawDataValuesNormlized)))
  # print(paste0("Number of cell types included: ",(table(rawDataValuesNormlized$CellType))))
  # print(paste0("Number of samples included: ",(table(rawDataValuesNormlized$Sample))))

  return(c(featureDataSubset,FeatureDataNormlized=list(rawDataValuesNormlized),FeatureDataNormlizedSummary=list(rawDataValuesNormlizedSummary)))
}




