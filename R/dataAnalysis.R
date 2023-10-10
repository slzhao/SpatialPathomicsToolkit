plotCorHeatmap=function(featureDataList,reductionName = paste0("ALL","_PCA"),PCs=1:3) {
  library(ComplexHeatmap)


  featureLoadings = data.frame(featureDataList[["Reductions"]][[reductionName]]$feature.loadings[,PCs])
  featureLoadingsMaxValue = max(abs(range(featureLoadings[,PCs])))
  featureLoadingsColorFun = circlize::colorRamp2(
    c(-featureLoadingsMaxValue, 0, featureLoadingsMaxValue),
    c("green", "white", "red")
  )

  corMatrix=(featureDataList[["FeatureCors"]][[reductionName]])
  featureNames=row.names(corMatrix)
  featureNamesToType=(featureDataList[["FeatureType"]][featureNames])
  featureNamesToTypeTable=table(featureNamesToType)

  col_list <- setNames(rep(list(featureLoadingsColorFun), ncol(featureLoadings)), colnames(featureLoadings))
  if (length(featureNamesToTypeTable)>1) {
    FeatureTypeColorFun=RColorBrewer::brewer.pal(length(featureNamesToTypeTable),"Set1")[1:length(featureNamesToTypeTable)]
    names(FeatureTypeColorFun)=names(featureNamesToTypeTable)
    col_list=c(col_list,FeatureType=list(FeatureTypeColorFun))

    annotationDf=data.frame(featureLoadings[featureNames,],FeatureType=featureNamesToType)
  } else {
    annotationDf=data.frame(featureLoadings)
  }
  ha = columnAnnotation(
    df = annotationDf,
    col = col_list,
    annotation_legend_param = list(PC1 = list(title = "PC")),
    show_legend = c(TRUE,c(rep(FALSE, ncol(featureLoadings)-1)))
  )

  p = Heatmap(
    corMatrix,
    top_annotation = ha,
    show_row_names = FALSE,
    show_column_names = FALSE,
    heatmap_legend_param = list(title = "Correlation")
  )
  draw(p,merge_legend = TRUE)

}


testFeatureDiff=function(featureDataList,dataName="FeatureDataNormlized",groupColumn="Sample",featureNames=NULL) {
  groupToTest=featureDataList[["SampleData"]][,groupColumn]

  if (dataName %in% names(featureDataList)) {
    dataToTest=featureDataList[[dataName]][row.names(featureDataList[["SampleData"]]),]
  } else {
    dataToTest=featureDataList[["Reductions"]][[dataName]]$cell.embeddings[row.names(featureDataList[["SampleData"]]),]
  }
  if (!is.null(featureNames)) {
    dataToTest=dataToTest[,featureNames]
  }

  resultAll=data.frame()
  #add c(pValue=pValue,groupMean=groupMean) to data frame resultAll
  for (i in 1:ncol(dataToTest)) {
    groupMean=tapply(dataToTest[,i],groupToTest,mean)
    pValue=kruskal.test(dataToTest[,i]~groupToTest)$p.value
    resultAll=rbind(resultAll,c(pValue=pValue,groupMean=groupMean))
  }


  resultAll=data.frame()
  for (i in 1:ncol(dataToTest)) {
    groupMean=tapply(dataToTest[,i],groupToTest,mean)
    pValue=kruskal.test(dataToTest[,i]~groupToTest)$p.value
    resultAll=rbind(resultAll,c(pValue=pValue,groupMean=groupMean))
  }
  colnames(resultAll)=c("pValue",names(groupMean))
  row.names(resultAll)=colnames(dataToTest)

  resultAll$pAdj=p.adjust(resultAll$pValue,method="BH")

  return(resultAll)
}

plotFeatureByGroup=function(featureDataList,featureNames,dataName="FeatureDataNormlized",groupColumn="Sample") {
  groupToTest=featureDataList[["SampleData"]][,groupColumn]

  if (dataName %in% names(featureDataList)) {
    dataToTest=featureDataList[[dataName]][row.names(featureDataList[["SampleData"]]),]
  } else {
    dataToTest=featureDataList[["Reductions"]][[dataName]]$cell.embeddings[row.names(featureDataList[["SampleData"]]),]
  }
  dataToTest=dataToTest[,featureNames]

  dataForPlot=data.frame(dataToTest,Group=groupToTest)

  #make dataForPlot wider except Group column
  dataForPlot=dataForPlot %>% pivot_longer(-Group)
  if (length(featureNames)==1) {
    dataForPlot$name=featureNames
  }
  p=ggplot(dataForPlot,aes(x=Group,y=value,fill=Group))+
    geom_boxplot()+facet_wrap(~name,scales="free_y")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+theme(legend.position = "none")

  return(p)
}



