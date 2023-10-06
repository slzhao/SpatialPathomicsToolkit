featureDataPca=function(featureDataList,dataName="FeatureDataNormlized") {

  featureData=featureDataList[[dataName]]
  #remove features with 0 SD
  featureData=featureData[,apply(featureData,2,sd)!=0]

  featureDataScaled=scale(featureData)

  pca.results <- prcomp(x = featureDataScaled,scale=FALSE,center = FALSE )
  feature.loadings <- pca.results$rotation
  cell.embeddings <- pca.results$x
  importance <- summary(pca.results)$importance

  row.names(cell.embeddings)=row.names(featureData)

  return(list(
    feature.loadings=feature.loadings,
    cell.embeddings=cell.embeddings,
    importance=importance
  ))
}

featureDataSparsePca=function(featureDataList,dataName="FeatureDataNormlized",
                              k=NULL, max_iter=1000,alpha = 1e-04) {
  library(sparsepca)

  featureData=featureDataList[[dataName]]
  #remove features with 0 SD
  featureData=featureData[,apply(featureData,2,sd)!=0]

  featureDataScaled=scale(featureData)

  rspca.results <- rspca(featureDataScaled, k=k,verbose = FALSE,
                         max_iter=max_iter,center=FALSE, scale=FALSE,alpha = alpha)
  row.names(rspca.results$loadings)=colnames(featureDataScaled)
  colnames(rspca.results$loadings)=paste0("PC",1:ncol(rspca.results$loadings))
  colnames(rspca.results$scores)= colnames(rspca.results$loadings)

  feature.loadings <- rspca.results$loadings
  cell.embeddings <- rspca.results$scores
  importance <- data.frame(summary(rspca.results))

  return(list(
    feature.loadings=feature.loadings,
    cell.embeddings=cell.embeddings,
    importance=importance
  ))
}



featureDataPcaAll=function(featureDataList,dataName="FeatureDataNormlized") {
  reductionName=paste0("ALL","_PCA")
  featureDataSubset=featureDataList
  pcaResult=featureDataPca(featureDataSubset,dataName=dataName)
  featureDataList[["Reductions"]][[reductionName]]=pcaResult
  corResult=cor((featureDataSubset[[dataName]]), method = "sp")
  featureDataList[["FeatureCors"]][[reductionName]]=corResult

  FeatureTypeAll=unique(featureDataList[["FeatureType"]])
  CellTypeAll=unique(featureDataList[["SampleData"]]$CellType)

  for (FeatureTypeOne in FeatureTypeAll) {
    reductionName=paste0(FeatureTypeOne,"_PCA")

    featureDataSubset=DataSelection(featureDataList,featureType=FeatureTypeOne)
    pcaResult=featureDataPca(featureDataSubset,dataName=dataName)
    featureDataList[["Reductions"]][[reductionName]]=pcaResult

    corResult=cor((featureDataSubset[[dataName]]), method = "sp")
    featureDataList[["FeatureCors"]][[reductionName]]=corResult

  }

  for (CellTypeOne in CellTypeAll) {
    reductionName=paste0(CellTypeOne,"_PCA")

    featureDataSubset=DataSelection(featureDataList,CellType =CellTypeOne)
    pcaResult=featureDataPca(featureDataSubset,dataName=dataName)
    featureDataList[["Reductions"]][[reductionName]]=pcaResult

    corResult=cor((featureDataSubset[[dataName]]), method = "sp")
    featureDataList[["FeatureCors"]][[reductionName]]=corResult
  }

  for (CellTypeOne in CellTypeAll) {
    for (FeatureTypeOne in FeatureTypeAll) {

      reductionName=paste0(CellTypeOne,"_",FeatureTypeOne,"_PCA")

      featureDataSubset=DataSelection(featureDataList,CellType =CellTypeOne,featureType=FeatureTypeOne)
      pcaResult=featureDataPca(featureDataSubset,dataName=dataName)
      featureDataList[["Reductions"]][[reductionName]]=pcaResult

      corResult=cor((featureDataSubset[[dataName]]), method = "sp")
      featureDataList[["FeatureCors"]][[reductionName]]=corResult
    }
  }

  return(featureDataList)
}


plotPcaBarplot=function(featureDataList,reductionName = paste0("ALL","_PCA"),PCs=1:10, loadingCut = NULL,textSize=10) {
  library("viridis")

  featureLoadings = featureDataList[["Reductions"]][[reductionName]]$feature.loadings
  if (ncol(featureLoadings)<=length(PCs)) {
    PCs=1:ncol(featureLoadings)
  }

  featureLoadings = featureLoadings[,PCs]

  if (!is.null(loadingCut)) { #loadingCut defined, color by number of genes >=loadingCut in this PC
    featureLoadingsSignificant = apply(featureLoadings, 2, function(x)
      length(which(abs(x) >= loadingCut)))
    colorLegend=paste0("# Of Features")
  } else { #loadingCut NOT defined, color by number of genes Max in this PC
    temp=apply(featureLoadings,1,function(x) which.max(abs(x)))
    featureLoadingsSignificant=structure(rep(0,length(PCs)),names=as.character(PCs))
    featureLoadingsSignificant[names(table(temp))]=as.vector(table(temp))
    colorLegend="# Of Features"
  }

  dataForPlot = data.frame(
    Variance = featureDataList[["Reductions"]][[reductionName]]$importance["Proportion of Variance",PCs],
    CumulativeVariance = featureDataList[["Reductions"]][[reductionName]]$importance["Cumulative Proportion",PCs],
    Name = factor(paste0("PC", PCs), levels = (paste0("PC", PCs))),
    NumberOfFeatures = featureLoadingsSignificant[PCs],
    PercentOfFeatures = featureLoadingsSignificant[PCs]/nrow(featureLoadings)
  )
  # scaling_factor <-
  #   max(dataForPlot$Variance) / max(dataForPlot$NumberOfFeatures)
  # scaling_factor <-
  #   max(dataForPlot$Variance) / max(dataForPlot$PercentOfFeatures)
  scaling_factor <-
    max(dataForPlot$Variance) / max(dataForPlot$CumulativeVariance)

  p = ggplot(dataForPlot) +
    geom_col(aes(x = Name, y = Variance,fill=NumberOfFeatures)) +
    geom_line(
      aes(
        x = Name,
        #        y = NumberOfFeatures * scaling_factor,
        y = CumulativeVariance * scaling_factor,
        group = 1
      ),
      color = 'red',
      linewidth = 1
    ) +
    scale_y_continuous(sec.axis = sec_axis( ~ . / scaling_factor, name = "Cumulative Proportion of Variance")) +
    labs(y = "Proportion of Variance", x = "") +
    theme_bw()
  p=    p + theme(
    axis.title = element_text(size = textSize),
    axis.text = element_text(size = textSize),
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    legend.text = element_text(size = textSize),
    legend.title = element_text(size = textSize)
  )

  #p=p+ guides(fill=guide_legend(title=colorLegend))+scale_fill_gradient(low="blue", high="yellow")
  #color_palette <- rev(RColorBrewer::brewer.pal(7, 'RdBu'))
  #p=p+ guides(fill=guide_legend(title=colorLegend))+scale_fill_gradientn(colors = color_palette)

  p=p+ guides(fill=guide_legend(title=colorLegend))+scale_fill_viridis()

  return(
    p
  )
}



plotPCsByCellType=function(featureDataAllNormlized,reductionName,reductionNameToComp=NULL,PCs=1:4) {
  if (is.null(reductionName)) {
    reductionName="ALL"
  }
  if (is.null(reductionNameToComp)) {
    reductionNameToComp=grep(paste0("_",reductionName),names(featureDataAllNormlized[["Reductions"]]),value=TRUE)
  }

  if (length(PCs)>ncol(featureDataAllNormlized[["Reductions"]][[reductionName]]$cell.embeddings)) {
    PCs=PCs[1:ncol(featureDataAllNormlized[["Reductions"]][[reductionName]]$cell.embeddings)]
  }
  pList=list()
  for (reductionNameToCompOne in reductionNameToComp) {
    for (i in seq_along(PCs)) {
      temp1=featureDataAllNormlized[["Reductions"]][[reductionName]]$cell.embeddings[,PCs[i]]
      temp2=featureDataAllNormlized[["Reductions"]][[reductionNameToCompOne]]$cell.embeddings[,PCs[i]]
      dataForPlot=data.frame(temp1[names(temp2)],temp2)
      names(dataForPlot)=c(reductionName,reductionNameToCompOne)
      p=ggplot(dataForPlot,aes(x=!!sym(reductionName),y=!!sym(reductionNameToCompOne)))+geom_point()+ggtitle(paste0(reductionNameToCompOne," PC",PCs[i]))
      p=p+ggtitle(paste0(reductionName," vs ",reductionNameToCompOne," PC",PCs[i]))+theme_bw()
      pList=c(pList,list(p))
    }
  }
  patchwork::wrap_plots(pList,ncol=length(PCs))
}

plotPCsLoading=function(featureLoadings,topN=10,ZeroLoadingToNA=TRUE) {
  if (ZeroLoadingToNA) {
    featureLoadings[featureLoadings==0]=NA
  }

  #extract topN features for each reduction
  if (nrow(featureLoadings)<topN) {
    selectedFeatures=row.names(featureLoadings)
  } else {
    temp = apply(featureLoadings, 2, function(x)
      names(sort(abs(x), decreasing = TRUE)[1:topN]))
    selectedFeatures=na.omit(unique(as.vector(temp)))
  }

  dataForPlot=featureLoadings[selectedFeatures,]
  #make dataForPlot into longer form by pivot_longer
  dataForPlot=as.data.frame(dataForPlot) %>% rownames_to_column(var = "feature") %>%
    pivot_longer(-feature,names_to = "Reduction",values_to = "Loading") %>% mutate(
      AbsLoading=abs(Loading),
      SignLoading=sign(Loading)
    )
  #make barplot in dataForPlot, facet by reduction column
  p=ggplot(dataForPlot,aes(x=AbsLoading,y=feature,fill=factor(SignLoading)))+geom_col()+
    facet_grid(.~Reduction)+theme_bw()
  #change the fill color in p to pretty colors
  p=p+scale_fill_manual(values=c("green","red"))

  return(p)
}


plotPCsLoadingByCellType=function(featureDataAllNormlized,reductionName,reductionNameToComp=NULL,topN=10,PCs=1:4) {
  if (is.null(reductionName)) {
    reductionName="ALL"
  }
  if (is.null(reductionNameToComp)) {
    reductionNameToComp=grep(paste0("_",reductionName),names(featureDataAllNormlized[["Reductions"]]),value=TRUE)
  }
  if (length(PCs)>ncol(featureDataAllNormlized[["Reductions"]][[reductionName]]$feature.loadings)) {
    PCs=PCs[1:ncol(featureDataAllNormlized[["Reductions"]][[reductionName]]$feature.loadings)]
  }

  pList=list()
  for (i in seq_along(PCs)) {
    featureLoadings=NULL
    featureLoadingsOne=featureDataAllNormlized[["Reductions"]][[reductionName]]$feature.loadings[,PCs[i]]
    featureLoadings=cbind(featureLoadings,featureLoadingsOne)
    for (reductionNameToCompOne in reductionNameToComp) {
      featureLoadingsOne=featureDataAllNormlized[["Reductions"]][[reductionNameToCompOne]]$feature.loadings[,PCs[i]]
      featureLoadings=cbind(featureLoadings,featureLoadingsOne[row.names(featureLoadings)])
    }
    colnames(featureLoadings)=c(reductionName,reductionNameToComp)

    p=plotPCsLoading(featureLoadings,topN=topN)
    # #extract topN features for each reduction
    # if (nrow(featureLoadings)<topN) {
    #   selectedFeatures=row.names(featureLoadings)
    # } else {
    #   temp = apply(featureLoadings, 2, function(x)
    #     names(sort(abs(x), decreasing = TRUE)[1:topN]))
    #   selectedFeatures=unique(as.vector(temp))
    # }
    #
    # dataForPlot=featureLoadings[selectedFeatures,]
    # #make dataForPlot into longer form by pivot_longer
    # dataForPlot=as.data.frame(dataForPlot) %>% rownames_to_column(var = "feature") %>%
    #   pivot_longer(-feature,names_to = "Reduction",values_to = "Loading") %>% mutate(
    #     AbsLoading=abs(Loading),
    #     SignLoading=sign(Loading)
    #   )
    # #make barplot in dataForPlot, facet by reduction column
    # p=ggplot(dataForPlot,aes(x=AbsLoading,y=feature,fill=factor(SignLoading)))+geom_col()+
    #   facet_grid(.~Reduction)+theme_bw()
    # #change the fill color in p to pretty colors
    # p=p+scale_fill_manual(values=c("green","red"))

    p=p+ggtitle(paste0("PC",PCs[i]))

    pList=c(pList,list(p))
  }
  patchwork::wrap_plots(pList,ncol=1)
}
