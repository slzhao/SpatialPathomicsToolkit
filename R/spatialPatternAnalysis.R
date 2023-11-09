require(tidyverse)
require(Seurat)
require(png)
PrepareFeaturesSeuratObj <- function(featureData,
                                     group = NULL,
                                     subset = NULL,
                                     assay = NULL,
                                     image = NULL,
                                     image.id = NULL,
                                     positionTable,
                                     project = NULL,
                                     scalefactors = NULL,
                                     spot.radius = NULL,
                                     Class = ifelse(is.null(image), "SlideSeq", "VisiumV1"), ...) {
  #overall data
  if(is.null(group)) {
    pca_id <- "ALL_PCA"
  } else {
    pca_id <- paste0(group, "_PCA")
  }

  meta <- featureData$SampleData
  countTable <- featureData$FeatureData
  countTable.norm <- featureData$FeatureDataNormlized
  cell_emb <- featureData$Reductions[[pca_id]]$cell.embeddings

  if(!is.null(subset)){
    meta <- meta[which(rownames(meta) %in% subset),]
    countTable <- countTable[which(rownames(countTable) %in% subset),]
    countTable.norm <- countTable.norm[which(rownames(countTable.norm) %in% subset),]
    positionTable <- positionTable[which(rownames(positionTable) %in% subset),]
    cell_emb <- cell_emb[which(rownames(cell_emb) %in% subset),]
  }
  object <- CreateSeuratObject(t(countTable),
                               project = project,
                               assay = assay)
  DefaultAssay(object) <- assay
  assay_key <- tolower(assay)
  assay_key <- str_remove_all(assay_key, "_")
  Key(object@assays[[assay]]) <- paste0(assay_key, "_")
  #load image
  if(is.null(image.id)){
    image.id <- "slice1"
  }
  if(is.null(image)) {
    object[['slice1']] =  new(
      Class = Class,
      assay = assay,
      coordinates = positionTable
    )
  } else {
    if(class(image) == "array" & dim(image)[3] == 3) {
      img <- image
      object[['slice1']] <- new(Class = Class,
                                assay = assay,
                                image = img,
                                scale.factors = scalefactors,
                                spot.radius = spot.radius,
                                coordinates = positionTable)
    } else {
      print("Error: Please provide existing H&E image file or
            the image array!\n")
    }
  }
  object@meta.data <- cbind(object@meta.data, meta)
  object@assays[[assay]]@data <- t(countTable.norm)
  #load reduction
  object@reductions[[pca_id]] <- CreateDimReducObject(
    embeddings = cell_emb,
    loadings = featureData$Reductions[[pca_id]]$feature.loadings,
    stdev = featureData$Reductions[[pca_id]]$importance["Standard deviation",],
    key = "PC_",
    assay = assay
  )
  #load each feature type
  FeatureTypes <- featureData$FeatureType
  for(i in 1:length(table(FeatureTypes))){
    type <- names(table(FeatureTypes)[i])
    features <- names(FeatureTypes[which(FeatureTypes == type)])
    sub_count <- countTable[,features]
    sub_count.norm <- countTable.norm[,features]
    if(is.null(group)) {
      assay_id <- paste0(assay, "_", type)
      pca_id <- paste0(type, "_PCA")
    } else {
      assay_id <- paste0(assay, "_", group, "_", type)
      pca_id <- paste0(group, "_", type, "_PCA")
    }
    object@assays[[assay_id]] <- CreateAssayObject(counts = t(sub_count))
    object@assays[[assay_id]]@data <- t(sub_count.norm)
    DefaultAssay(object) <- assay_id
    assay_key <- tolower(assay_id)
    assay_key <- str_remove_all(assay_key, "_")
    Key(object@assays[[assay_id]]) <- paste0(assay_key, "_")

    cell_emb <- featureData$Reductions[[pca_id]]$cell.embeddings
    if(!is.null(subset)){
      cell_emb <- cell_emb[which(rownames(cell_emb) %in% subset),]
    }

    object@reductions[[pca_id]] <- CreateDimReducObject(
      embeddings = cell_emb,
      loadings = featureData$Reductions[[pca_id]]$feature.loadings,
      stdev = featureData$Reductions[[pca_id]]$importance["Standard deviation",],
      key = "PC_",
      assay = assay_id
    )
  }
  return(object)
}

PrepareFeaturesSeuratObjVisium <- function(featureData,
                                           st.object,
                                           project,
                                           patch.size = 512,
                                           image.id = NULL,
                                           group = NULL,
                                           subset = NULL,
                                           assay = NULL,
                                           matchTable = NULL,
                                           ...) {
  #overall data
  if(is.null(image.id)){
    image.id <- "slice1"
  }
  img <- st.object@images[[image.id]]@image
  scalefactors <- st.object@images[[image.id]]@scale.factors
  spot.radius <- st.object@images[[image.id]]@scale.factors$spot
  patch.scale.factor <- patch.size / (2 * spot.radius)

  colnames(matchTable) <- c("V1", "V2")
  #prepare coordinates
  meta <- featureData$SampleData
  coord_st <- st.object@images[[image.id]]@coordinates
  coord_st$barcode <- row.names(coord_st)
  coord_tmp <- meta %>%
    rownames_to_column("cellid") %>%
    mutate(cellid = str_replace_all(cellid, "_", "-")) %>%
    dplyr::filter(cellid %in% subset) %>%
    left_join(., matchTable, by = c("cellid" = "V2")) %>%
    left_join(., coord_st, by = c("V1" = "barcode")) %>%
    mutate_at(c("ObjectNumber", "Location_Center_X", "Location_Center_Y",
                "tissue", "row", "col", "imagerow", "imagecol"), as.numeric) %>%
    mutate(new_row = row + ObjectNumber/10,
           new_col = col + ObjectNumber/10,
           new_pxl_row = imagerow + (Location_Center_X - patch.size / 2)/patch.scale.factor,
           new_pxl_col = imagecol + (Location_Center_Y - patch.size / 2)/patch.scale.factor) %>%
    dplyr::select(cellid, tissue, new_row, new_col, new_pxl_row, new_pxl_col)
  colnames(coord_tmp) <- c("barcode", "in_tissue", "array_row", "array_col",
                           "pxl_row_in_fullres", "pxl_col_in_fullres")
  coord_new <- coord_tmp %>%
    column_to_rownames("barcode")
  colnames(coord_new) <- c("tissue", "row", "col", "imagerow", "imagecol")
  object <- PrepareFeaturesSeuratObj(featureData = featureData,
                                     group = group,
                                     subset = subset,
                                     assay = assay,
                                     image = img,
                                     positionTable = coord_new,
                                     scalefactors = scalefactors,
                                     spot.radius = spot.radius,
                                     project = project)

  return(object)
}

MergeMatrices <- function(matrix1,
                          matrix2,
                          mergeTable = NULL) {
  colnames(mergeTable) <- c("V1", "V2")
  mtx1.tmp <- matrix1[,colnames(matrix1) %in% mergeTable$V1]
  mtx2.tmp <- as.data.frame(t(matrix2)) %>%
    rownames_to_column("ID") %>%
    left_join(., mergeTable, by = c("ID" = "V2")) %>%
    gather("feature", "value", rownames(matrix2)) %>%
    group_by(feature, V1) %>%
    dplyr::summarize(mean=mean(value), n = n(), median = median(value), .groups = 'drop') %>%
    dplyr::select(feature, V1, mean) %>%
    spread(V1, mean) %>%
    column_to_rownames("feature")
  mtx2.tmp <- as.matrix(mtx2.tmp)
  return(list(mtx1.tmp, mtx2.tmp))
}

MergeObjects <- function(object1,
                         object2,
                         mergeTable = NULL) {
  object <- object1
  slot1 <- names(object1@assays)[1]
  mtx1 <- as.matrix(object1@assays[[slot1]]@counts)
  colnames(mergeTable) <- c("V1", "V2")
  cells_to_keep <- unique(mergeTable[,1])
  object <- object[, Cells(object) %in% cells_to_keep]
  assay_slots <- names(object2@assays)
  pca_slots <- names(object2@reductions)
  for(i in 1:length(assay_slots)){
    slot2 <- assay_slots[i]
    mtx2 <- as.matrix(object2@assays[[slot2]]@counts)
    mtx2.norm <- as.matrix(object2@assays[[slot2]]@data)

    new.mtx2 <- MergeMatrices(mtx1, mtx2, mergeTable)
    new.mtx2.norm <- MergeMatrices(mtx1, mtx2.norm, mergeTable)
    object@assays[[slot2]] <- CreateAssayObject(counts = new.mtx2[[2]])
    object@assays[[slot2]]@data <- new.mtx2.norm[[2]]
    Key(object@assays[[slot2]]) <- Key(object2@assays[[slot2]])
  }
  for(i in 1:length(pca_slots)){
    pca.slot2 <- pca_slots[i]
    mtx2 <- as.matrix(t(object2@reductions[[pca.slot2]]@cell.embeddings))

    new.mtx2 <- MergeMatrices(mtx1, mtx2, mergeTable)
    new_cell_emb <- new.mtx2[[2]]
    #object@reductions[[pca.slot2]]$cell.embeddings <- new_cell_emb
    object@reductions[[pca.slot2]] <- CreateDimReducObject(
      embeddings = t(new_cell_emb),
      loadings = object2@reductions[[pca.slot2]]@feature.loadings,
      stdev = object2@reductions[[pca.slot2]]@stdev,
      key = "PC_",
      assay = slot2
    )
  }
  return(object)
}















