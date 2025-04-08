read_matrice_matrix <- function(h5data) {
  if ("data" %in% names(h5data) && "indices" %in% names(h5data) && "indptr" %in% names(h5data)) {
    i <- h5data[["indices"]][]
    p <- h5data[["indptr"]][]
    x <- h5data[["data"]][]

    X <- if ("csr_matrix" %in% getEncodingType(h5data)) {
      Matrix::sparseMatrix(i = i, p = p, x = x, dims = rev(getShape(h5data)), index1 = FALSE)
    } else {
      Matrix::t(Matrix::sparseMatrix(j = i, p = p, x = x, dims = getShape(h5data), index1 = FALSE))
    }

    return(X)
  } else {
    return(h5data[])
  }
}

# name <- names(h5layer)[[1]]
# h5data <- group_h5layer
read_matrix <- function(h5layer, useBPcells = FALSE, cellNames = NULL, geneNames = NULL) {
  data_list <- list()
  for (name in names(h5layer)) {
    group_h5layer <- h5layer[[name]]
    if ("array" %in% getEncodingType(h5layer)) {
      data_list[[name]] <- group_h5layer[, ]
    } else {
      if (useBPcells) {
        data_list[[name]] <- as(read_matrice_matrix(group_h5layer), "IterableMatrix")
      } else {
        data_list[[name]] <- read_matrice_matrix(h5data = group_h5layer)
      }
    }
  }

  # dim(data_list[[1]])
  # dim(data_list[[2]])
  # rbind(data_list[[1]], data_list[[2]])
  all_data <- do.call(cbind, data_list)
  if (!is.null(cellNames) & ncol(all_data) == length(cellNames)) {
    colnames(all_data) <- cellNames
  }
  if (!is.null(geneNames) & nrow(all_data) == length(geneNames)) {
    rownames(all_data) <- geneNames
  }
  return(all_data)
}


h5_to_X <- function(h5, assay = "RNA", layer = "rawdata", useBPcells = FALSE, useDisk = TRUE, cellNames = NULL, geneNames = NULL) {
  if (is.null(cellNames)) {
    cellNames <- h5[["names_obs"]][]
  }

  if (is.null(geneNames)) {
    if (layer == "rawdata") {
      geneNames <- h5[["var/rawvar/_index"]][]
    } else {
      geneNames <- h5[["names_var"]][]
    }
  }

  h5layer <- h5[["assay"]][[assay]][["layers"]][[layer]]
  if (length(names(h5layer)) > 200) {
    print("The number of layers is too large. using BPCells to load the data.")
    useBPcells <- TRUE
  }
  if (useBPcells) {
    if (!require(BPCells)) {
      stop("The number of layers is too large. But the BPCells package not installed.")
    }
  }

  varType <- if (layer == "rawdata") "rawvar" else "var"
  X <- read_matrix(h5layer, useBPcells = useBPcells, cellNames = cellNames, geneNames = geneNames)
  if (useDisk & useBPcells) {
    diskFile <- sprintf("./.BPCells_dont_delete/%s/%s/%s", sub(".h5", "", basename(h5$filename)), assay, layer)
    if (dir.exists(diskFile)) {
      print("The data has been cached in disk. Delete the cache file to resave the data.")
      unlink(diskFile, recursive = TRUE)
    }
    write_matrix_dir(X, diskFile)
    X <- open_matrix_dir(diskFile)
  }
  return(X)
}

h5_to_DF <- function(h5data) {
  print("Starting h5_to_DF function...")

  # Check if the _index exists and has dimensions
  if (!"_index" %in% names(h5data)) {
    stop("The h5data object does not contain an '_index' dataset.")
  }
  if (h5data[["_index"]]$dims == 0) {
    print("h5data has no dimensions (dims == 0). Returning NULL.")
    return(NULL)
  }

  # Retrieve row names from _index
  rownamesStr <- h5data[["_index"]][]
  print(sprintf("Retrieved %d row names from '_index'.", length(rownamesStr)))
  
  # Initialize empty data frame for the result
  df <- data.frame()
  
  # Loop through each dataset in the group
  for (name in names(h5data)) {
    print(paste("Processing dataset:", name))
    if (name == "_index") {
      next  # row names already processed
    }
    
    # Determine encoding type with extra error handling
    encodingType <- tryCatch({
      getEncodingType(h5data[[name]])
    }, error = function(e) {
      print(paste("Error determining encoding type for", name, ":", e$message))
      return(NULL)
    })
    
    if (is.null(encodingType)) {
      print(paste("Skipping", name, "due to missing encoding type."))
      next
    }
    
    # Process the dataset based on its encoding type
    if ("categorical" %in% encodingType) {
      print(paste("Dataset", name, "is categorical."))
      values_attr <- h5data[[name]]
      labelName <- tryCatch({
        values_attr[["categories"]][]
      }, error = function(e) {
        print(paste("Error retrieving categories for", name, ":", e$message))
        return(NULL)
      })
      if (is.null(labelName)) next
      
      values <- tryCatch({
        values_attr[["codes"]][]
      }, error = function(e) {
        print(paste("Error retrieving codes for", name, ":", e$message))
        return(NULL)
      })
      if (is.null(values)) next
      values <- factor(as.integer(values), labels = labelName)
    } else if (encodingType %in% c("array", "string-array")) {
      print(paste("Dataset", name, "is an array."))
      values <- tryCatch({
        h5data[[name]][]
      }, error = function(e) {
        print(paste("Error retrieving array data for", name, ":", e$message))
        return(NULL)
      })
      if (is.null(values)) next
    } else {
      print(paste("Skipping", name, "due to unknown encoding type:", encodingType))
      next
    }
    
    # Check if the length of the extracted values matches the row names length
    if (length(values) != length(rownamesStr)) {
      warning(sprintf("Length mismatch in dataset '%s': extracted %d values vs row names length %d.", 
                      name, length(values), length(rownamesStr)))
    }
    
    # Try to add the data as a column to the data frame
    tryCatch({
      df <- addDF(df, setNames(data.frame(values), name), "col")
      print(paste("Added dataset", name, "to the data frame."))
    }, error = function(e) {
      print(paste("Error adding column", name, "to data frame:", e$message))
      # Continue with the next dataset
    })
  }
  
  print(sprintf("Data frame assembled with %d rows and %d columns.", nrow(df), ncol(df)))
  print(sprintf("Row names length is %d.", length(rownamesStr)))
  
  # Before assigning row names, double-check that the lengths match
  if (nrow(df) != length(rownamesStr)) {
    stop(sprintf("Mismatch: number of rows in df (%d) does not match row names length (%d)! Check source datasets for inconsistencies.", 
                 nrow(df), length(rownamesStr)))
  }
  rownames(df) <- rownamesStr
  
  # Optionally reorder columns if the 'column-order' attribute exists
  if ("column-order" %in% hdf5r::h5attr_names(h5data)) {
    print("Reordering columns based on 'column-order' attribute...")
    colnamesOrder <- tryCatch({
      hdf5r::h5attr(h5data, "column-order")
    }, error = function(e) {
      print("Error retrieving 'column-order' attribute.")
      return(NULL)
    })
    if (!is.null(colnamesOrder)) {
      print(paste("New column order:", paste(colnamesOrder, collapse = ", ")))
      df <- df[, colnamesOrder, drop = FALSE]
    }
  }
  
  print("Returning final data frame from h5_to_DF.")
  return(df)
}





openH5 <- function(FileName) {
  if (!hdf5r::is_hdf5(FileName)) {
    stop("File is not a hdf5 file.")
  }
  h5 <- hdf5r::H5File$new(FileName, "r")
  if (!("var" %in% names(h5))) {
    stop("var not found.")
  }
  if (!("obs" %in% names(h5))) {
    stop("obs not found.")
  }
  if (!("assay" %in% names(h5))) {
    stop("assay not found.")
  }
  return(h5)
}

# 读取 uns 数据
h5_to_uns <- function(h5obj) {
  data <- list()

  # 检查 h5obj 是否为组
  if (is(h5obj, "H5Group")) {
    keys <- names(h5obj)
    for (key in keys) {
      # 使用 tryCatch 跳过报错项
      tryCatch(
        {
          data[[key]] <- h5_to_uns(h5obj[[key]])
        },
        error = function(e) {
          message(paste("message: Error processing key:", key, "-", e$message))
        }
      )
    }
  } else { # h5obj 是数据集
    item <- tryCatch(
      {
        h5obj[]
      },
      error = function(e) {
        message(paste("message: Error reading dataset:", e$message))
        NULL
      }
    )

    if (!is.null(item)) {
      if (is.raw(item[1])) {
        item_text <- rawToChar(item)
        # 尝试将内容作为 JSON 解析，否则保持原样
        data <- tryCatch(
          {
            jsonlite::fromJSON(item_text)
          },
          error = function(e) {
            message("message: Error parsing JSON - returning raw text")
            item_text
          }
        )
      } else {
        data <- item
      }
    }
  }
  return(data)
}


h5_to_obs <- function(h5, obsName = "obs") {
  return(h5_to_DF(h5[["obs"]]))
}

sce_add_h5_to_graphs <- function(sce, h5, cellNames, graphsName = "graphs") {
  for (name in names(h5[[graphsName]])) {
    tryCatch(
      {
        Graphmt <- read_matrice_matrix(h5[[graphsName]][[name]])
        rownames(Graphmt) <- cellNames
        colnames(Graphmt) <- cellNames
        sce@graphs[[name]] <- Seurat::as.Graph(Graphmt)
      },
      error = function(e) {
        print(paste0("Reading graph ", name, ": ", e$message, "failed."))
      }
    )
  }
  return(sce)
}

sce_add_h5_to_var <- function(sce, h5, assay, varName = "var", SeuratVersion = checkSeuratVersion()) {
  print("Starting sce_add_h5_to_var function...")
  
  # Debug input arguments
  print(paste("Seurat version:", SeuratVersion))
  print(paste("Assay name:", assay))
  print(paste("Variable name in HDF5:", varName))
  
  # Check if h5[varName] exists
  if (!varName %in% names(h5)) {
    stop(paste("Variable", varName, "not found in HDF5 object."))
  }
  
  # Check if rawvar exists within varName
  if (!"rawvar" %in% names(h5[[varName]])) {
    stop(paste("rawvar not found within", varName, "in HDF5 object."))
  }
  
  print("Attempting to convert rawvar to a data frame...")
  varData <- tryCatch({
    h5_to_DF(h5[[varName]][["rawvar"]])
  }, error = function(e) {
    print("Error converting rawvar to data frame:")
    print(e)
    return(NULL)
  })
  
  if (is.null(varData)) {
    stop("Failed to convert rawvar to a data frame. Exiting function.")
  }
  
  print("Successfully converted rawvar to data frame.")
  print(paste("Dimensions of varData:", paste(dim(varData), collapse = " x ")))

  # Add meta data based on Seurat version
  if (SeuratVersion == 5) {
    print("Adding metadata to Seurat object (version 5)...")
    sce@assays[[assay]] <- tryCatch({
      Seurat::AddMetaData(sce@assays[[assay]], varData)
    }, error = function(e) {
      print("Error adding metadata:")
      print(e)
      stop("Failed to add metadata in Seurat version 5.")
    })
  } else if (SeuratVersion == 4) {
    print("Adding metadata to Seurat object (version 4)...")
    sce@assays[[assay]]@meta.features <- tryCatch({
      varData
    }, error = function(e) {
      print("Error adding metadata:")
      print(e)
      stop("Failed to add metadata in Seurat version 4.")
    })
  } else {
    stop("Unsupported Seurat version. Supported versions are 4 and 5.")
  }
  
  print("Metadata successfully added to Seurat object.")
  return(sce)
}


sce_add_h5_to_reductions <- function(sce, h5, cellNames, assay = "RNA", reductionsName = "reductions") {
  name <- names(h5[[reductionsName]])[[1]]
  for (name in names(h5[[reductionsName]])) {
    tryCatch(
      {
        matrix <- h5[[reductionsName]][[name]][, ]
        colnames(matrix) <- cellNames
        sce@reductions[[name]] <- Seurat::CreateDimReducObject(
          embeddings = t(matrix),
          key = paste0(name, "_"),
          assay = assay
        )
      },
      error = function(e) {
        print(paste0("Reading reductions ", name, ": ", e$message, "failed."))
      }
    )
  }
  return(sce)
}

#' loadH5
#'
#' @param FileName
#' @param assay
#' @param SeuratVersion
#' @param calData
#' @param calScale
#' @param calFeatures
#' @param readType
#'
#' @return
#' @export
#'
#' @examples
loadH5 <- function(FileName,
                   assay = "RNA",
                   SeuratVersion = checkSeuratVersion(),
                   image_name = "Spatial",
                   useBPcells = FALSE,
                   useDisk = TRUE,
                   calData = TRUE,
                   calScale = FALSE,
                   calFeatures = FALSE,
                   group_by = NULL,
                   readType = "Seurat") {
  options(warn = -1)
  if (!readType %in% c("Seurat", "monocle2", "monocle3", "cellchat", "SingleCellExperiment")) {
    stop("readType should be one of Seurat, monocle2, monocle3")
  }

  if (readType == "monocle2") {
    if (!require(monocle)) {
      stop("The monocle package is not installed, please install it first.")
    }
  }

  if (readType == "monocle3") {
    if (!require(monocle3)) {
      stop("The monocle3 package is not installed, please install it first.")
    }
  }

  if (readType == "cellchat") {
    if (!require(CellChat)) {
      stop("The cellchat package is not installed, please install it first.")
    }
    if (is.null(group_by)) {
      stop("group_by should be specified for cellchat data.")
    }
  }

  if (readType == "SingleCellExperiment") {
    if (!require(SingleCellExperiment)) {
      stop("The SingleCellExperiment package is not installed, please install it first.")
    }
  }

  h5 <- openH5(FileName)
  tryCatch(
    {
      sce <- h5_to_seurat(h5, assay, SeuratVersion,
        image_name = image_name,
        useBPcells = useBPcells, useDisk = useDisk, calData = calData, calScale = calScale, calFeatures = calFeatures
      )
    },
    error = function(e) {
      print(e)
    },
    finally = {
      h5$close_all()
    }
  )

  if (readType == "Seurat") {
    return(sce)
  } else if (readType == "monocle2") {
    library(monocle)
    cds <- Seurat_to_Monocle2(sce, assay = assay, SeuratVersion = SeuratVersion)
    return(cds)
  } else if (readType == "monocle3") {
    library(monocle3)
    cds <- Seurat_to_Monocle3(sce, assay = assay, SeuratVersion = SeuratVersion)
    return(cds)
  } else if (readType == "cellchat") {
    loadlayers <- "data"
    if (SeuratVersion == 5) {
      if (!("data" %in% names(sce@assays[[assay]]@layers)) & !calData) {
        print("No data found in the h5 file and calData is FALSE, using the counts in the SeuratObject.")
        loadlayers <- "counts"
      }
    }
    cellchatdata <- Seurat_to_cellchat(sce, assay = assay, SeuratVersion = SeuratVersion, group_by = group_by, loadlayers = loadlayers)
    return(cellchatdata)
  } else if (readType == "SingleCellExperiment") {
    SingleCellData <- Seurat::as.SingleCellExperiment(sce)
    return(SingleCellData)
  } else {
    return(sce)
  }
}

Seurat_to_cellchat <- function(
    sce,
    assay = "RNA",
    SeuratVersion = checkSeuratVersion(),
    group_by = "orig.ident",
    loadlayers = "data") {
  if (SeuratVersion == 4) {
    data <- SeuratObject::GetAssayData(object = sce, assay = assay, slot = loadlayers)
  } else if (SeuratVersion == 5) {
    data <- SeuratObject::GetAssayData(object = sce, assay = assay, layer = loadlayers)
  }
  meta <- sce@meta.data
  cellchat <- createCellChat(object = data, meta = meta, group.by = group_by)
  return(cellchat)
}

Seurat_to_Monocle3 <- function(
    sce,
    assay = "RNA",
    SeuratVersion = checkSeuratVersion()) {
  if (SeuratVersion == 4) {
    data <- SeuratObject::GetAssayData(object = sce, assay = assay, slot = "counts")
  } else if (SeuratVersion == 5) {
    data <- SeuratObject::GetAssayData(object = sce, assay = assay, layer = "counts")
  }

  cell_metadata <- sce@meta.data
  gene_annotation <- data.frame(gene_short_name = rownames(sce))
  rownames(gene_annotation) <- rownames(sce)

  cds <- monocle3::new_cell_data_set(data,
    cell_metadata = cell_metadata,
    gene_metadata = gene_annotation
  )
  return(cds)
}

Seurat_to_Monocle2 <- function(
    sce,
    assay = "RNA",
    SeuratVersion = checkSeuratVersion()) {
  if (SeuratVersion == 4) {
    data <- SeuratObject::GetAssayData(object = sce, assay = assay, slot = "counts")
  } else if (SeuratVersion == 5) {
    data <- SeuratObject::GetAssayData(object = sce, assay = assay, layer = "counts")
  }

  data <- as(as.matrix(data), "sparseMatrix")
  pd <- new("AnnotatedDataFrame", data = sce@meta.data)
  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  fd <- new("AnnotatedDataFrame", data = fData)

  cds <- monocle::newCellDataSet(data,
    phenoData = pd,
    featureData = fd,
    lowerDetectionLimit = 0.1,
    expressionFamily = negbinomial.size()
  )
  return(cds)
}

h5_to_seurat <- function(h5,
                         assay = "RNA",
                         SeuratVersion = checkSeuratVersion(),
                         image_name = "Spatial",
                         useBPcells = FALSE,
                         useDisk = TRUE,
                         calData = TRUE,
                         calScale = FALSE,
                         calFeatures = FALSE) {
  # Retrieve expected cell and gene names from the h5 file
  cellNames <- h5[["names_obs"]][]
  geneNames <- h5[["names_var"]][]
  allgeneNames <- h5[["var/rawvar/_index"]][]

  print(sprintf("Assay in file: %s", names(h5[["assay"]])))
  print("Reading raw data...")

  # Read raw data; note that rawX should have rows corresponding to genes and columns to cells
  rawX <- h5_to_X(h5, assay, layer = "rawdata", useBPcells = useBPcells, useDisk = useDisk, cellNames = cellNames, geneNames = allgeneNames)
  print(sprintf("Raw data dimensions: %d rows x %d columns", nrow(rawX), ncol(rawX)))
  print(sprintf("Number of cell names: %d", length(cellNames)))
  print(sprintf("Number of gene names (from rawvar index): %d", length(allgeneNames)))
  
  # Check for consistency in dimensions
  if (ncol(rawX) != length(cellNames)) {
    stop(sprintf("Column mismatch: raw data has %d columns but cellNames has length %d.", ncol(rawX), length(cellNames)))
  }
  if (nrow(rawX) != length(allgeneNames)) {
    stop(sprintf("Row mismatch: raw data has %d rows but gene names has length %d.", nrow(rawX), length(allgeneNames)))
  }
  
  # Create the Seurat object using the raw counts
  sce <- Seurat::CreateSeuratObject(counts = rawX, assay = assay, project = "scanpy")
  print("Seurat object created:")
  print(sce)
  
  # Normalize counts if requested
  print("Normalization counts...")
  if (calData) {
    sce <- Seurat::NormalizeData(sce)
  }
  
  # If you want to add features from the h5 file (or calculate new features)
  print("Add Feature...")
  if (calFeatures) {
    sce <- Seurat::FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
  } else if ("data" %in% names(h5[["assay"]][[assay]][["layers"]])) {
    print("Loading features from h5 file")
    Seurat::VariableFeatures(sce) <- rownames(h5_to_DF(h5[["var"]][["var"]]))
  } else {
    print("No feature information available in h5 file; please set calFeatures = TRUE to calculate features.")
  }
  
  # (Other steps such as scale data, adding meta data, obs, reductions, etc., go here)
  print("Add obs...")
  sce <- Seurat::AddMetaData(sce, h5_to_obs(h5, "obs"))
  
  # … Further processing of graphs, reductions, uns, images, etc.
  
  return(sce)
}

h5_to_images <- function(images, assay, image_name = NULL, cellNames = NULL) {
  library(Seurat)
  coords <- as.data.frame(images[["coords"]][, ])
  if (!is.null(cellNames)) {
    rownames(coords) <- cellNames
  }
  colnames(coords) <- c("imagerow", "imagecol")
  image <- images[["image"]][, , ]
  scale_factors <- images[["scale_factors"]]
  fiducial <- scale_factors[["fiducial"]][]
  hires <- scale_factors[["hires"]][]
  lowres <- scale_factors[["lowres"]][]
  spot <- scale_factors[["spot"]][]

  scale.factors <- scalefactors(
    spot = spot,
    fiducial = fiducial,
    hires = hires,
    lowres = lowres
  )
  fov <- CreateFOV(as.data.frame(coords),
    type = "centroids", radius = spot,
    assay = assay, key = Key(image_name, quiet = TRUE)
  )
  visium.fov <- new(
    Class = "VisiumV2", boundaries = fov@boundaries,
    molecules = fov@molecules, assay = fov@assay, key = fov@key,
    image = image, scale.factors = scale.factors
  )
  # visium.fov@project.name <- image_name
  return(visium.fov)
}
