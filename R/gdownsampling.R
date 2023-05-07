#' Gene Expression Down-sampling for Matrix
#'
#' Down-sample gene expression based on matrix
#' @param mat The input matrix that will be down-sampled
#' @param mean The mean parameter for down-sampling
#' @param sd The standard deviation parameter for down-sampling
#' @return The matrix after down-sampling
#' @examples
#' library(scRNAseq)
#' sce = ReprocessedAllenData("tophat_counts")
#' counts = assay(sce, 'tophat_counts')
#' d_counts = gdownsample_mat(counts, mean=0.3, sd=0.05)
#' @export
gdownsample_mat = function(mat, mean, sd){
  beta = mean/sd^2
  alpha = beta*mean
  shrinkages = rgamma(ncol(mat), alpha, beta)
  mat2 = matrix(0, ncol = ncol(mat), nrow = nrow(mat))
  for (i in 1:ncol(mat)){
    libSize2 = sum(mat[, i])*shrinkages[i]
    prop = exp(log(mat[, i]/sum(mat[, i])))
    mat2[, i] = rmultinom(n=1,size=libSize2,prob=prop)
    cat("Now is ", i,"\n")
  }
  rownames(mat2) = rownames(mat)
  colnames(mat2) = colnames(mat)
  return(mat2)
}
#' Gene Expression Down-sampling for Seurat or SingleCellExperiment Object
#'
#' Down-sample gene expression based on Seurat or SingleCellExperiment object
#' @param object The input Seurat or SingleCellExperiment object for down-sampling
#' @param user_assay The input assay for down-sampling. For Seurat object, if not specified, down-sampling will be performed based on the active assay. For SingleCellExperiment object, an assay has to be specified.
#' @param returnMode The return mode for down-sampling. Default is "new": creating a new down-sampled object while keep the original object. Other two available options are 1. "add": add a new down-sampled assay. 2. "replace": add a new down-sampled assay and delete the original assay.
#' @param mean The mean parameter for down-sampling
#' @param sd The standard deviation parameter for down-sampling
#' @return A Seurat/SingleCellExperiment object or none based on the return mode parameter.
#' @examples
#' library(scRNAseq)
#' sce = ReprocessedAllenData("tophat_counts")
#' d_sce = gdownsample_obj(sce,user_assay="tophat_counts", returnMode="new", mean=0.3, sd=0.05)
#' @export
gdownsample_obj = function(object, user_assay=NULL, returnMode="new",
                        mean, sd){
  if(class(object)=="Seurat"){
    # if user specify an assay and it exists in Seurat object
    if(!is.null(user_assay)&&user_assay%in%names(object@assays)){
      assay_name = user_assay
      assay = object[[user_assay]]
      mat = assay@counts
    }else{
      message("No input assay name. Use activate assay for Seurat object instead.")
      assay_name = object@active.assay
      assay = object[[assay_name]]
      mat = assay@counts
    }

    ## replace mat here
    mat2 = gdownsample_mat(mat, mean, sd)

    if(returnMode=="add"){
      new_assay = assay
      new_assay@counts = mat2
      new_assay@key = "downsampled_"
      object[['downsampled']] = new_assay
      object@active.assay = 'downsampled'
    } else if(returnMode=="replace"){
      new_assay = assay
      new_assay@counts = mat2
      new_assay@key = "downsampled_"
      object[['downsampled']] = new_assay
      object[[assay_name]] = NULL
      object@active.assay = 'downsampled'
    } else{
      new_object = object
      new_assay = assay
      new_assay@counts = mat2
      new_assay@key = "downsampled_"
      new_object[['downsampled']] = new_assay
      new_object@active.assay = 'downsampled'
      new_object[[assay_name]] = NULL
      return(new_object)
    }
  }

  if(class(object)=="SingleCellExperiment"){
    # if user specify an assay and it exists in SCE object
    if(!is.null(user_assay)&&user_assay%in%assayNames(object)){
      mat = assay(object, user_assay)
    }
    # else stop function
    else{
      stop("You must input an existing assay for SingleCellExperiment object!")
    }

    mat2 = gdownsample_mat(mat, mean, sd)

    if(returnMode=="add"){
      assay(object, "downsampled") = mat2
    } else if(returnMode=="replace"){
      assay(object, "downsampled") = mat2
      assay(object, user_assay) = NULL
    } else{
      new_object = object
      assay(new_object, "downsampled") = mat2
      assay(new_object, user_assay) = NULL
      return(new_object)
    }
  }
}
