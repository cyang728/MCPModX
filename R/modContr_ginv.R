#' Calculate Optimal Contrasts for Multiple Models Using Generalized Inverse
#'
#' This function calculates the optimal contrasts for multiple models.
#'
#' @param means A matrix of means.
#' @param W A weight matrix. Default is NULL.
#' @param Sinv A matrix representing the generalized inverse of the covariance matrix. Default is NULL.
#' @param placAdj Logical, whether to adjust for placebo (default is FALSE).
#' @param type Character, either "unconstrained" or "constrained".
#' @param direction Character, direction of the effect ("increasing" or "decreasing").
#'
#' @return A matrix of optimal contrasts.
#' @export
#'
modContr_ginv = function(means, W = NULL, Sinv = NULL, placAdj = FALSE,
                         type, direction){
  ## call optC on matrix
  ## check whether constant shape was specified and remove (can happen for linInt model)
  if(!placAdj){
    ind = apply(means, 2, function(x){
      length(unique(x)) > 1
    })
  } else { ## placAdj
    ind = apply(means, 2, function(x){
      any(x != 0)
    })
  }
  if(all(!ind))
    stop("All models correspond to a constant shape, no optimal contrasts calculated.")
  if(any(!ind)){
    nam = colnames(means)[!ind]
    namsC = paste(nam, collapse = ", ")
    if(length(nam) == 1){
      message("The ", namsC, " model has a constant shape, cannot
calculate optimal contrasts for this shape.")
    } else {
      message("The ", namsC, " models have a constant shape, cannot
calculate optimal contrasts for these shapes.")
    }
    means = means[,ind, drop=FALSE]
  }

  if(is.null(Sinv))
    Sinv = ginv(W)
  if(type == "unconstrained"){
    out = apply(means, 2, optC_ginv, Sinv = Sinv, placAdj = placAdj)
  } else { # type == "constrained"
    out = apply(means, 2, constOptC_ginv, Sinv = Sinv,
                placAdj = placAdj, direction = direction)
  }
  if(!is.matrix(out)){ ## can happen for placAdj=T and only 1 act dose
    nam = names(out)
    out = matrix(out, nrow = 1)
    colnames(out) = nam
  }
  out
}
