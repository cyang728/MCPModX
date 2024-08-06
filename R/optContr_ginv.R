#' Calculate Optimal Contrasts and Critical Value Using Generalized Inverse
#'
#' This function calculates optimal contrasts and critical values for specified models and doses.
#'
#' @param models An object of class `Mods` representing the dose-response models.
#' @param doses A numeric vector of dose levels. Default is taken from `models`.
#' @param w A weight vector. Default is NULL.
#' @param S A covariance matrix. Default is NULL.
#' @param placAdj Logical, whether to adjust for placebo (default is FALSE).
#' @param type Character, either "unconstrained" or "constrained". Default is "unconstrained".
#'
#' @return A list containing the contrast matrix, mean response matrix, and correlation matrix.
#' @export
#'
optContr_ginv = function(models, doses, w, S, placAdj = FALSE,
                         type = c("unconstrained", "constrained")){
  ## calculate optimal contrasts and critical value
  if(!(inherits(models, "Mods")))
    stop("models needs to be of class Mods")
  if(missing(doses))
    doses = attr(models, "doses")
  scal = attr(models, "scal")
  off = attr(models, "off")
  nodes = attr(models, "doses")
  direction = unique(attr(models, "direction"))
  if(length(direction) > 1)
    stop("need to provide either \"increasing\" or \"decreasing\" as direction to optContr")
  mu = getResp(models, doses)
  if(placAdj){
    mu0 = getResp(models, 0)
    mu = mu-matrix(mu0[1,], byrow = TRUE,
                   nrow=nrow(mu), ncol=ncol(mu))
  }
  type = match.arg(type)
  if(type == "constrained"){
    avail = requireNamespace("quadprog", quietly = TRUE)
    if(!avail)
      stop("Need suggested package quadprog to calculate constrained contrasts")
  }
  if(any(doses == 0) & placAdj)
    stop("If placAdj == TRUE there should be no placebo group in \"doses\"")
  ## check for n and vCov arguments
  if(!xor(missing(w), missing(S)))
    stop("Need to specify exactly one of \"w\" or \"S\"")
  if(!missing(w)){
    if(length(w) == 1){ # assume equal weights
      S = Sinv = diag(length(doses))
    } else {
      if(length(w) != length(doses))
        stop("w needs to be of length 1 or of the same length as doses")
      S = diag(1/w)
      Sinv = diag(w)
    }
  } else {
    if(!is.matrix(S))
      stop("S needs to be a matrix")
    Sinv = ginv(S)
  }
  contMat = modContr_ginv(mu, Sinv=Sinv, placAdj = placAdj,
                          type = type, direction = direction)
  rownames(contMat) = doses
  corMat = cov2cor(t(contMat) %*% S %*% contMat)
  res = list(contMat = contMat, muMat = mu, corMat = corMat)
  attr(res, "type") = type
  attr(res, "placAdj") = placAdj
  class(res) = "optContr"
  res
}
