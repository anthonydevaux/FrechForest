#' Function to compute individual random effects using hlme output parameters
#'
#' @param model output object from hlme function
#' @param data data to compute random effect
#'
#' @return
#' @export
#'
#' @importFrom lcmm fixef
#'
#' @examples
predRE <- function(model, data){

  # a revoir pour gerer les NA

  subject <- model$call$subject
  beta <- lcmm::fixef(model)[[2]]

  # Variance-covariance matrix of the random-effects

  B <- matrix(0, ncol = sum(model$idea0), nrow = sum(model$idea0))
  colnames(B) <- model$Xnames[model$idea0==1]
  rownames(B) <- model$Xnames[model$idea0==1]
  B[upper.tri(B,diag=TRUE)] <- model$best[(model$N[1]+model$N[2]+1):(model$N[1]+model$N[2]+model$N[3])]
  B <- t(B)
  B[upper.tri(B,diag=TRUE)] <- model$best[(model$N[1]+model$N[2]+1):(model$N[1]+model$N[2]+model$N[3])]

  se <- tail(model$best, n = 1)^2 # residual variance error
  Z <- model.matrix(as.formula(model$call$random), data) # random design matrix
  X <- model.matrix(as.formula(model$call$fixed[-2]), data) # fixed design matrix
  Y <- na.omit(data[,as.character(model$call$fixed)[2], drop = FALSE]) # outcome

  bi <- matrix(NA, nrow = length(unique(data$id)), ncol = ncol(B),
               dimnames = list(unique(data$id), colnames(Z)))

  for (id in unique(data$id)){

    row.id <- which(data$id==id)
    Zi <- Z[row.id, , drop = FALSE]
    Xi <- X[row.id, ]
    Yi <- Y[row.id, ]
    Vi <- Zi%*%B%*%t(Zi) + se*diag(length(row.id))
    b <- B%*%t(Zi)%*%solve(Vi)%*%(Yi-Xi%*%beta)

    bi[rownames(bi)==id,] <- b

  }

  return(list(bi = bi))

}
