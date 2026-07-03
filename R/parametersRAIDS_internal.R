


validateParamRAIDS <- function(parameters) {
  pRAIDSNames <- c("reference", # HGDP1kg
                      "studyDF",
                      "studyDFSyn",
                      "pedStudy",
                      "studyType",
                      "genoSource",
                      "blockTypeId",
                      "genome",
                      "chrInfo",
                      "profileFile",
                      "profileFileGeno",
                      "pathProfileGDS",
                      "fileReferenceGDS",
                      "fileReferenceAnnotGDS",
                      "inferenceType",
                      "batch",
                      "prefix",
                      "nbSim",
                      "offset",
                      "minCov",
                      "minProb",
                      "seqError",
                      "seqErrorSyn",
                      "pRecomb",
                      "np",
                      "listPos",
                      "syntheticRefDF",
                      "pruningMethod",
                      "slideWindowMaxBP",
                      "thresholdLD",
                      "specificSNV",
                      "genoType",
                      "phaseType",
                      "phase",
                      "PCAmissingRate",
                      "PCAalgorithm",
                      "eigenCount",
                      "eigenCountSyn",
                      "kList",
                      "pcaList",
                      "fieldPopInRef",
                      "fieldPopInfAnc",
                      "fieldSubPop",
                      "verbose")
  if (!is.numeric(params$n_iter) || length(params$n_iter) != 1 ||
      params$n_iter <= 0 || params$n_iter %% 1 != 0) {
    stop("`n_iter` must be a single positive integer.")
  }

  if (!is.numeric(params$alpha) || length(params$alpha) != 1 ||
      params$alpha <= 0 || params$alpha >= 1) {
    stop("`alpha` must be a single number between 0 and 1.")
  }

  allowed_models <- c("linear", "logistic", "poisson")
  if (!is.character(params$model_type) || length(params$model_type) != 1 ||
      !params$model_type %in% allowed_models) {
    stop(sprintf("`model_type` must be one of: %s.",
                  paste(allowed_models, collapse = ", ")))
  }

  if (!is.character(params$covariates) || length(params$covariates) == 0) {
    stop("`covariates` must be a non-empty character vector.")
  }

  if (!is.character(params$output_dir) || length(params$output_dir) != 1) {
    stop("`output_dir` must be a single string.")
  }

  invisible(TRUE)
}



