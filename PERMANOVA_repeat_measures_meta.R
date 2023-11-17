PERMANOVA_repeat_measures_meta <- function(
  D,
  study, study_longitudinal,
  subject, subject_data = NULL,
  sample_data = NULL,
  metadata_order = c(names(subject_data), names(sample_data)),
  permutations=999, ncores=1)
{
  
  # Make sure D is a dist object
  if (class(D) != "dist") {
    stop("D must be a dist object")
  }
  
  # check sample identifiers in D, study, and subject
  if(nrow(as.matrix(D)) != length(study) |
     nrow(as.matrix(D)) != length(subject))
    stop("sample number from D, study, and subject must match!")
  if(any(rownames(as.matrix(D)) != names(study)) |
     any(rownames(as.matrix(D)) != names(subject)))
    stop("sample identifiers from D, study, and subject must match!")
  
  # check study data
  if(!is.character(study))
    stop("study identifiers must be character!")
  if(!is.logical(study_longitudinal))
    stop("study_longitudinal must be TRUE/FALSE!")
  if(length(unique(study)) != length(study_longitudinal))
    stop("study number from study and study_longitudinal must match!")
  if(!setequal(study, names(study_longitudinal)))
    stop("study identifiers from study and study_longitudinal must match!")
  
  # subject_data and sample_data cannot both be missing
  if(is.null(subject_data) & is.null(sample_data)) {
    stop("At least one of subject_data and sample_data must be provided!")
  }
  
  # check subject data, create if missing
  if(!is.character(subject))
    stop("subject must be character!")
  if(is.null(subject_data)) {
    subject_data <- data.frame(place_holder_subject = rep(1, length(unique(subject))),
                               row.names = unique(subject))
  }
  if(length(unique(subject)) != nrow(subject_data))
    stop("subject number from subject and subject_data must match!")
  if(!setequal(subject, rownames(subject_data)))
    stop("subject identifiers from subject and subject_data must match!")
  
  # study identifiers and subject identifiers cannot overlap for the funciton to work
  if(length(intersect(study, subject)) > 0)
    stop("study identifiers and subject identifiers cannot overlap!")
  # check that subject are study-specific
  if(any(apply(table(study, subject) > 0, 2, sum) > 1))
    stop("subject must be study-specific!")
  # map studies to subjects, and create block variable to permute
  # subject-specific data with
  study_subject <- study[!duplicated(subject)]
  names(study_subject) <- subject[!duplicated(subject)]
  study_permute <- study_subject[rownames(subject_data)]
  
  # create sample data if missing
  if(is.null(sample_data)) {
    sample_data <- data.frame(place_holder_sample = rep(1, nrow(as.matrix(D))),
                              row.names = rownames(as.matrix(D)))
  }
  if(nrow(as.matrix(D)) != nrow(sample_data))
    stop("sample number from D and sample_data must match!")
  if(any(rownames(as.matrix(D)) != rownames(sample_data)))
    stop("sample identifiers from D and sample_data must match!")
  # create subject identifiers to permute sample-specific data with
  # these are just subject identifiers for longitudinal data
  # and study identifiers for cross-sectional data
  subject_permute <- ifelse(study_longitudinal[study],
                            subject,
                            study)
  names(subject_permute) <- names(study)
  
  # Ensure no metadata overlap between subject_data and sample_data
  if(length(intersect(names(subject_data), names(sample_data))) > 0) {
    stop("metadata is repeated across subject_data and sample_data")
  }
  
  # Ensure that metadata_order only contains stuff in subject_data and sample_data
  if(length(setdiff(metadata_order, union(names(subject_data), names(sample_data)))) > 0) {
    stop("metadata_order contains metadata not in subject_data and sample_data!")
  }
  
  
  # Warn on some suspicious input
  # persample <- apply(sample_data, 1, function(x)is.factor(x) && !any(duplicated(x)))
  # if (any(persample)) {
  #   warning(sprintf("%s in sample_data has one DOF per sample.", colnames(sample_data)[which(persample)[1]]))
  # }
  # if (length(unique(subject)) < nrow(subject_data)) {
  #   warning("Not all subject have a sample associated with them. Block permutations will still be performed over the full set of subject - if this is not desired, subset subject_data to only the subject which appear in the data.")
  # }
  # if (!any(duplicated(subject))) {
  #   warning("subject contains no duplicated elements")
  # }
  library(permute)
  mtdat <- cbind(subject_data[subject,,drop=F], sample_data)
  # Error out on NA metadata rather than allowing adonis to error out with a totally
  # nonsensical error message
  if (any(is.na(mtdat[, metadata_order, drop=F]))) {
    stop("Some metadata is NA! adonis does not support any NA in the metadata")
  }
  # Test statistic from non-permuted data
  ad <- vegan::adonis(D ~ ., permutations=0, data=mtdat[, metadata_order, drop=F])
  R2 <- ad$aov.tab$R2
  names(R2) <- rownames(ad$aov.tab)
  
  doParallel::registerDoParallel(ncores)
  nullsamples <- foreach::`%dopar%`(
    foreach::foreach(i = seq_len(permutations), .combine = cbind),
    {
      subject.i <- shuffle(nrow(subject_data), 
                           control = how(blocks = study_permute))
      sample.i <- shuffle(nrow(sample_data), 
                          control = how(blocks = subject_permute))
      i.subject_data <- subject_data[subject.i,,drop=F]
      rownames(i.subject_data) <- rownames(subject_data)
      mtdat <- cbind(i.subject_data[subject,,drop=F],
                     sample_data[sample.i,,drop=F])
      perm.ad <- vegan::adonis(D ~ ., permutations=0, data=mtdat[, metadata_order, drop=F])
      return(perm.ad$aov.tab$R2)
    })
  doParallel::stopImplicitCluster()
  
  
  # For residuals, test the other direction (i.e. p-value of all covariates)
  n <- length(R2)
  R2[n-1] <- 1 - R2[n-1]
  nullsamples[n-1,] <- 1 - nullsamples[n-1,]
  
  # P value calculation similar to adonis's
  exceedances <- rowSums(nullsamples > R2)
  P <- (exceedances + 1) / (permutations + 1)
  
  P[n] <- NA    # No p-values for "Total"
  ad$aov.tab$`Pr(>F)` <- P
  
  return (ad)
}
