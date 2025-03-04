#-------------------------------------------------------------------------------
#' Features selection based on Mutual Information (MI)
#'
#' @description Select variables by computing mutual information between the
#' variables and a variable of interest.
#'
#' @param input_data [a data frame or a matrix, required]
#'
#' Expected layout is samples as rows and variables as columns.
#'
#' @param var_of_interest_name [a string, optional, NULL by default]
#'
#' The name of the variable of interest. If not NULL, the variable name
#' must be one the variables of \emph{input_data}.
#'
#' @param var_of_interest_values [a vector, NULL by default]
#'
#' The values for the variable of interest. When the the variable of interest
#' is not in the input data (when \emph{var_of_interest_name} is set to NULL),
#' a list of values have to be supplied for the variable of interest.
#' The length of the vector must be equal to the number of samples in
#' \emph{input_data}.
#'
#' @param n_features [an integer, optional, NULL by default]
#'
#' The number of features to select. If NULL, no filter is applied and
#' all the input variables are returned with their mutual information values.
#'
#' @param skip_cheks [a boolean, FALSE by default]
#'
#' Before computing MI between the variable of interest and the features,
#' \emph{input_data} is checked to filter out constant features and rows full
#' of NAs. When the \emph{input_data} does not need such filtering,
#' these checks can be skipped to speed up the process.
#'
#' @param n_threads [a positive integer, optional, 1 by default]
#'
#' When set greater than 1, n_threads parallel threads will be used for
#' computation. Make sure your compiler is compatible with openmp
#' if you wish to use multithreading.
#'
#' @param verbose [a boolean, TRUE by default]
#'
#' TRUE to show the progress, FALSE to turn off the display.
#'
#' @return A named vector with features as names and MI in bits as values.
#' The vector is sorted in decreasing order and, if requested, filtered on
#' the top \emph{n_features} variables.
#'
#' @export
#-------------------------------------------------------------------------------
selectFeatures_v2 <- function (input_data, var_of_interest_name=NULL,
                               var_of_interest_values=NULL, n_features=NULL,
                               skip_cheks=F, n_threads=1, verbose=T)
{
  # LN_2 equivalent to constant for the function
  #
  LN_2 <- log(2)
  #
  # Check inputs
  #
  if (  ( ! is.data.frame (input_data) )
        && ( ! is.matrix(input_data) )
        && ( ! inherits(input_data, "Matrix") ) )
    miic:::miic_error  ("parameters", "the input data must be a data frame or a matrix.")
  if ( (ncol (input_data) <= 0) || (nrow (input_data) <= 0) )
    miic:::miic_error  ("parameters", "the input data is empty.")
  if ( is.data.frame (input_data) )
    # Ensure we have a true data frame, e.g. not a tibble
    input_data <- as.data.frame (input_data)
  
  if ( is.null (var_of_interest_name) && is.null (var_of_interest_values) )
    miic:::miic_error  ("parameters", "the name of the variable of interest",
                        " or a list of values must be supplied.")
  if ( ! is.null (var_of_interest_name) )
  {
    # if ( miic:::test_param_wrong_string (var_of_interest_name, colnames(input_data) ) )
    #   stop ("The variable of interest name is incorrect.")
    #
    # Extract variable of interest from the input data and remove it
    #
    # var_of_interest_values <- input_data[, var_of_interest_name]
    # input_data <- input_data[, colnames(input_data) != var_of_interest_name,
    #                          drop=F]
    # if (ncol (input_data) <= 0)
    #   miic:::miic_error  ("parameters", "the input data contain only the variable of interest.")
  }
  else
  {
    if ( (  ( ! is.factor(var_of_interest_values) )
            && ( ! is.vector(var_of_interest_values) ) )
         || ( length (var_of_interest_values) != nrow (input_data) ) )
      miic:::miic_error  ("parameters",
                          "the variable of interest values are incorrect.")
  }
  
  if (  ( ! is.null (n_features) )
        && miic:::test_param_wrong_int ( n_features, 0, ncol(input_data) ) )
  {
    miic:::miic_warning ("parameters",
                         "the number of features to select is incorect and will be ignored.")
    n_features <- NULL
  }
  
  skip_cheks <- miic:::check_param_logical (skip_cheks, "skip checks", F)
  n_threads <- miic:::check_param_int (n_threads, "number of threads", 1)
  verbose <- miic:::check_param_logical (verbose, "verbose", T)
  #
  # The bin size controls the number of features evaluated in one go
  #
  if ( (ncol (input_data) < 750) || (nrow(input_data) <= 2000) )
    bin_size = 100
  else  if (nrow(input_data) <= 4000)
    bin_size = 75
  else  if (nrow(input_data) <= 7000)
    bin_size = 50
  else
    bin_size = 20
  #
  # Prepare the returned value: a vector with all features set with 0 MI
  #
  n_vars = ncol (input_data)
  v_mis = rep (0, n_vars)
  names (v_mis) = colnames (input_data)
  #
  # Compute the mutual information by group of bin_size variables using miic
  #
  start_idx <- 1
  while (start_idx < n_vars)
  {
    end_idx <- min (start_idx + bin_size - 1, n_vars)
    # print(paste0 ("From ", start_idx, " to ", end_idx, " (n_vars=", n_vars, ")") )
    if (verbose)
      cat (paste0 ("Computing MI: ",
                   format (round ( ((start_idx-1) / n_vars) * 100, 2), nsmall=2), " %   \r") )
    
    data_loop <- input_data [, start_idx:end_idx, drop=FALSE]
    if ( ! is.data.frame(data_loop) )
      data_loop = as.data.frame (data_loop)
    
    if (!skip_cheks)
    {
      # Remove rows full of NAs and constant variables
      # (would generate warnings if sent to miic function)
      #
      count_vals = unlist (apply (data_loop, MARGIN=2, FUN=function (x) {
        length (unique (x[!is.na(x)] ) ) }) )
      data_loop = data_loop[, count_vals >= 2, drop=F]
      
      count_nas = apply (data_loop, MARGIN=1, FUN=function(x) { sum (is.na(x) ) } )
      data_loop = data_loop[ count_nas < ncol(data_loop), , drop=F]
    }
    
    # print (paste0 ("nrow: ", nrow (data_loop),
    #               ", ncol: ", ncol (data_loop) ) )
    #
    if ( (nrow (data_loop) > 0) && (ncol (data_loop) > 0) )
    {
      so <- data.frame ("var_names" = c (colnames(data_loop), "var_interest"),
                        "is_consequence" = c (rep(1, ncol(data_loop)), 0),
                        stringsAsFactors = FALSE)
      data_loop$var_interest <- var_of_interest_values
      miic_res <- miic (data_loop, state_order=so,
                        orientation=F, latent="no", n_threads=n_threads, verbose=0)
      miic_res <- miic_res$summary
      rownames (miic_res) <- NULL
      rownames (miic_res)[miic_res$x != "var_interest"] <- (
        miic_res[miic_res$x != "var_interest", "x"] )
      rownames (miic_res)[miic_res$y != "var_interest"] <- (
        miic_res[miic_res$y != "var_interest", "y"] )
      v_mis[rownames(miic_res)] = round (
        (miic_res$info_shifted / miic_res$n_xy_ai) / LN_2, 6)
    }
    start_idx <- start_idx + bin_size
  }
  
  if (verbose)
    cat ("Computing MI: 100 %   \r")
  
  v_mis <- sort (v_mis, decreasing=TRUE)
  if ( is.null (n_features) )
  {
    if (verbose)
      cat (paste0 (length(v_mis), " features evaluated   \n") )
  }
  else
  {
    if (verbose)
      cat (paste0 (length(v_mis), " features evaluated, top ",
                   min( n_features:length(v_mis) ), " returned\n") )
    v_mis <- v_mis[1:min( n_features:length(v_mis) )]
  }
  return (v_mis)
}