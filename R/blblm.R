#' @import purrr
#' @import stats
#' @import parallel
#' @importFrom magrittr %>%
#' @importFrom utils count.fields
#' @importFrom utils capture.output
#' @importFrom utils read.csv
#' @details
#' Linear Regression with Little Bag of Bootstraps
"_PACKAGE"


## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
utils::globalVariables(c("."))


#' @param formula the user specifies the independent and dependent variables in their analysis
#' @param data the data frame or series of csv files with data
#' @param m an integer indicating how many pieces to split the data into before bootstrap analysis
#' @param B an integer indicating the number of bootstrapping procedures to perform
#' @param parallelization a logical value indicating whether or not to use parallelization
#' @param spec an integer to specify the number of cores to use while parallelizing
#' @param model_choice the model to be used, only lm1 and glm1 are supported
#'
#' @return an object of class blblm, a bootstrapped linear model fitted to the provided data
#' @example
#' blblm(Sepal.Length~Sepal.Width, iris, m=20, B = 10000, parallelization = TRUE, spec = 3, model_choice = "lm1")
#' @export
blblm <- function(formula,
                  data,
                  m = 10, B = 5000,
                  parallelization = FALSE,
                  spec = 1,
                  model_choice = "lm1") {
  if (is.data.frame(data)) {
    data_list <- split_data(data, m)
    estimates <- map(data_list,
                    ~ lm_each_subsample(
                      formula = formula,
                      data = .,
                      n = nrow(data),
                      B = B,
                      model_choice = model_choice
                    ))
  } else {
    # total number of observations over multiple files
    n = sum(sapply(data, function(x) length(count.fields(x))))
    estimates <- map(
      data,
      ~ lm_each_subsample(
        formula = formula,
        data = .,
        n = n,
        B = B,
        parallelization = parallelization,
        spec = spec,
        model_choice = model_choice
      )
    )
  }
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblm"
  invisible(res)
}


#' split data into m parts of approximated equal sizes
split_data <- function(data, m) {
  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}

#' compute the estimates
lm_each_subsample <-
  function(formula, data, n, B, parallelization = FALSE, spec = 1, model_choice = "lm1") {
    if (parallelization) {
      cl <- makeCluster(spec)
      clusterExport(cl, c("lm_each_boot", "lm1", "blbcoef", "blbsigma"))
      out <- parSapply(cl, seq_len(B), function(i) {
        lm_each_boot(formula, data, n, model_choice)
      })
      stopCluster(cl)
    } else {
      out = replicate(B, lm_each_boot(formula, data, n, model_choice), simplify = FALSE)
    }
    return(out)
  }

#' compute the regression estimates for a blb dataset
lm_each_boot <- function(formula, data, n, model_choice) {
  if (!is.data.frame(data)){
    data = read.csv(data)
  }
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  if (model_choice == "lm1"){
    lm1(formula, data, freqs)
  } else {
    glm1(formula, data, freqs)
  }
}

#' estimate the regression estimates based on the given number of repetitions
lm1 <- function(formula, data, freqs) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wrong variable from the global scope.
  environment(formula) <- environment()
  fit <- lm(formula, data, weights = freqs)
  #fit <- lm_cpp(formula, data, weights = freqs)
  list(coef = blbcoef(fit), sigma = blbsigma(fit))
}

#' estimate the LOGISTIC regression estimates based on the given number of repetitions
glm1 <- function(formula, data, freqs) {
  environment(formula) <- environment()
  fit <- glm(formula, data, weights = freqs, family = "binomial")
  #fit <- lm_cpp(formula, data, weights = freqs)
  list(coef = blbcoef(fit), sigma = blbsigma(fit))
}


#' compute the coefficients from fit

blbcoef <- function(fit) {
  coef(fit)
}


#' compute sigma from fit

blbsigma <- function(fit) {
  p <- fit$rank
  y <- model.extract(fit$model, "response")
  e <- fitted(fit) - y
  w <- fit$weights
  sqrt(sum(w * (e^2)) / (sum(w) - p))
}


#' @export
#' @method print blblm
#' @param x an object of class blblm
#' @return prints the blblm object's formula to the console
print.blblm <- function(x, ...) {
  cat("blblm model:", capture.output(x$formula))
  cat("\n")
}


#' @export
#' @method sigma blblm
#' @param object an object of class blblm
#' @param confidence logical value (default false) wherein the user chooses whether to return the confidence interval or just the sigma
#' @param level a double value between 0 and 1 signifying the alpha level desired
#' @return the sigma of a blblm object and, if confidence = TRUE, its confidence interval
sigma.blblm <- function(object, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  sigma <- mean(map_dbl(est, ~ mean(map_dbl(., "sigma"))))
  if (confidence) {
    alpha <- 1 - 0.95
    limits <- est %>%
      map_mean(~ quantile(map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
      set_names(NULL)
    return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
  } else {
    return(sigma)
  }
}

#' @export
#' @method coef blblm
#' @param object an object of class blblm
#' @return the (Intercept) and coefficient of the independent variables inside the blblm model created previously
coef.blblm <- function(object, ...) {
  est <- object$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
}


#' @export
#' @method confint blblm
#' @param object an object of class blblm
#' @param parm
confint.blblm <- function(object, parm = NULL, level = 0.95, ...) {
  if (is.null(parm)) {
    parm <- attr(terms(object$formula), "term.labels")
  }
  alpha <- 1 - level
  est <- object$estimates
  out <- map_rbind(parm, function(p) {
    map_mean(est, ~ map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)))
  })
  if (is.vector(out)) {
    out <- as.matrix(t(out))
  }
  dimnames(out)[[1]] <- parm
  out
}

#' @export
#' @method predict blblm
#' @param object an object of class blblm
#' @param new_data the new data to apply the blblm model to
#' @param confidence a boolean indicating whether the confidence interval is desired
#' @param level a double between 0 and 1 indicating alpha level
#' @return fitted values for the new_data according to the blblm estimates
predict.blblm <- function(object, new_data, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  X <- model.matrix(reformulate(attr(terms(object$formula), "term.labels")), new_data)
  if (confidence) {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>%
      apply(1, mean_lwr_upr, level = level) %>%
      t())
  } else {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>% rowMeans())
  }
}


mean_lwr_upr <- function(x, level = 0.95) {
  alpha <- 1 - level
  c(fit = mean(x), quantile(x, c(alpha / 2, 1 - alpha / 2)) %>% set_names(c("lwr", "upr")))
}

map_mean <- function(.x, .f, ...) {
  (map(.x, .f, ...) %>% reduce(`+`)) / length(.x)
}

map_cbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(cbind)
}

map_rbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(rbind)
}
