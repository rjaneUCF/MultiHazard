#' @noRd 
revTransform <-
  function (x, data, qu, th = 0, sigma = 1, xi = 0, method = "mixture") {
    if (!is.element(method, c("mixture", "empirical")))
      stop("method should be 'mixture' or 'empirical'")
    
    n <- length(data)
    probs <- (1:n)/(n + 1)
    px <- sapply(x, function(x, p) p[abs(x - p) == min(abs(x - p))][1], p = probs) # take 1st item in case of ties
    px <- as.integer(round(px * (1 + n)))
    res <- sort(data)[px]
    if (method == "mixture") {
      res[res > th] <- u2gpd(x[res > th], p=1-qu, th = th, sigma = sigma, xi = xi)
    }
    res[order(x)] <- sort(res)
    res
  }
