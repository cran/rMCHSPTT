# ============================================================
# CHAIN SAMPLING INSPECTION PLAN - ChSP-1
# ============================================================


#' Chain Sampling Inspection Plan (ChSP-1)
#'
#' Calculates the minimum required sample size for a
#' Chain Sampling Inspection Plan (ChSP-1) under a
#' time-truncated life test.
#'
#' The probability of failure before the termination time is
#' supplied directly by the user through `p`. Therefore,
#' the function is distribution-free and can be used with
#' different lifetime distributions.
#'
#' For a sample of size `n`, let `P0` denote the probability
#' of observing zero defective units and `P1` denote the
#' probability of observing exactly one defective unit:
#'
#' \deqn{
#' P_0=(1-p)^n
#' }{
#' P0=(1-p)^n
#' }
#'
#' and
#'
#' \deqn{
#' P_1=np(1-p)^{n-1}.
#' }{
#' P1=np(1-p)^(n-1)
#' }
#'
#' According to the attribute-case ChSP-1 plan, the probability
#' of acceptance is
#'
#' \deqn{
#' P_a=P_0+P_1P_0^i.
#' }{
#' Pa=P0+P1*P0^i
#' }
#'
#' where `i` is the number of preceding samples considered
#' in the chain.
#'
#' The function determines the minimum sample size `n`
#' satisfying the consumer's risk condition
#'
#' \deqn{
#' P_a\leq\beta.
#' }{
#' Pa<=beta
#' }
#'
#' The lifetime distribution is not embedded in the function.
#' The user calculates the failure probability `p` from the
#' required lifetime distribution and supplies it to the
#' function.
#'
#' @param p User-defined probability of failure before the
#'   termination time. It must lie strictly between 0 and 1.
#' @param a Termination ratio, defined as a = t/theta0.
#'   It must be a positive numeric value.
#' @param b Quality ratio, defined as b = theta/theta0.
#'   It must be a positive numeric value.
#' @param i Number of preceding samples considered in the
#'   chain. It must be a positive integer.
#' @param beta Consumer's risk. It must be a value between
#'   0 and 1.
#' @param n_max Maximum sample size to be searched.
#'   It must be a positive integer.
#'
#' @return A data frame containing `a`, `b`, `p`, `i`, `beta`,
#'   the probabilities `P0` and `P1`, the required sample size
#'   `n`, and the probability of acceptance `Pa`.
#'
#' @examples
#'
#' # ----------------------------------------------------------
#' # Example 1: User-defined failure probabilities
#' # ----------------------------------------------------------
#'
#' p <- c(0.05, 0.10, 0.15, 0.20)
#'
#' chsp1_asip(
#'   p = p,
#'   a = c(0.5, 1, 1.5, 2),
#'   b = 1,
#'   i = 3,
#'   beta = 0.25
#' )
#'
#'
#' # ----------------------------------------------------------
#' # Example 2: Weibull distribution
#' # ----------------------------------------------------------
#'
#' shape <- 2
#' b <- 1
#'
#' a <- c(
#'   0.5, 0.75, 1, 1.25,
#'   1.5, 1.75, 2
#' )
#'
#' p <- 1 - exp(
#'   -((a / b)^shape)
#' )
#'
#' chsp1_asip(
#'   p = p,
#'   a = a,
#'   b = b,
#'   i = 3,
#'   beta = 0.25
#' )
#'
#'
#' # ----------------------------------------------------------
#' # Example 3: Generalized Exponential distribution
#' # Scale parameterization
#' # ----------------------------------------------------------
#'
#' alpha <- 2
#' b <- 1
#'
#' a <- c(
#'   0.5, 0.75, 1, 1.25,
#'   1.5, 1.75, 2
#' )
#'
#' p <- (
#'   1 - exp(-a / b)
#' )^alpha
#'
#' chsp1_asip(
#'   p = p,
#'   a = a,
#'   b = b,
#'   i = 3,
#'   beta = 0.25
#' )
#'
#' @export


chsp1_asip <- function(
    p,
    a,
    b,
    i = 1,
    beta = 0.25,
    n_max = 10000) {

  # ----------------------------------------------------------
  # Input checking
  # ----------------------------------------------------------

  if (!is.numeric(p) || any(!is.finite(p))) {
    stop("'p' must be a numeric value or numeric vector.")
  }

  if (any(p <= 0 | p >= 1)) {
    stop("'p' must lie strictly between 0 and 1.")
  }

  if (!is.numeric(a) ||
      any(!is.finite(a)) ||
      any(a <= 0)) {

    stop("'a' must contain positive finite values.")
  }

  if (!is.numeric(b) ||
      any(!is.finite(b)) ||
      any(b <= 0)) {

    stop("'b' must contain positive finite values.")
  }

  if (!is.numeric(i) ||
      length(i) != 1 ||
      !is.finite(i) ||
      i <= 0 ||
      i != floor(i)) {

    stop("'i' must be a positive integer.")
  }

  if (!is.numeric(beta) ||
      length(beta) != 1 ||
      !is.finite(beta) ||
      beta <= 0 ||
      beta >= 1) {

    stop("'beta' must be between 0 and 1.")
  }

  if (!is.numeric(n_max) ||
      length(n_max) != 1 ||
      !is.finite(n_max) ||
      n_max <= 0 ||
      n_max != floor(n_max)) {

    stop("'n_max' must be a positive integer.")
  }

  i <- as.integer(i)
  n_max <- as.integer(n_max)


  # ----------------------------------------------------------
  # Make a and b compatible with p
  # ----------------------------------------------------------

  if (length(a) == 1) {
    a <- rep(a, length(p))
  }

  if (length(b) == 1) {
    b <- rep(b, length(p))
  }

  if (length(a) != length(p) ||
      length(b) != length(p)) {

    stop("'a', 'b', and 'p' must have compatible lengths.")
  }


  # ----------------------------------------------------------
  # Find minimum sample size
  # ----------------------------------------------------------

  find_n <- function(p_value) {

    for (n in seq_len(n_max)) {

      P0 <- (1 - p_value)^n

      P1 <- n *
        p_value *
        (1 - p_value)^(n - 1)

      Pa <- P0 + P1 * P0^i

      if (Pa <= beta) {

        return(
          c(
            P0 = P0,
            P1 = P1,
            n = n,
            Pa = Pa
          )
        )
      }
    }

    stop(
      "No sample size satisfies Pa <= beta within n_max."
    )
  }


  # ----------------------------------------------------------
  # Calculate n for each p
  # ----------------------------------------------------------

  result <- t(
    sapply(
      p,
      find_n
    )
  )


  # ----------------------------------------------------------
  # Final output
  # ----------------------------------------------------------

  output <- data.frame(
    a = a,
    b = b,
    p = p,
    i = i,
    beta = beta,
    P0 = result[, "P0"],
    P1 = result[, "P1"],
    n = result[, "n"],
    Pa = result[, "Pa"]
  )

  return(output)
}


# ============================================================
# MODIFIED CHAIN SAMPLING INSPECTION PLAN - MChSP-1
# ============================================================


#' Modified Chain Sampling Inspection Plan (MChSP-1)
#'
#' Calculates the minimum required sample size for a
#' Modified Chain Sampling Inspection Plan (MChSP-1)
#' under a time-truncated life test.
#'
#' The probability of failure before the termination time is
#' supplied directly by the user through `p`. Therefore,
#' the function is distribution-free.
#'
#' Let
#'
#' \deqn{
#' P_0=(1-p)^n
#' }{
#' P0=(1-p)^n
#' }
#'
#' denote the probability of observing zero defective units
#' and
#'
#' \deqn{
#' P_1=np(1-p)^{n-1}
#' }{
#' P1=np(1-p)^(n-1)
#' }
#'
#' denote the probability of observing exactly one defective
#' unit.
#'
#' For the attribute-case MChSP-1 plan, the probability of
#' acceptance is
#'
#' \deqn{
#' P_a=P_0^{i+1}+iP_0^iP_1.
#' }{
#' Pa=P0^(i+1)+i*P0^i*P1
#' }
#'
#' where `i` is the number of preceding samples considered
#' in the chain.
#'
#' The function determines the minimum sample size `n`
#' satisfying
#'
#' \deqn{
#' P_a\leq\beta.
#' }{
#' Pa<=beta
#' }
#'
#' @param p User-defined probability of failure before the
#'   termination time. It must lie strictly between 0 and 1.
#' @param a Termination ratio, defined as a = t/theta0.
#'   It must be a positive numeric value.
#' @param b Quality ratio, defined as b = theta/theta0.
#'   It must be a positive numeric value.
#' @param i Number of preceding samples considered in the
#'   chain. It must be a positive integer.
#' @param beta Consumer's risk. It must be a value between
#'   0 and 1.
#' @param n_max Maximum sample size to be searched.
#'   It must be a positive integer.
#'
#' @return A data frame containing `a`, `b`, `p`, `i`, `beta`,
#'   the probabilities `P0` and `P1`, the required sample size
#'   `n`, and the probability of acceptance `Pa`.
#'
#' @examples
#'
#' # ----------------------------------------------------------
#' # Example 1: User-defined failure probabilities
#' # ----------------------------------------------------------
#'
#' p <- c(0.05, 0.10, 0.15, 0.20)
#'
#' mchsp1_asip(
#'   p = p,
#'   a = c(0.5, 1, 1.5, 2),
#'   b = 1,
#'   i = 3,
#'   beta = 0.25
#' )
#'
#'
#' # ----------------------------------------------------------
#' # Example 2: Weibull distribution
#' # ----------------------------------------------------------
#'
#' shape <- 2
#' b <- 1
#'
#' a <- c(
#'   0.5, 0.75, 1, 1.25,
#'   1.5, 1.75, 2
#' )
#'
#' p <- 1 - exp(
#'   -((a / b)^shape)
#' )
#'
#' mchsp1_asip(
#'   p = p,
#'   a = a,
#'   b = b,
#'   i = 3,
#'   beta = 0.25
#' )
#'
#'
#' # ----------------------------------------------------------
#' # Example 3: Generalized Exponential distribution
#' # Scale parameterization
#' # ----------------------------------------------------------
#'
#' alpha <- 2
#' b <- 1
#'
#' a <- c(
#'   0.5, 0.75, 1, 1.25,
#'   1.5, 1.75, 2
#' )
#'
#' p <- (
#'   1 - exp(-a / b)
#' )^alpha
#'
#' mchsp1_asip(
#'   p = p,
#'   a = a,
#'   b = b,
#'   i = 3,
#'   beta = 0.25
#' )
#'
#' @export


mchsp1_asip <- function(
    p,
    a,
    b,
    i = 1,
    beta = 0.25,
    n_max = 10000) {

  # ----------------------------------------------------------
  # Input checking
  # ----------------------------------------------------------

  if (!is.numeric(p) || any(!is.finite(p))) {
    stop("'p' must be a numeric value or numeric vector.")
  }

  if (any(p <= 0 | p >= 1)) {
    stop("'p' must lie strictly between 0 and 1.")
  }

  if (!is.numeric(a) ||
      any(!is.finite(a)) ||
      any(a <= 0)) {

    stop("'a' must contain positive finite values.")
  }

  if (!is.numeric(b) ||
      any(!is.finite(b)) ||
      any(b <= 0)) {

    stop("'b' must contain positive finite values.")
  }

  if (!is.numeric(i) ||
      length(i) != 1 ||
      !is.finite(i) ||
      i <= 0 ||
      i != floor(i)) {

    stop("'i' must be a positive integer.")
  }

  if (!is.numeric(beta) ||
      length(beta) != 1 ||
      !is.finite(beta) ||
      beta <= 0 ||
      beta >= 1) {

    stop("'beta' must be between 0 and 1.")
  }

  if (!is.numeric(n_max) ||
      length(n_max) != 1 ||
      !is.finite(n_max) ||
      n_max <= 0 ||
      n_max != floor(n_max)) {

    stop("'n_max' must be a positive integer.")
  }

  i <- as.integer(i)
  n_max <- as.integer(n_max)


  # ----------------------------------------------------------
  # Make a and b compatible with p
  # ----------------------------------------------------------

  if (length(a) == 1) {
    a <- rep(a, length(p))
  }

  if (length(b) == 1) {
    b <- rep(b, length(p))
  }

  if (length(a) != length(p) ||
      length(b) != length(p)) {

    stop("'a', 'b', and 'p' must have compatible lengths.")
  }


  # ----------------------------------------------------------
  # Find minimum sample size
  # ----------------------------------------------------------

  find_n <- function(p_value) {

    for (n in seq_len(n_max)) {

      P0 <- (1 - p_value)^n

      P1 <- n *
        p_value *
        (1 - p_value)^(n - 1)

      Pa <- P0^(i + 1) +
        i * P0^i * P1

      if (Pa <= beta) {

        return(
          c(
            P0 = P0,
            P1 = P1,
            n = n,
            Pa = Pa
          )
        )
      }
    }

    stop(
      "No sample size satisfies Pa <= beta within n_max."
    )
  }


  # ----------------------------------------------------------
  # Calculate n for each p
  # ----------------------------------------------------------

  result <- t(
    sapply(
      p,
      find_n
    )
  )


  # ----------------------------------------------------------
  # Final output
  # ----------------------------------------------------------

  output <- data.frame(
    a = a,
    b = b,
    p = p,
    i = i,
    beta = beta,
    P0 = result[, "P0"],
    P1 = result[, "P1"],
    n = result[, "n"],
    Pa = result[, "Pa"]
  )

  return(output)
}


# ============================================================
# MODIFIED CHAIN SAMPLING INSPECTION PLAN - MChSP
# ============================================================


#' Modified Chain Sampling Inspection Plan (MChSP)
#'
#' Calculates the minimum required sample size for a
#' Modified Chain Sampling Inspection Plan (MChSP) for
#' inspection by attributes under a time-truncated life test.
#'
#' The probability of failure before the termination time is
#' supplied directly by the user through `p`. Therefore,
#' the function is distribution-free.
#'
#' The MChSP plan for attributes is specified by the three
#' integers `(n, c, i)`, where `n` is the sample size, `c`
#' is the acceptance criterion and `i` is the number of
#' preceding lots considered.
#'
#' Let `D` denote the number of defective units in a sample.
#' Under a binomial model,
#'
#' \deqn{
#' D\sim Binomial(n,p).
#' }{
#' D ~ Binomial(n,p).
#' }
#'
#' Define
#'
#' \deqn{
#' m=P(D\leq c).
#' }{
#' m=P(D<=c).
#' }
#'
#' The acceptance probability of the MChSP plan is
#'
#' \deqn{
#' P_a=m\left[m^i+i m^{i-1}(1-m)\right].
#' }{
#' Pa=m[m^i+i*m^(i-1)*(1-m)].
#' }
#'
#' The function determines the minimum sample size `n`
#' satisfying
#'
#' \deqn{
#' P_a\leq\beta.
#' }{
#' Pa<=beta.
#' }
#'
#' For `c = 0`, the MChSP plan is related to the MChSP-1
#' plan.
#'
#' @param p User-defined probability of failure before the
#'   termination time. It must lie strictly between 0 and 1.
#' @param a Termination ratio, defined as a = t/theta0.
#'   It must be a positive numeric value.
#' @param b Quality ratio, defined as b = theta/theta0.
#'   It must be a positive numeric value.
#' @param i Number of preceding lots considered in the
#'   chain. It must be a positive integer.
#' @param beta Consumer's risk. It must be a value between
#'   0 and 1.
#' @param c Acceptance criterion. It must be a non-negative
#'   integer.
#' @param n_max Maximum sample size to be searched.
#'   It must be a positive integer.
#'
#' @return A data frame containing `a`, `b`, `p`, `i`, `beta`,
#'   `c`, `m`, the required sample size `n`, and the
#'   probability of acceptance `Pa`.
#'
#' @examples
#'
#' # ----------------------------------------------------------
#' # Example 1: User-defined failure probabilities
#' # ----------------------------------------------------------
#'
#' p <- c(0.05, 0.10, 0.15, 0.20)
#'
#' mchsp_asip(
#'   p = p,
#'   a = c(0.5, 1, 1.5, 2),
#'   b = 1,
#'   i = 3,
#'   beta = 0.25,
#'   c = 0
#' )
#'
#'
#' # ----------------------------------------------------------
#' # Example 2: MChSP with c = 1
#' # ----------------------------------------------------------
#'
#' mchsp_asip(
#'   p = p,
#'   a = c(0.5, 1, 1.5, 2),
#'   b = 1,
#'   i = 3,
#'   beta = 0.25,
#'   c = 1
#' )
#'
#'
#' # ----------------------------------------------------------
#' # Example 3: Weibull distribution
#' # ----------------------------------------------------------
#'
#' shape <- 2
#' b <- 1
#'
#' a <- c(
#'   0.5, 0.75, 1, 1.25,
#'   1.5, 1.75, 2
#' )
#'
#' p <- 1 - exp(
#'   -((a / b)^shape)
#' )
#'
#' mchsp_asip(
#'   p = p,
#'   a = a,
#'   b = b,
#'   i = 3,
#'   beta = 0.25,
#'   c = 1
#' )
#'
#'
#' # ----------------------------------------------------------
#' # Example 4: Generalized Exponential distribution
#' # Scale parameterization
#' # ----------------------------------------------------------
#'
#' alpha <- 2
#' b <- 1
#'
#' a <- c(
#'   0.5, 0.75, 1, 1.25,
#'   1.5, 1.75, 2
#' )
#'
#' p <- (
#'   1 - exp(-a / b)
#' )^alpha
#'
#' mchsp_asip(
#'   p = p,
#'   a = a,
#'   b = b,
#'   i = 3,
#'   beta = 0.25,
#'   c = 1
#' )
#'
#' @importFrom stats pbinom
#' @export


mchsp_asip <- function(
    p,
    a,
    b,
    i = 1,
    beta = 0.25,
    c = 0,
    n_max = 10000) {

  # ----------------------------------------------------------
  # Input checking
  # ----------------------------------------------------------

  if (!is.numeric(p) || any(!is.finite(p))) {
    stop("'p' must be a numeric value or numeric vector.")
  }

  if (any(p <= 0 | p >= 1)) {
    stop("'p' must lie strictly between 0 and 1.")
  }

  if (!is.numeric(a) ||
      any(!is.finite(a)) ||
      any(a <= 0)) {

    stop("'a' must contain positive finite values.")
  }

  if (!is.numeric(b) ||
      any(!is.finite(b)) ||
      any(b <= 0)) {

    stop("'b' must contain positive finite values.")
  }

  if (!is.numeric(i) ||
      length(i) != 1 ||
      !is.finite(i) ||
      i <= 0 ||
      i != floor(i)) {

    stop("'i' must be a positive integer.")
  }

  if (!is.numeric(beta) ||
      length(beta) != 1 ||
      !is.finite(beta) ||
      beta <= 0 ||
      beta >= 1) {

    stop("'beta' must be between 0 and 1.")
  }

  if (!is.numeric(c) ||
      length(c) != 1 ||
      !is.finite(c) ||
      c < 0 ||
      c != floor(c)) {

    stop("'c' must be a non-negative integer.")
  }

  if (!is.numeric(n_max) ||
      length(n_max) != 1 ||
      !is.finite(n_max) ||
      n_max <= 0 ||
      n_max != floor(n_max)) {

    stop("'n_max' must be a positive integer.")
  }

  i <- as.integer(i)
  c <- as.integer(c)
  n_max <- as.integer(n_max)


  # ----------------------------------------------------------
  # Make a and b compatible with p
  # ----------------------------------------------------------

  if (length(a) == 1) {
    a <- rep(a, length(p))
  }

  if (length(b) == 1) {
    b <- rep(b, length(p))
  }

  if (length(a) != length(p) ||
      length(b) != length(p)) {

    stop("'a', 'b', and 'p' must have compatible lengths.")
  }


  # ----------------------------------------------------------
  # Find minimum sample size
  # ----------------------------------------------------------

  find_n <- function(p_value) {

    for (n in seq_len(n_max)) {

      if (c > n) {
        next
      }

      # Probability of accepting one lot

      m <- stats::pbinom(
        q = c,
        size = n,
        prob = p_value
      )

      Pa <- m * (
        m^i +
          i * m^(i - 1) *
          (1 - m)
      )

      if (Pa <= beta) {

        return(
          c(
            m = m,
            n = n,
            Pa = Pa
          )
        )
      }
    }

    stop(
      "No sample size satisfies Pa <= beta within n_max."
    )
  }


  # ----------------------------------------------------------
  # Calculate n for each p
  # ----------------------------------------------------------

  result <- t(
    sapply(
      p,
      find_n
    )
  )


  # ----------------------------------------------------------
  # Final output
  # ----------------------------------------------------------

  output <- data.frame(
    a = a,
    b = b,
    p = p,
    i = i,
    beta = beta,
    c = c,
    m = result[, "m"],
    n = result[, "n"],
    Pa = result[, "Pa"]
  )

  return(output)
}


# ============================================================
# MULTIPLE DEPENDENT SAMPLING PLAN - MDS
# ============================================================


#' Multiple Dependent Sampling Plan (MDS) for Attributes
#'
#' Calculates the minimum required sample size for an MDS plan
#' for inspection by attributes under a time-truncated life test.
#'
#' The probability of failure before the termination time is
#' supplied directly by the user through `p`. Therefore, the
#' sampling-plan calculation is distribution-free.
#'
#' The MDS plan is specified by `(n, c1, c2, i)`.
#'
#' Let D denote the number of defective units in a sample.
#' Under the binomial model,
#'
#' \deqn{
#' D\sim Binomial(n,p).
#' }{
#' D ~ Binomial(n,p).
#' }
#'
#' Define
#'
#' \deqn{
#' A=P(D\leq c_1)
#' }{
#' A=P(D<=c1)
#' }
#'
#' and
#'
#' \deqn{
#' B=P(c_1<D\leq c_2).
#' }{
#' B=P(c1<D<=c2)
#' }
#'
#' The new MDS acceptance probability is
#'
#' \deqn{
#' P_a=A+B A^m,
#' }{
#' Pa=A+B*A^m
#' }
#'
#' where `m` is the number of preceding samples. In this
#' implementation, the existing argument `i` is retained for
#' compatibility and represents `m`.
#'
#' Equivalently,
#'
#' \deqn{
#' P_a =
#' \sum_{D=0}^{c_1}{n\choose D}p^D(1-p)^{n-D}
#' +
#' \left[
#' \sum_{D=c_1+1}^{c_2}{n\choose D}p^D(1-p)^{n-D}
#' \right]
#' \left[
#' \sum_{D=0}^{c_1}{n\choose D}p^D(1-p)^{n-D}
#' \right]^m.
#' }{
#' Pa = A + B*A^m
#' }
#'
#' The function determines the minimum sample size `n`
#' satisfying
#'
#' \deqn{
#' P_a\leq\beta.
#' }{
#' Pa<=beta
#' }
#'
#' @param p User-defined probability of failure before the
#'   termination time. It must lie strictly between 0 and 1.
#' @param a Termination ratio, defined as a = t/theta0.
#'   It must be a positive numeric value.
#' @param b Quality ratio, defined as b = theta/theta0.
#'   It must be a positive numeric value.
#' @param i Number of preceding samples considered in the
#'   MDS plan. In the new formulation, `i` represents the
#'   exponent `m` in `P_a=A+B A^m`.
#' @param beta Consumer's risk. It must be a value between
#'   0 and 1.
#' @param c1 First acceptance number. It must be a
#'   non-negative integer.
#' @param c2 Second acceptance number. It must be an integer
#'   greater than or equal to `c1`.
#' @param n_max Maximum sample size to be searched.
#'   It must be a positive integer.
#'
#' @return A data frame containing `a`, `b`, `p`, `i`, `beta`,
#'   `c1`, `c2`, the probabilities `A` and `B`, the required
#'   sample size `n`, and the probability of acceptance `Pa`.
#'
#' @examples
#'
#' # ----------------------------------------------------------
#' # Example 1: User-defined failure probabilities
#' # ----------------------------------------------------------
#'
#' p <- c(0.05, 0.10, 0.15, 0.20)
#'
#' mds_asip(
#'   p = p,
#'   a = c(0.5, 1, 1.5, 2),
#'   b = 1,
#'   i = 3,
#'   beta = 0.25,
#'   c1 = 0,
#'   c2 = 1
#' )
#'
#'
#' # ----------------------------------------------------------
#' # Example 2: Weibull distribution
#' # ----------------------------------------------------------
#'
#' shape <- 2
#' b <- 1
#'
#' a <- c(
#'   0.5, 0.75, 1, 1.25,
#'   1.5, 1.75, 2
#' )
#'
#' p <- 1 - exp(
#'   -((a / b)^shape)
#' )
#'
#' mds_asip(
#'   p = p,
#'   a = a,
#'   b = b,
#'   i = 3,
#'   beta = 0.25,
#'   c1 = 0,
#'   c2 = 1
#' )
#'
#'
#' # ----------------------------------------------------------
#' # Example 3: Generalized Exponential distribution
#' # Scale parameterization
#' # ----------------------------------------------------------
#'
#' alpha <- 2
#' b <- 1
#'
#' a <- c(
#'   0.5, 0.75, 1, 1.25,
#'   1.5, 1.75, 2
#' )
#'
#' p <- (
#'   1 - exp(-a / b)
#' )^alpha
#'
#' mds_asip(
#'   p = p,
#'   a = a,
#'   b = b,
#'   i = 3,
#'   beta = 0.25,
#'   c1 = 0,
#'   c2 = 1
#' )
#'
#' @importFrom stats pbinom
#' @export


mds_asip <- function(
    p,
    a,
    b,
    i = 1,
    beta = 0.25,
    c1 = 0,
    c2 = 1,
    n_max = 10000) {

  # ----------------------------------------------------------
  # Input checking
  # ----------------------------------------------------------

  if (!is.numeric(p) || any(!is.finite(p))) {
    stop("'p' must be a numeric value or numeric vector.")
  }

  if (any(p <= 0 | p >= 1)) {
    stop("'p' must lie strictly between 0 and 1.")
  }

  if (!is.numeric(a) ||
      any(!is.finite(a)) ||
      any(a <= 0)) {

    stop("'a' must contain positive finite values.")
  }

  if (!is.numeric(b) ||
      any(!is.finite(b)) ||
      any(b <= 0)) {

    stop("'b' must contain positive finite values.")
  }

  if (!is.numeric(i) ||
      length(i) != 1 ||
      !is.finite(i) ||
      i <= 0 ||
      i != floor(i)) {

    stop("'i' must be a positive integer.")
  }

  if (!is.numeric(beta) ||
      length(beta) != 1 ||
      !is.finite(beta) ||
      beta <= 0 ||
      beta >= 1) {

    stop("'beta' must be between 0 and 1.")
  }

  if (!is.numeric(c1) ||
      length(c1) != 1 ||
      !is.finite(c1) ||
      c1 < 0 ||
      c1 != floor(c1)) {

    stop("'c1' must be a non-negative integer.")
  }

  if (!is.numeric(c2) ||
      length(c2) != 1 ||
      !is.finite(c2) ||
      c2 < c1 ||
      c2 != floor(c2)) {

    stop("'c2' must be an integer greater than or equal to 'c1'.")
  }

  if (!is.numeric(n_max) ||
      length(n_max) != 1 ||
      !is.finite(n_max) ||
      n_max <= 0 ||
      n_max != floor(n_max)) {

    stop("'n_max' must be a positive integer.")
  }

  i <- as.integer(i)
  c1 <- as.integer(c1)
  c2 <- as.integer(c2)
  n_max <- as.integer(n_max)


  # ----------------------------------------------------------
  # Make a and b compatible with p
  # ----------------------------------------------------------

  if (length(a) == 1) {
    a <- rep(a, length(p))
  }

  if (length(b) == 1) {
    b <- rep(b, length(p))
  }

  if (length(a) != length(p) ||
      length(b) != length(p)) {

    stop("'a', 'b', and 'p' must have compatible lengths.")
  }


  # ----------------------------------------------------------
  # Find minimum sample size
  # New MDS formulation:
  # Pa = A + B*A^i
  # ----------------------------------------------------------

  find_n <- function(p_value) {

    for (n in seq_len(n_max)) {

      if (c2 > n) {
        next
      }

      # Probability of D <= c1
      A <- stats::pbinom(
        q = c1,
        size = n,
        prob = p_value
      )

      # Probability of c1 < D <= c2
      B <- stats::pbinom(
        q = c2,
        size = n,
        prob = p_value
      ) - A

      # New MDS acceptance probability
      Pa <- A + B * A^i

      if (Pa <= beta) {

        return(
          c(
            A = A,
            B = B,
            n = n,
            Pa = Pa
          )
        )
      }
    }

    stop(
      "No sample size satisfies Pa <= beta within n_max."
    )
  }


  # ----------------------------------------------------------
  # Calculate n for each p
  # ----------------------------------------------------------

  result <- t(
    sapply(
      p,
      find_n
    )
  )


  # ----------------------------------------------------------
  # Final output
  # ----------------------------------------------------------

  output <- data.frame(
    a = a,
    b = b,
    p = p,
    i = i,
    beta = beta,
    c1 = c1,
    c2 = c2,
    A = result[, "A"],
    B = result[, "B"],
    n = result[, "n"],
    Pa = result[, "Pa"]
  )

  return(output)
}


# ============================================================
# COMPARISON OF CHAIN SAMPLING INSPECTION PLANS
# ChSP-1, MChSP-1, MDS AND MChSP
# ============================================================


#' Compare ChSP-1, MChSP-1, MDS and MChSP
#'
#' Compares the minimum required sample sizes of four
#' attribute sampling inspection plans:
#' ChSP-1, MChSP-1, MDS and MChSP under a time-truncated
#' life test.
#'
#' The probability of failure before the termination time is
#' supplied directly by the user through `p`. Therefore,
#' the sampling-plan calculations are distribution-free and
#' can be used with different lifetime distributions.
#'
#' All four plans are designed using the single consumer's
#' risk constraint
#'
#' \deqn{
#' P_a\leq\beta.
#' }{
#' Pa<=beta.
#' }
#'
#' For ChSP-1, the probability of acceptance is
#'
#' \deqn{
#' P_a=P_0+P_1P_0^i,
#' }{
#' Pa=P0+P1*P0^i,
#' }
#'
#' where
#'
#' \deqn{
#' P_0=(1-p)^n
#' }{
#' P0=(1-p)^n
#' }
#'
#' and
#'
#' \deqn{
#' P_1=np(1-p)^{n-1}.
#' }{
#' P1=np(1-p)^(n-1).
#' }
#'
#' For MChSP-1, the probability of acceptance is
#'
#' \deqn{
#' P_a=P_0^{i+1}+iP_0^iP_1.
#' }{
#' Pa=P0^(i+1)+i*P0^i*P1.
#' }
#'
#' For MDS, define
#'
#' \deqn{
#' A=P(D\leq c_1)
#' }{
#' A=P(D<=c1)
#' }
#'
#' and
#'
#' \deqn{
#' B=P(c_1<D\leq c_2).
#' }{
#' B=P(c1<D<=c2).
#' }
#'
#' The probability of acceptance satisfies
#'
#' \deqn{
#' P_a=A+B(P_a)^i.
#' }{
#' Pa=A+B(Pa)^i.
#' }
#'
#' For MChSP, let
#'
#' \deqn{
#' m=P(D\leq c).
#' }{
#' m=P(D<=c).
#' }
#'
#' The probability of acceptance is
#'
#' \deqn{
#' P_a=m[m^i+i m^{i-1}(1-m)].
#' }{
#' Pa=m[m^i+i*m^(i-1)*(1-m)].
#' }
#'
#' The function searches for the minimum sample size satisfying
#' the consumer's risk constraint for each plan.
#'
#' @param p User-defined probability of failure before the
#'   termination time. It must lie strictly between 0 and 1.
#' @param a Termination ratio, defined as a = t/theta0.
#'   It must be a positive numeric value.
#' @param b Quality ratio, defined as b = theta/theta0.
#'   It must be a positive numeric value.
#' @param i_ChSP1 Number of preceding samples or lots considered
#'   in the chain for ChSP-1. It must be a positive integer.
#' @param i_MChSP1 Number of preceding samples or lots considered
#'   in the chain for MChSP-1. It must be a positive integer.
#' @param i_MDS Number of preceding samples or lots considered
#'   in the chain for MDS. It must be a positive integer.
#' @param i_MChSP Number of preceding samples or lots considered
#'   in the chain for MChSP. It must be a positive integer.
#' @param beta Consumer's risk. It must be a value between
#'   0 and 1.
#' @param c1 First acceptance number for the MDS plan.
#'   It must be a non-negative integer.
#' @param c2 Second acceptance number for the MDS plan.
#'   It must be an integer greater than or equal to `c1`.
#' @param c Acceptance number for the MChSP plan. It must be
#'   a non-negative integer.
#' @param n_max Maximum sample size to be searched. It must
#'   be a positive integer.
#'
#' @return A data frame containing the design parameters,
#'   the user-specified chain lengths, minimum sample sizes and
#'   probabilities of acceptance for ChSP-1, MChSP-1, MDS and
#'   MChSP.
#'
#' @examples
#'
#' # ----------------------------------------------------------
#' # Example 1: User-defined failure probabilities
#' # ----------------------------------------------------------
#'
#' p <- c(
#'   0.02, 0.04, 0.06, 0.08,
#'   0.10, 0.12, 0.14, 0.16
#' )
#'
#' compare_chain_plans(
#'   p = p,
#'   a = c(
#'     0.5, 0.75, 1, 1.25,
#'     1.5, 1.75, 2, 2.25
#'   ),
#'   b = 1,
#'   i_ChSP1 = 1,
#'   i_MChSP1 = 2,
#'   i_MDS = 3,
#'   i_MChSP = 4,
#'   beta = 0.25,
#'   c1 = 0,
#'   c2 = 1,
#'   c = 1
#' )
#'
#'
#' # ----------------------------------------------------------
#' # Example 2: Weibull distribution
#' # ----------------------------------------------------------
#'
#' shape <- 2
#' b <- 1
#'
#' a <- c(
#'   0.5, 0.75, 1, 1.25,
#'   1.5, 1.75, 2
#' )
#'
#' p <- 1 - exp(
#'   -((a / b)^shape)
#' )
#'
#' compare_chain_plans(
#'   p = p,
#'   a = a,
#'   b = b,
#'   i_ChSP1 = 1,
#'   i_MChSP1 = 2,
#'   i_MDS = 3,
#'   i_MChSP = 4,
#'   beta = 0.25,
#'   c1 = 0,
#'   c2 = 1,
#'   c = 1
#' )
#'
#'
#' # ----------------------------------------------------------
#' # Example 3: Generalized Exponential distribution
#' # Scale parameterization
#' # ----------------------------------------------------------
#'
#' alpha <- 2
#' b <- 1
#'
#' a <- c(
#'   0.5, 0.75, 1, 1.25,
#'   1.5, 1.75, 2
#' )
#'
#' p <- (
#'   1 - exp(-a / b)
#' )^alpha
#'
#' compare_chain_plans(
#'   p = p,
#'   a = a,
#'   b = b,
#'   i_ChSP1 = 1,
#'   i_MChSP1 = 2,
#'   i_MDS = 3,
#'   i_MChSP = 4,
#'   beta = 0.25,
#'   c1 = 0,
#'   c2 = 1,
#'   c = 1
#' )
#'
#' @importFrom stats pbinom
#' @export


compare_chain_plans <- function(
    p,
    a,
    b,
    i_ChSP1 = 1,
    i_MChSP1 = 1,
    i_MDS = 1,
    i_MChSP = 1,
    beta = 0.25,
    c1 = 0,
    c2 = 1,
    c = 1,
    n_max = 10000) {


  # ----------------------------------------------------------
  # Input checking
  # ----------------------------------------------------------

  if (!is.numeric(p) || any(!is.finite(p))) {
    stop("'p' must be a numeric value or numeric vector.")
  }

  if (any(p <= 0 | p >= 1)) {
    stop("'p' must lie strictly between 0 and 1.")
  }

  if (!is.numeric(a) ||
      any(!is.finite(a)) ||
      any(a <= 0)) {

    stop("'a' must contain positive finite values.")
  }

  if (!is.numeric(b) ||
      any(!is.finite(b)) ||
      any(b <= 0)) {

    stop("'b' must contain positive finite values.")
  }

  validate_i <- function(i_value, name) {

    if (!is.numeric(i_value) ||
        length(i_value) != 1 ||
        !is.finite(i_value) ||
        i_value <= 0 ||
        i_value != floor(i_value)) {

      stop(
        paste0("'", name, "' must be a positive integer.")
      )
    }
  }

  validate_i(i_ChSP1, "i_ChSP1")
  validate_i(i_MChSP1, "i_MChSP1")
  validate_i(i_MDS, "i_MDS")
  validate_i(i_MChSP, "i_MChSP")

  if (!is.numeric(beta) ||
      length(beta) != 1 ||
      !is.finite(beta) ||
      beta <= 0 ||
      beta >= 1) {

    stop("'beta' must be between 0 and 1.")
  }

  if (!is.numeric(c1) ||
      length(c1) != 1 ||
      !is.finite(c1) ||
      c1 < 0 ||
      c1 != floor(c1)) {

    stop("'c1' must be a non-negative integer.")
  }

  if (!is.numeric(c2) ||
      length(c2) != 1 ||
      !is.finite(c2) ||
      c2 < c1 ||
      c2 != floor(c2)) {

    stop(
      "'c2' must be an integer greater than or equal to 'c1'."
    )
  }

  if (!is.numeric(c) ||
      length(c) != 1 ||
      !is.finite(c) ||
      c < 0 ||
      c != floor(c)) {

    stop("'c' must be a non-negative integer.")
  }

  if (!is.numeric(n_max) ||
      length(n_max) != 1 ||
      !is.finite(n_max) ||
      n_max <= 0 ||
      n_max != floor(n_max)) {

    stop("'n_max' must be a positive integer.")
  }


  i_ChSP1 <- as.integer(i_ChSP1)
  i_MChSP1 <- as.integer(i_MChSP1)
  i_MDS <- as.integer(i_MDS)
  i_MChSP <- as.integer(i_MChSP)
  c1 <- as.integer(c1)
  c2 <- as.integer(c2)
  c <- as.integer(c)
  n_max <- as.integer(n_max)


  # ----------------------------------------------------------
  # Make a and b compatible with p
  # ----------------------------------------------------------

  if (length(a) == 1) {
    a <- rep(a, length(p))
  }

  if (length(b) == 1) {
    b <- rep(b, length(p))
  }

  if (length(a) != length(p) ||
      length(b) != length(p)) {

    stop("'a', 'b', and 'p' must have compatible lengths.")
  }


  # ----------------------------------------------------------
  # MDS
  # New formulation:
  # Pa = A + B*A^i_MDS
  # ----------------------------------------------------------

  # ----------------------------------------------------------
  # ChSP-1
  # ----------------------------------------------------------

  find_chsp1 <- function(p_value) {

    i <- i_ChSP1

    for (n in seq_len(n_max)) {

      P0 <- (1 - p_value)^n

      P1 <- n *
        p_value *
        (1 - p_value)^(n - 1)

      Pa <- P0 +
        P1 * P0^i

      if (Pa <= beta) {

        return(
          c(
            n = n,
            Pa = Pa
          )
        )
      }
    }

    stop(
      "ChSP-1: no sample size satisfies Pa <= beta."
    )
  }


  # ----------------------------------------------------------
  # MChSP-1
  # ----------------------------------------------------------

  find_mchsp1 <- function(p_value) {

    i <- i_MChSP1

    for (n in seq_len(n_max)) {

      P0 <- (1 - p_value)^n

      P1 <- n *
        p_value *
        (1 - p_value)^(n - 1)

      Pa <- P0^(i + 1) +
        i * P0^i * P1

      if (Pa <= beta) {

        return(
          c(
            n = n,
            Pa = Pa
          )
        )
      }
    }

    stop(
      "MChSP-1: no sample size satisfies Pa <= beta."
    )
  }


  # ----------------------------------------------------------
  # MDS
  # ----------------------------------------------------------

  find_mds <- function(p_value) {

    i <- i_MDS

    for (n in seq_len(n_max)) {

      if (c2 > n) {
        next
      }

      # Probability of D <= c1
      A <- stats::pbinom(
        q = c1,
        size = n,
        prob = p_value
      )

      # Probability of c1 < D <= c2
      B <- stats::pbinom(
        q = c2,
        size = n,
        prob = p_value
      ) - A

      # New MDS acceptance probability:
      # Pa = A + B*A^i
      Pa <- A + B * A^i

      if (Pa <= beta) {

        return(
          c(
            n = n,
            Pa = Pa
          )
        )
      }
    }

    stop(
      "MDS: no sample size satisfies Pa <= beta."
    )
  }


  # ----------------------------------------------------------
  # MChSP
  # ----------------------------------------------------------

  find_mchsp <- function(p_value) {

    i <- i_MChSP

    for (n in seq_len(n_max)) {

      if (c > n) {
        next
      }

      m <- stats::pbinom(
        q = c,
        size = n,
        prob = p_value
      )

      Pa <- m * (
        m^i +
          i * m^(i - 1) *
          (1 - m)
      )

      if (Pa <= beta) {

        return(
          c(
            n = n,
            Pa = Pa
          )
        )
      }
    }

    stop(
      "MChSP: no sample size satisfies Pa <= beta."
    )
  }


  # ----------------------------------------------------------
  # Calculate all four plans
  # ----------------------------------------------------------

  chsp_result <- t(
    sapply(
      p,
      find_chsp1
    )
  )

  mchsp1_result <- t(
    sapply(
      p,
      find_mchsp1
    )
  )

  mds_result <- t(
    sapply(
      p,
      find_mds
    )
  )

  mchsp_result <- t(
    sapply(
      p,
      find_mchsp
    )
  )


  # ----------------------------------------------------------
  # Final output
  # ----------------------------------------------------------

  output <- data.frame(

    a = a,
    b = b,
    p = p,

    i_ChSP1 = i_ChSP1,
    i_MChSP1 = i_MChSP1,
    i_MDS = i_MDS,
    i_MChSP = i_MChSP,
    beta = beta,

    c1 = c1,
    c2 = c2,
    c = c,

    n_ChSP1 =
      chsp_result[, "n"],

    Pa_ChSP1 =
      chsp_result[, "Pa"],

    n_MChSP1 =
      mchsp1_result[, "n"],

    Pa_MChSP1 =
      mchsp1_result[, "Pa"],

    n_MDS =
      mds_result[, "n"],

    Pa_MDS =
      mds_result[, "Pa"],

    n_MChSP =
      mchsp_result[, "n"],

    Pa_MChSP =
      mchsp_result[, "Pa"]
  )

  return(output)
}


# ============================================================
# PLOT COMPARISON OF FOUR CHAIN SAMPLING PLANS
# ============================================================


#' Plot Comparison of ChSP-1, MChSP-1, MDS and MChSP
#'
#' Produces a comparison plot of the minimum sample sizes
#' required by ChSP-1, MChSP-1, MDS and MChSP.
#'
#' The probability of failure before the termination time is
#' supplied by the user through `p`.
#'
#' The x-axis represents the termination ratio
#' \eqn{a=t/theta0} and the y-axis represents the minimum
#' required sample size.
#'
#' All plans satisfy the common consumer's risk constraint
#'
#' \deqn{
#' P_a\leq\beta.
#' }{
#' Pa<=beta.
#' }
#'
#' @param p User-defined probability of failure before the
#'   termination time. It must lie strictly between 0 and 1.
#' @param a Termination ratio, defined as a = t/theta0.
#'   It must be a positive numeric value.
#' @param b Quality ratio, defined as b = theta/theta0.
#'   It must be a positive numeric value.
#' @param i_ChSP1 Number of preceding samples or lots considered
#'   in the chain for ChSP-1. It must be a positive integer.
#' @param i_MChSP1 Number of preceding samples or lots considered
#'   in the chain for MChSP-1. It must be a positive integer.
#' @param i_MDS Number of preceding samples or lots considered
#'   in the chain for MDS. It must be a positive integer.
#' @param i_MChSP Number of preceding samples or lots considered
#'   in the chain for MChSP. It must be a positive integer.
#' @param beta Consumer's risk. It must be a value between
#'   0 and 1.
#' @param c1 First acceptance number for MDS.
#' @param c2 Second acceptance number for MDS.
#' @param c Acceptance number for MChSP.
#' @param n_max Maximum sample size to be searched.
#' @param ylim Numeric vector of length two specifying the limits
#'   of the y-axis. Default is `NULL`.
#'
#' @return The comparison data frame is returned invisibly.
#'   The function also produces the comparison plot.
#'
#' @examples
#'
#' # ----------------------------------------------------------
#' # Example 1: User-defined failure probabilities
#' # ----------------------------------------------------------
#'
#' p <- c(
#'   0.02, 0.04, 0.06, 0.08,
#'   0.10, 0.12, 0.14, 0.16
#' )
#'
#' a <- c(
#'   0.5, 0.75, 1, 1.25,
#'   1.5, 1.75, 2, 2.25
#' )
#'
#' plot_chain_plans(
#'   p = p,
#'   a = a,
#'   b = 1,
#'   i_ChSP1 = 1,
#'   i_MChSP1 = 2,
#'   i_MDS = 3,
#'   i_MChSP = 4,
#'   beta = 0.25,
#'   c1 = 0,
#'   c2 = 1,
#'   c = 1
#' )
#'
#'
#' # ----------------------------------------------------------
#' # Example 2: Weibull distribution
#' # ----------------------------------------------------------
#'
#' shape <- 2
#' b <- 1
#'
#' a <- c(
#'   0.5, 0.75, 1, 1.25,
#'   1.5, 1.75, 2
#' )
#'
#' p <- 1 - exp(
#'   -((a / b)^shape)
#' )
#'
#' plot_chain_plans(
#'   p = p,
#'   a = a,
#'   b = b,
#'   i_ChSP1 = 1,
#'   i_MChSP1 = 2,
#'   i_MDS = 3,
#'   i_MChSP = 4,
#'   beta = 0.25,
#'   c1 = 0,
#'   c2 = 1,
#'   c = 1
#' )
#'
#'
#' # ----------------------------------------------------------
#' # Example 3: Generalized Exponential distribution
#' # Scale parameterization
#' # ----------------------------------------------------------
#'
#' alpha <- 2
#' b <- 0.5
#'
#' a <- c(
#'   0.5, 0.75, 1, 1.25,
#'   1.5, 1.75, 2
#' )
#'
#' p <- (
#'   1 - exp(-a / b)
#' )^alpha
#'
#' plot_chain_plans(
#'   p = p,
#'   a = a,
#'   b = b,
#'   i_ChSP1 = 1,
#'   i_MChSP1 = 2,
#'   i_MDS = 3,
#'   i_MChSP = 4,
#'   beta = 0.25,
#'   c1 = 0,
#'   c2 = 1,
#'   c = 1
#' )
#'
#' @importFrom graphics plot lines legend
#' @export


plot_chain_plans <- function(
    p,
    a,
    b,
    i_ChSP1 = 1,
    i_MChSP1 = 1,
    i_MDS = 1,
    i_MChSP = 1,
    beta = 0.25,
    c1 = 0,
    c2 = 1,
    c = 1,
    n_max = 10000,
    ylim = NULL) {


  # ----------------------------------------------------------
  # Calculate comparison
  # ----------------------------------------------------------

  result <- compare_chain_plans(
    p = p,
    a = a,
    b = b,
    i_ChSP1 = i_ChSP1,
    i_MChSP1 = i_MChSP1,
    i_MDS = i_MDS,
    i_MChSP = i_MChSP,
    beta = beta,
    c1 = c1,
    c2 = c2,
    c = c,
    n_max = n_max
  )


  # ----------------------------------------------------------
  # Plot ChSP-1
  # ----------------------------------------------------------

  plot(
    result$a,
    result$n_ChSP1,
    type = "o",
    pch = 16,
    lwd = 2,
    xlab = "Termination Ratio (a)",
    ylab = "Minimum Sample Size",
    ylim = ylim,
    main = paste(
      "Comparison of Chain Sampling Plans",
      "\nConsumer's Risk =", beta,
      ", i = (",
      i_ChSP1, ", ",
      i_MChSP1, ", ",
      i_MDS, ", ",
      i_MChSP, ")"
    )
  )


  # ----------------------------------------------------------
  # Add MChSP-1
  # ----------------------------------------------------------

  lines(
    result$a,
    result$n_MChSP1,
    type = "o",
    pch = 17,
    lwd = 2
  )


  # ----------------------------------------------------------
  # Add MDS
  # ----------------------------------------------------------

  lines(
    result$a,
    result$n_MDS,
    type = "o",
    pch = 15,
    lwd = 2
  )


  # ----------------------------------------------------------
  # Add MChSP
  # ----------------------------------------------------------

  lines(
    result$a,
    result$n_MChSP,
    type = "o",
    pch = 18,
    lwd = 2
  )


  # ----------------------------------------------------------
  # Legend
  # ----------------------------------------------------------

  legend(
    "topright",
    legend = c(
      "ChSP-1",
      "MChSP-1",
      "MDS",
      "MChSP"
    ),
    pch = c(
      16,
      17,
      15,
      18
    ),
    lty = 1,
    lwd = 2,
    bty = "n"
  )


  # ----------------------------------------------------------
  # Return result
  # ----------------------------------------------------------

  invisible(result)
}


 #-----------------------------------------------------------

# ============================================================
# MChSP: OC VALUES
# ============================================================

#' Operating Characteristic Values for an MChSP Plan
#'
#' Calculates operating characteristic (OC) values for a
#' Modified Chain Sampling Inspection Plan (MChSP).
#'
#' The sample size is first determined at the design condition
#' using `p_design` and the constraint `Pa <= beta`. The resulting
#' sample size is then kept fixed while `p_oc` is used to calculate
#' the probability of acceptance for the quality ratios in `b_oc`.
#'
#' Let
#'
#' \deqn{m=P(D\leq c)}
#' \deqn{P_a=m[m^i+i m^{i-1}(1-m)].}
#'
#' The resulting `Pa` values are the OC values for the specified
#' quality ratios.
#'
#' @param p_design Failure probability at the design condition.
#'   It must lie strictly between 0 and 1.
#' @param a Termination ratio, defined as `a = t/theta0`.
#' @param p_oc Matrix or data frame of failure probabilities
#'   used for OC calculation. Rows correspond to `a` and columns
#'   correspond to `b_oc`.
#' @param b_oc Quality ratios used for OC calculation.
#' @param i Number of preceding lots considered in the MChSP plan.
#' @param beta Consumer's risk. It must be between 0 and 1.
#' @param c Acceptance criterion. It must be a non-negative integer.
#' @param n_max Maximum sample size searched at the design condition.
#'
#' @return A data frame containing `a`, `b`, `p`, `i`, `beta`, `c`,
#'   fixed `n`, `ASN`, `m`, and `Pa`. The `Pa` values are the OC
#'   values.
#'
#' @examples
#' a <- c(0.5, 1, 1.5, 2)
#'
#' p_design <- c(0.05, 0.10, 0.15, 0.20)
#'
#' b_oc <- 2:12
#'
#' p_oc <- outer(
#'   a, b_oc,
#'   function(a, b) pmin(0.95, a / (10 * b))
#' )
#'
#' mchsp_oc(
#'   p_design = p_design,
#'   a = a,
#'   p_oc = p_oc,
#'   b_oc = b_oc,
#'   i = 3,
#'   beta = 0.25,
#'   c = 0
#' )
#'
#' @importFrom stats pbinom
#' @export

mchsp_oc <- function(
    p_design,
    a,
    p_oc,
    b_oc,
    i = 1,
    beta = 0.25,
    c = 0,
    n_max = 10000) {

  # ----------------------------------------------------------
  # Input checking
  # ----------------------------------------------------------

  if (!is.numeric(p_design) ||
      any(!is.finite(p_design))) {
    stop("'p_design' must be numeric and finite.")
  }

  if (any(p_design <= 0 | p_design >= 1)) {
    stop("'p_design' must lie strictly between 0 and 1.")
  }

  if (!is.numeric(a) ||
      any(!is.finite(a)) ||
      any(a <= 0)) {
    stop("'a' must contain positive finite values.")
  }

  if (!is.numeric(p_oc) ||
      any(!is.finite(p_oc))) {
    stop("'p_oc' must be numeric and finite.")
  }

  if (any(p_oc <= 0 | p_oc >= 1)) {
    stop("All values of 'p_oc' must lie strictly between 0 and 1.")
  }

  if (!is.numeric(b_oc) ||
      any(!is.finite(b_oc)) ||
      any(b_oc <= 0)) {
    stop("'b_oc' must contain positive finite values.")
  }

  if (!is.numeric(i) ||
      length(i) != 1 ||
      !is.finite(i) ||
      i <= 0 ||
      i != floor(i)) {
    stop("'i' must be a positive integer.")
  }

  if (!is.numeric(beta) ||
      length(beta) != 1 ||
      !is.finite(beta) ||
      beta <= 0 ||
      beta >= 1) {
    stop("'beta' must be between 0 and 1.")
  }

  if (!is.numeric(c) ||
      length(c) != 1 ||
      !is.finite(c) ||
      c < 0 ||
      c != floor(c)) {
    stop("'c' must be a non-negative integer.")
  }

  if (!is.numeric(n_max) ||
      length(n_max) != 1 ||
      !is.finite(n_max) ||
      n_max <= 0 ||
      n_max != floor(n_max)) {
    stop("'n_max' must be a positive integer.")
  }

  p_oc <- as.matrix(p_oc)

  if (length(p_design) != length(a)) {
    stop("'p_design' and 'a' must have the same length.")
  }

  if (nrow(p_oc) != length(a)) {
    stop("The number of rows of 'p_oc' must equal the length of 'a'.")
  }

  if (ncol(p_oc) != length(b_oc)) {
    stop("The number of columns of 'p_oc' must equal the length of 'b_oc'.")
  }

  i <- as.integer(i)
  c <- as.integer(c)
  n_max <- as.integer(n_max)

  # ----------------------------------------------------------
  # Find design sample size
  # ----------------------------------------------------------

  find_n <- function(p_value) {

    for (n in seq_len(n_max)) {

      if (c > n) {
        next
      }

      m <- stats::pbinom(
        q = c,
        size = n,
        prob = p_value
      )

      Pa <- m * (
        m^i +
          i * m^(i - 1) * (1 - m)
      )

      if (is.finite(Pa) && Pa <= beta) {
        return(n)
      }
    }

    stop("No sample size satisfies Pa <= beta within n_max.")
  }

  n_design <- sapply(p_design, find_n)

  # ----------------------------------------------------------
  # Calculate OC values using fixed design sample size
  # ----------------------------------------------------------

  output_list <- vector(
    "list",
    length(a) * length(b_oc)
  )

  counter <- 1

  for (j in seq_along(a)) {

    n <- n_design[j]

    for (k in seq_along(b_oc)) {

      p_value <- p_oc[j, k]

      m <- stats::pbinom(
        q = c,
        size = n,
        prob = p_value
      )

      Pa <- m * (
        m^i +
          i * m^(i - 1) * (1 - m)
      )

      output_list[[counter]] <- data.frame(
        a = a[j],
        b = b_oc[k],
        p = p_value,
        i = i,
        beta = beta,
        c = c,
        n = as.integer(n),
        ASN = as.integer(n),
        m = m,
        Pa = Pa
      )

      counter <- counter + 1
    }
  }

  result <- do.call(rbind, output_list)
  rownames(result) <- NULL

  result
}


# ============================================================
# Plot MChSP OC Curve
# ============================================================

#' Plot OC Values for an MChSP Plan
#'
#' Produces a base R OC curve using the `Pa` values returned by
#' `mchsp_oc()`.
#'
#' Different termination ratios `a` are distinguished using
#' different plotting symbols and line types.
#'
#' @param x Output from `mchsp_oc()`.
#' @param ... Additional graphical arguments passed to `plot()`.
#'
#' @return Invisibly returns the supplied data frame.
#'
#' @examples
#' a <- c(0.5, 1, 1.5, 2)
#' p_design <- c(0.05, 0.10, 0.15, 0.20)
#' b_oc <- 2:12
#'
#' p_oc <- sapply(
#'   b_oc,
#'   function(b) a / (10 * b)
#' )
#'
#' result <- mchsp_oc(
#'   p_design = p_design,
#'   a = a,
#'   p_oc = p_oc,
#'   b_oc = b_oc,
#'   i = 3,
#'   beta = 0.25,
#'   c = 0
#' )
#'
#' plot_mchsp_oc(result)
#'
#' @export

plot_mchsp_oc <- function(x, ...) {

  if (!is.data.frame(x)) {
    stop("'x' must be a data frame returned by 'mchsp_oc()'.")
  }

  required <- c("a", "b", "Pa")

  if (!all(required %in% names(x))) {
    stop("'x' must contain columns 'a', 'b', and 'Pa'.")
  }

  a_values <- unique(x$a)
  n_a <- length(a_values)

  pch_values <- c(
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
    11, 12, 13, 14, 15, 16, 17, 18
  )

  lty_values <- c(1, 2, 3, 4, 5, 6)

  pch_values <- rep(
    pch_values,
    length.out = n_a
  )

  lty_values <- rep(
    lty_values,
    length.out = n_a
  )

  first <- a_values[1]

  current <- x[x$a == first, ]

  plot(
    current$b,
    current$Pa,
    type = "b",
    pch = pch_values[1],
    lty = lty_values[1],
    lwd = 1.2,
    ylim = c(0, 1),
    xlab = "Quality ratio (b)",
    ylab = "Probability of acceptance (Pa)",
    main = "OC Curve of MChSP Plan",
    ...
  )

  if (n_a > 1) {

    for (j in 2:n_a) {

      current <- x[x$a == a_values[j], ]

      graphics::lines(
        current$b,
        current$Pa,
        type = "b",
        pch = pch_values[j],
        lty = lty_values[j],
        lwd = 1.2
      )
    }
  }

  graphics::legend(
    "bottomright",
    legend = paste0("a = ", a_values),
    pch = pch_values,
    lty = lty_values,
    lwd = 1.2,
    bty = "n"
  )

  invisible(x)
}

