# ============================================================
# TESTS FOR CHAIN SAMPLING INSPECTION PLANS
# ChSP-1, MChSP-1, MDS (new formulation), and MChSP
# ============================================================

library(testthat)


# ------------------------------------------------------------
# 1. ChSP-1
# ------------------------------------------------------------

test_that("ChSP-1 satisfies Pa <= beta", {

  p <- c(0.05, 0.10, 0.15)

  result <- chsp1_asip(
    p = p,
    a = c(0.5, 1, 1.5),
    b = 1,
    i = 3,
    beta = 0.25,
    n_max = 1000
  )

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 3)
  expect_true(all(result$Pa <= 0.25))
  expect_true(all(result$n >= 1))
})


# ------------------------------------------------------------
# 2. MChSP-1
# ------------------------------------------------------------

test_that("MChSP-1 satisfies Pa <= beta", {

  p <- c(0.05, 0.10, 0.15)

  result <- mchsp1_asip(
    p = p,
    a = c(0.5, 1, 1.5),
    b = 1,
    i = 3,
    beta = 0.25,
    n_max = 1000
  )

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 3)
  expect_true(all(result$Pa <= 0.25))
  expect_true(all(result$n >= 1))
})


# ------------------------------------------------------------
# 3. MDS -- NEW FORMULATION
#
#     Pa = A + B*A^m
#
# In the package, i is retained as the argument representing m.
# ------------------------------------------------------------

test_that("MDS uses the new Pa = A + B*A^i formulation", {

  p <- 0.10
  beta <- 0.25
  i <- 3
  c1 <- 0
  c2 <- 1

  result <- mds_asip(
    p = p,
    a = 1,
    b = 1,
    i = i,
    beta = beta,
    c1 = c1,
    c2 = c2,
    n_max = 1000
  )

  n <- result$n

  A <- pbinom(
    q = c1,
    size = n,
    prob = p
  )

  B <- pbinom(
    q = c2,
    size = n,
    prob = p
  ) - A

  # Equation from the new MDS method:
  # Pa = A + B*A^i
  expected_Pa <- A + B * A^i

  expect_equal(
    result$A,
    A,
    tolerance = 1e-12
  )

  expect_equal(
    result$B,
    B,
    tolerance = 1e-12
  )

  expect_equal(
    result$Pa,
    expected_Pa,
    tolerance = 1e-12
  )

  expect_lte(
    result$Pa,
    beta + 1e-12
  )
})


# ------------------------------------------------------------
# 4. MDS -- MINIMUM SAMPLE SIZE
# ------------------------------------------------------------

test_that("MDS returns the minimum n satisfying Pa <= beta", {

  p <- 0.10
  beta <- 0.25
  i <- 3
  c1 <- 0
  c2 <- 1

  result <- mds_asip(
    p = p,
    a = 1,
    b = 1,
    i = i,
    beta = beta,
    c1 = c1,
    c2 = c2,
    n_max = 1000
  )

  n <- result$n

  expect_lte(result$Pa, beta)

  if (n > c2) {

    A_previous <- pbinom(
      q = c1,
      size = n - 1,
      prob = p
    )

    B_previous <- pbinom(
      q = c2,
      size = n - 1,
      prob = p
    ) - A_previous

    Pa_previous <- A_previous +
      B_previous * A_previous^i

    expect_gt(
      Pa_previous,
      beta
    )
  }
})


# ------------------------------------------------------------
# 5. MDS -- VECTOR INPUT
# ------------------------------------------------------------

test_that("MDS accepts vector p, a and scalar b", {

  p <- c(0.05, 0.10, 0.15, 0.20)
  a <- c(0.5, 1, 1.5, 2)

  result <- mds_asip(
    p = p,
    a = a,
    b = 1,
    i = 3,
    beta = 0.25,
    c1 = 0,
    c2 = 1,
    n_max = 1000
  )

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), length(p))
  expect_equal(result$p, p)
  expect_equal(result$a, a)
  expect_true(all(result$Pa <= 0.25))
  expect_true(all(result$n >= 2))
})


# ------------------------------------------------------------
# 6. MDS -- DIFFERENT c1 AND c2
# ------------------------------------------------------------

test_that("MDS works for general c1 and c2", {

  result <- mds_asip(
    p = 0.10,
    a = 1,
    b = 1,
    i = 2,
    beta = 0.25,
    c1 = 1,
    c2 = 2,
    n_max = 1000
  )

  n <- result$n

  A <- pbinom(1, n, 0.10)
  B <- pbinom(2, n, 0.10) - A

  expect_equal(
    result$Pa,
    A + B * A^2,
    tolerance = 1e-12
  )

  expect_lte(result$Pa, 0.25 + 1e-12)
})


# ------------------------------------------------------------
# 7. MDS -- INPUT VALIDATION
# ------------------------------------------------------------

test_that("MDS checks invalid inputs", {

  expect_error(
    mds_asip(
      p = 0.10,
      a = 1,
      b = 1,
      i = 0,
      beta = 0.25,
      c1 = 0,
      c2 = 1
    )
  )

  expect_error(
    mds_asip(
      p = 0.10,
      a = 1,
      b = 1,
      i = 2,
      beta = 0.25,
      c1 = 2,
      c2 = 1
    )
  )

  expect_error(
    mds_asip(
      p = 0.10,
      a = 1,
      b = 1,
      i = 2,
      beta = 0.25,
      c1 = -1,
      c2 = 1
    )
  )

  expect_error(
    mds_asip(
      p = 1.10,
      a = 1,
      b = 1,
      i = 2,
      beta = 0.25,
      c1 = 0,
      c2 = 1
    )
  )
})


# ------------------------------------------------------------
# 8. MChSP
# ------------------------------------------------------------

test_that("MChSP satisfies Pa <= beta", {

  result <- mchsp_asip(
    p = c(0.05, 0.10, 0.15),
    a = c(0.5, 1, 1.5),
    b = 1,
    i = 3,
    beta = 0.25,
    c = 1,
    n_max = 1000
  )

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 3)
  expect_true(all(result$Pa <= 0.25))
})


# ------------------------------------------------------------
# 9. COMPARISON FUNCTION
#
# The MDS part of compare_chain_plans() must use the same
# new formulation Pa = A + B*A^i_MDS.
# ------------------------------------------------------------

test_that("comparison uses the new MDS formulation", {

  p <- c(0.05, 0.10, 0.15)

  comparison <- compare_chain_plans(
    p = p,
    a = c(0.5, 1, 1.5),
    b = 1,
    i_ChSP1 = 1,
    i_MChSP1 = 2,
    i_MDS = 3,
    i_MChSP = 4,
    beta = 0.25,
    c1 = 0,
    c2 = 1,
    c = 1,
    n_max = 1000
  )

  mds <- mds_asip(
    p = p,
    a = c(0.5, 1, 1.5),
    b = 1,
    i = 3,
    beta = 0.25,
    c1 = 0,
    c2 = 1,
    n_max = 1000
  )

  expect_equal(
    comparison$n_MDS,
    mds$n
  )

  expect_equal(
    comparison$Pa_MDS,
    mds$Pa,
    tolerance = 1e-12
  )

  expect_true(
    all(comparison$Pa_MDS <= 0.25)
  )
})


# ------------------------------------------------------------
# 10. COMPARISON OUTPUT
# ------------------------------------------------------------

test_that("comparison returns all four plans", {

  result <- compare_chain_plans(
    p = c(0.05, 0.10),
    a = c(0.5, 1),
    b = 1,
    i_ChSP1 = 1,
    i_MChSP1 = 2,
    i_MDS = 3,
    i_MChSP = 4,
    beta = 0.25,
    c1 = 0,
    c2 = 1,
    c = 1,
    n_max = 1000
  )

  required_columns <- c(
    "n_ChSP1",
    "Pa_ChSP1",
    "n_MChSP1",
    "Pa_MChSP1",
    "n_MDS",
    "Pa_MDS",
    "n_MChSP",
    "Pa_MChSP"
  )

  expect_true(
    all(required_columns %in% names(result))
  )

  expect_true(all(result$Pa_ChSP1 <= 0.25))
  expect_true(all(result$Pa_MChSP1 <= 0.25))
  expect_true(all(result$Pa_MDS <= 0.25))
  expect_true(all(result$Pa_MChSP <= 0.25))
})


# ------------------------------------------------------------
# 11. PLOT FUNCTION
# ------------------------------------------------------------

test_that("comparison plot runs without error", {

  test_file <- tempfile(
    pattern = "chain_plan_plot_",
    fileext = ".pdf"
  )

  pdf(test_file)

  result <- plot_chain_plans(
    p = c(0.05, 0.10, 0.15),
    a = c(0.5, 1, 1.5),
    b = 1,
    i_ChSP1 = 1,
    i_MChSP1 = 2,
    i_MDS = 3,
    i_MChSP = 4,
    beta = 0.25,
    c1 = 0,
    c2 = 1,
    c = 1,
    n_max = 1000
  )

  dev.off()

  expect_true(
    file.exists(test_file)
  )

  expect_s3_class(
    result,
    "data.frame"
  )

  unlink(test_file)
})

#--------------------------
#   NEW

# ============================================================
# Tests for MChSP OC values
# ============================================================

test_that("mchsp_oc returns the expected structure", {

  a <- c(0.5, 1, 1.5, 2)

  p_design <- c(
    0.05, 0.10, 0.15, 0.20
  )

  b_oc <- 2:5

  p_oc <- sapply(
    b_oc,
    function(b) a / (10 * b)
  )

  result <- mchsp_oc(
    p_design = p_design,
    a = a,
    p_oc = p_oc,
    b_oc = b_oc,
    i = 3,
    beta = 0.25,
    c = 0
  )

  expect_s3_class(result, "data.frame")

  expect_equal(
    names(result),
    c(
      "a", "b", "p", "i", "beta",
      "c", "n", "ASN", "m", "Pa"
    )
  )

  expect_equal(
    nrow(result),
    length(a) * length(b_oc)
  )
})


test_that("MChSP OC uses the fixed design sample size", {

  a <- c(0.5, 1, 1.5, 2)

  p_design <- c(
    0.05, 0.10, 0.15, 0.20
  )

  b_oc <- 2:5

  p_oc <- sapply(
    b_oc,
    function(b) a / (10 * b)
  )

  result <- mchsp_oc(
    p_design = p_design,
    a = a,
    p_oc = p_oc,
    b_oc = b_oc,
    i = 3,
    beta = 0.25,
    c = 0
  )

  for (j in seq_along(a)) {

    n_values <- result$n[result$a == a[j]]

    expect_length(
      unique(n_values),
      1
    )

    expect_equal(
      unique(result$ASN[result$a == a[j]]),
      unique(n_values)
    )
  }
})


test_that("MChSP OC probability follows the stated formula", {

  a <- c(0.5, 1)

  p_design <- c(0.05, 0.10)

  b_oc <- c(2, 4)

  p_oc <- matrix(
    c(
      0.02, 0.01,
      0.04, 0.02
    ),
    nrow = 2,
    byrow = TRUE
  )

  result <- mchsp_oc(
    p_design = p_design,
    a = a,
    p_oc = p_oc,
    b_oc = b_oc,
    i = 3,
    beta = 0.25,
    c = 0
  )

  for (r in seq_len(nrow(result))) {

    m <- result$m[r]

    expected_pa <- m * (
      m^3 +
        3 * m^2 * (1 - m)
    )

    expect_equal(
      result$Pa[r],
      expected_pa,
      tolerance = 1e-12
    )
  }
})


test_that("MChSP OC validates p_design and p_oc", {

  a <- c(0.5, 1)

  b_oc <- c(2, 4)

  p_oc <- matrix(
    c(
      0.02, 0.01,
      0.04, 0.02
    ),
    nrow = 2,
    byrow = TRUE
  )

  expect_error(
    mchsp_oc(
      p_design = c(0, 0.1),
      a = a,
      p_oc = p_oc,
      b_oc = b_oc
    ),
    "p_design"
  )

  expect_error(
    mchsp_oc(
      p_design = c(0.05, 0.1),
      a = a,
      p_oc = matrix(
        c(0.02, 1.1, 0.04, 0.02),
        nrow = 2
      ),
      b_oc = b_oc
    ),
    "p_oc"
  )
})


test_that("MChSP OC validates matrix dimensions", {

  expect_error(
    mchsp_oc(
      p_design = c(0.05, 0.10),
      a = c(0.5, 1),
      p_oc = matrix(
        c(0.02, 0.01, 0.04),
        nrow = 1
      ),
      b_oc = c(2, 4)
    ),
    "rows"
  )

  expect_error(
    mchsp_oc(
      p_design = c(0.05, 0.10),
      a = c(0.5, 1),
      p_oc = matrix(
        c(0.02, 0.01, 0.04, 0.02),
        nrow = 2
      ),
      b_oc = c(2, 4, 6)
    ),
    "columns"
  )
})


test_that("MChSP OC rejects invalid plan parameters", {

  a <- c(0.5, 1)

  p_design <- c(0.05, 0.10)

  b_oc <- c(2, 4)

  p_oc <- matrix(
    c(0.02, 0.01, 0.04, 0.02),
    nrow = 2,
    byrow = TRUE
  )

  expect_error(
    mchsp_oc(
      p_design = p_design,
      a = a,
      p_oc = p_oc,
      b_oc = b_oc,
      i = 0
    ),
    "i"
  )

  expect_error(
    mchsp_oc(
      p_design = p_design,
      a = a,
      p_oc = p_oc,
      b_oc = b_oc,
      beta = 1
    ),
    "beta"
  )

  expect_error(
    mchsp_oc(
      p_design = p_design,
      a = a,
      p_oc = p_oc,
      b_oc = b_oc,
      c = -1
    ),
    "c"
  )
})


# ============================================================
# END OF TEST FILE
# ============================================================

