# Byte-level snapshots of the printed record.
#
# `test-print-methods.R` asserts the printing contract structurally --- a
# `<class>` header, a line budget, no closure, environment, or whole digest ---
# and says why: a structural assertion holds across BLAS implementations where
# a byte-for-byte one may not. This file supplies the other half of the
# contract for the classes a first-hour user actually holds. `design/api-tiers.md`
# lists twenty-one core-tier exports; every value class they return is pinned
# here, so a change to how a record renders has to be argued for in a diff
# rather than discovered by a reader.
#
# The rule that keeps these stable is that nothing snapshotted here contains a
# number produced by linear algebra. Declaration records --- domains, frames,
# spaces, relations, pairings, plans --- print counts, identifiers, and content
# signatures, all of which are fixed by the inputs. The four view classes print
# a numeric preview of a computed result, so for those it is the deterministic
# one-line `format()` that is pinned, and the preview itself stays under the
# structural assertions next door.

snapshot_fixture <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) {
      return(cached)
    }
    example <- example_fmri_effects()
    relation <- example$fit$relation
    pairing <- cross_partitions(relation, independence = "independent")
    plan <- plan_geometry(relation, example$frame, pairing)
    cached <<- list(
      example = example,
      relation = relation,
      pairing = pairing,
      plan = plan,
      geometry = materialize_geometry(plan)
    )
    cached
  }
})

# Declarations ---------------------------------------------------------------

test_that("the worked example and its domain print exactly", {
  fixture <- snapshot_fixture()
  expect_snapshot(print(fixture$example))
  expect_snapshot(print(fixture$example$domain))
  expect_snapshot(print(abstract_domain(4, id = "snapshot:abstract:v1")))
})

test_that("frame specifications and a compiled frame print exactly", {
  fixture <- snapshot_fixture()
  expect_snapshot(print(whole_brain()))
  expect_snapshot(print(voxelwise()))
  expect_snapshot(print(searchlights(radius = 4)))
  expect_snapshot(print(regions(rep(c("a", "b"), length.out = 4))))
  expect_snapshot(print(fixture$example$frame))
})

test_that("the effect space prints exactly", {
  expect_snapshot(print(effect_space(c("face", "body", "house", "tool"),
    basis_id = "snapshot:conditions:v1")))
})

test_that("a relation, its fit, and a pairing print exactly", {
  fixture <- snapshot_fixture()
  expect_snapshot(print(fixture$relation))
  expect_snapshot(print(fixture$example$fit))
  expect_snapshot(print(fixture$pairing))
})

test_that("a geometry plan and a materialized geometry print exactly", {
  fixture <- snapshot_fixture()
  expect_snapshot(print(fixture$plan))
  expect_snapshot(print(fixture$geometry))
})

test_that("the compute policy prints exactly", {
  expect_snapshot(print(compute_policy()))
})

# Capabilities and refusals --------------------------------------------------

test_that("sampling capabilities print exactly, granted and withheld", {
  fixture <- snapshot_fixture()
  expect_snapshot(print(sampling_capabilities(fixture$plan, fixture$example$fit)))
  expect_snapshot(print(sampling_capabilities(fixture$plan)))
})

test_that("a capability refusal prints exactly", {
  domain <- abstract_domain(2, id = "snapshot:refusal:v1")
  relation <- relation(
    list(a = matrix(1:4, 2), b = matrix(2:5, 2)),
    effects = c("x", "y"), domain = domain
  )
  refusal <- catch_refusal(rdm_sampling_covariance(
    plan_geometry(relation, compile_frame(whole_brain(), domain),
      cross_partitions(relation, independence = "independent")),
    relation, target = "null", at = 1L
  ))
  expect_snapshot(print(refusal))
})

# Views ----------------------------------------------------------------------
#
# The printed body of each of these is a numeric preview, which is why only the
# one-line summary is pinned. `test-print-methods.R` holds the preview to the
# structural contract, and the numbers themselves are the subject of the view
# tests rather than of a printing test.

test_that("the view classes summarize themselves exactly", {
  fixture <- snapshot_fixture()
  expect_snapshot(cat(format(contrast_energy(fixture$geometry,
    fixture$example$contrast))))
  expect_snapshot(cat(format(rdm(fixture$geometry))))
  expect_snapshot(cat(format(rsa(fixture$plan,
    models = list(category = fixture$example$model_rdm)))))
  expect_snapshot(cat(format(rdm_sampling_covariance(fixture$plan,
    fixture$example$fit, target = "null", at = 1L))))
})
