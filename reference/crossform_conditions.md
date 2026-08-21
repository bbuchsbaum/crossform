# Conditions raised by crossform

Every failure crossform raises is a classed condition, so a caller can
branch on the *cause* rather than matching message prose. There are four
classes, and each answers a different question about what went wrong.

## Value

These are condition classes, not functions; there is nothing to call.
[`catch_refusal()`](https://bbuchsbaum.github.io/crossform/reference/catch_refusal.md)
captures a refusal as a value.

## The classes

- `effect_input_error`:

  An argument has the wrong type, shape, or value — a character vector
  where one string was expected, a matrix holding `NA`, a negative
  count. The caller fixes it by passing something else. May carry `$arg`
  (the argument name), `$received` (a short description of what
  arrived), and `$expected` (what was wanted); see "The three optional
  fields" below.

- `effect_contract_error`:

  Two objects disagree, or an object disagrees with its own recorded
  identity: a result whose receipt names a different scientific plan, a
  frame compiled against a different neural domain, a signature that no
  longer matches the fields it summarizes. Each object may be
  individually well formed; they simply do not belong together. The fix
  is to pass objects that were built from the same declarations, not to
  repair either one.

- `effect_invariant_error`:

  The package computed something impossible. Nothing the caller passed
  can explain it, so the message asks for a bug report rather than a
  change of input.

- `effect_capability_refusal`:

  Not an error in the input: the requested interpretation cannot be
  earned from the objects supplied — for example, analytic standard
  errors from a fit that did not retain the residual channel they
  require. Refusals carry `$capability`, `$namespace`, all unmet
  `$reasons`, and concrete `$remedies`. See
  [`catch_refusal()`](https://bbuchsbaum.github.io/crossform/reference/catch_refusal.md).

The first three inherit from `effect_error`, so
`tryCatch(expr, effect_error = ...)` catches any of them while letting a
refusal through. All four inherit from `error`, so ordinary
[`try()`](https://rdrr.io/r/base/try.html) and
`tryCatch(expr, error = ...)` behave exactly as before.

## The three optional fields

`$arg`, `$received`, and `$expected` are a convenience, not a contract.
They are filled in by the shared argument guards (the `.check_*` family
behind most exported entry points) and by the hand-written checks at the
entry points that know the three values, which covers the argument
errors a caller provokes in ordinary use. Every other site — internal
checks whose failure is about a whole record rather than one argument,
and checks whose "expected" is a paragraph rather than a phrase — leaves
them `NULL` and says everything in
[`conditionMessage()`](https://rdrr.io/r/base/conditions.html).

So branch on the **class** first, which is always there, and treat the
fields as something to report when present:

    if (!is.null(condition$arg)) {
      message("bad argument `", condition$arg, "`: ", condition$received)
    }

The same three fields exist on `effect_contract_error` and
`effect_invariant_error` under the same rule; a contract error that
compares two identities typically fills `$received` and `$expected` with
the two shortened signatures.

## See also

[`catch_refusal()`](https://bbuchsbaum.github.io/crossform/reference/catch_refusal.md)

Other conditions:
[`catch_refusal()`](https://bbuchsbaum.github.io/crossform/reference/catch_refusal.md)

## Examples

``` r
domain <- abstract_domain(2, id = "conditions-example")
relation <- relation(
  list(a = matrix(1:4, 2), b = matrix(2:5, 2)),
  effects = c("x", "y"), domain = domain
)
plan <- plan_geometry(
  relation, compile_frame(whole_brain(), domain),
  cross_partitions(relation, independence = "independent")
)

# Branch on the class. Three weights for two effects is an input error.
tryCatch(
  contrast_energy(plan, c(1, -1, 0)),
  effect_input_error = function(e) conditionMessage(e),
  effect_contract_error = function(e) "these objects do not belong together"
)
#> [1] "`weights` has 3 values but the relation declares 2 effects (`x`, `y`). Supply one weight per effect, or name the weights to have them aligned for you."

# A guard-raised error additionally names the argument and the value.
guarded <- tryCatch(
  effect_space(c("x", "y"), basis_id = 42),
  effect_input_error = function(e) e
)
guarded$arg
#> [1] "basis_id"
guarded$received
#> [1] "42"
guarded$expected
#> [1] "one nonempty identifier"

# Not every site fills them, so test before you use them. Here the whole
# explanation is in the message.
bare <- tryCatch(abstract_domain(2, id = ""), effect_input_error = identity)
is.null(bare$arg)
#> [1] TRUE
conditionMessage(bare)
#> [1] "Domain `id` must be one nonempty identifier."

# `effect_error` catches input, contract, and invariant failures alike,
# while letting a capability refusal through.
tryCatch(abstract_domain(-1), effect_error = function(e) "caught")
#> [1] "caught"
```
