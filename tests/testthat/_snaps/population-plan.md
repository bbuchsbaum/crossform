# a population plan prints its estimand-bearing choices

    Code
      print(plan)
    Output
      <effect_population_plan>
        subjects:      3 (s01, s02, s03)
        group nodes:   2 + sink
        sink:          empty in every subject (full native coverage)
        transport:     budget, anatomical
        model:         ~age -> 2 columns, rank 2
        normalization: unit_budget
        fit:           OLS (subject-constant weights), transport then fit
        estimand:      population-sha256:<digest>...
        signature:     sha256:<digest>...

# the printed sink line reports partial coverage

    Code
      print(plan)
    Output
      <effect_population_plan>
        subjects:      2 (s01, s02)
        group nodes:   2 + sink
        sink:          present in 2 of 2 subjects, worst 50.0% of territory
        transport:     budget, anatomical
        model:         ~1 -> 1 column, rank 1
        normalization: none
        fit:           OLS (subject-constant weights), transport then fit
        estimand:      population-sha256:<digest>...
        signature:     sha256:<digest>...

