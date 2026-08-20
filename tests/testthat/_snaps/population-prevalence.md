# the printed record is stable

    Code
      print(prevalence)
    Output
      <effect_population_prevalence>
        layer:        latent descriptive (fractions over participants)
        ledger:       transported_total (component "total")
        participants: 4 (s01, s02, s03, s04)
        group nodes:  3 (sink excluded)
        queries:      face-house, face-tool
        threshold:    value > 0 (ledger units, exact)
        sign:         median 0.75 (range 0.5 to 1)
        alignment:    median 1 (range 0.5 to 1), leave-one-out reference
        readout:      2 queries; bank Gram off identity by 3, not Frobenius
        coverage:     minimum 2 of 4 participants contributing; no floor declared
        frame:        undeclared, conservative
        transport:    budget, anatomical, cross-fit not declared
        estimand:     population-sha256:1b0f2bba5cf2...
        reading:      latent descriptive layer; not for inference
        a cell at which nothing reproduces reports a fraction near 0.5, not near
          0: thresholding a signed crossvalidated estimate keeps the sign and
          discards the magnitude, and the sign is the noisy part.
        descriptive only. No standard error, interval or p-value is attached to a
          count of participants, and none follows from one;
          population_uncertainty() is the separate inferential layer.
        next:         as.data.frame(x, measure), x$coverage, x$sign$resolved

---

    Code
      format(prevalence)
    Output
      [1] "<effect_population_prevalence: transported_total, 4 participants, 3 group nodes, threshold 0>"

