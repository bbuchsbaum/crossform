# first-journey objects print compact summaries

    Code
      print(geometry)
    Output
      <effect_geometry>
        effects:      2 x 2
        measurements: 2
        components:   total + coherent + configuration
        storage:      memory
        estimate:     signed cross-generalized form; PSD not assumed

---

    Code
      print(query)
    Output
      <effect_view>
        measurements: 2
       measurement view1
                m1    -2
                m2    -2

---

    Code
      print(contrast_view)
    Output
      <effect_contrast_view>
        measurements: 2
        contrast:     a 1, b -1
       measurement signed coherent configuration total coherence_fraction
                m1   -1.0   0.8358        0.3358 1.172             0.7134
                m2   -0.5   1.1716        1.7574 2.929             0.4000
        coherence_fraction: 2 of 2 valid; NA where coherent and configuration are
          not a nonnegative partition

---

    Code
      print(rdm_view)
    Output
      <effect_rdm_view>
        measurements: 2
       measurement a - b
                m1 1.172
                m2 2.929

---

    Code
      print(rsa_view)
    Output
      <effect_rsa_view>
        measurements: 2
        models:       separation
        intercept:    omitted
        pairs:        1 distance over 2 effects
        estimate:     OLS coefficients on the total component
       measurement separation
                m1      1.172
                m2      2.929
        next:         as.data.frame(x), plot(x, terms = ...)

---

    Code
      print(spectrum_view)
    Output
      <effect_spectrum_view>
        measurements: 2
       measurement root1  root2
                m1 3.732 0.2679
                m2 8.674 1.3258

# crossnobis views print and coerce without losing signs

    Code
      print(view)
    Output
      <effect_crossnobis_view>
        measurements: 3
        contrast:     a 1, b -1
       measurement crossnobis
                 1        1.1
                 2        0.9
                 3        0.0

# relation fits print capabilities rather than nested sources

    Code
      print(fit)
    Output
      <effect_relation_fit>
        effects:      1
        features:     3
        partitions:   1
        residuals:    1/1 partitions
        covariance:   1/1 partitions

