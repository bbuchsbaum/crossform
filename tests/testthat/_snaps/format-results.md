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
       measurement signed  coherent configuration    total coherence_fraction
                m1   -1.0 0.8357864     0.3357864 1.171573          0.7133883
                m2   -0.5 1.1715729     1.7573593 2.928932          0.4000000

---

    Code
      print(rdm_view)
    Output
      <effect_rdm_view>
        measurements: 2
       measurement    a - b
                m1 1.171573
                m2 2.928932

---

    Code
      print(rsa_view)
    Output
      <effect_rsa_view>
        measurements: 2
       measurement separation
                m1   1.171573
                m2   2.928932

---

    Code
      print(spectrum_view)
    Output
      <effect_spectrum_view>
        measurements: 2
       measurement    root1     root2
                m1 3.732051 0.2679492
                m2 8.674235 1.3257654

# crossnobis views print and coerce without losing signs

    Code
      print(view)
    Output
      <effect_crossnobis_view>
        measurements: 3
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

