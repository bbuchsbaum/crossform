# an exactly zero form moves no mass and earns no share

    Code
      print(latent)
    Output
      <effect_latent_geometry>
        measurements: 1
        component:    total
        projection:   psd_projection (eigenvalue truncation at zero)
        moved_mass:   none, no source root was negative
        n_eff:        undefined; no measurement kept positive mass
        reading:      latent descriptive layer; not for inference
       measurement n_eff moved_mass moved_share root1 root2 root3
              zero    NA          0           0     0     0     0
        next:         x$cumulative, x$projection

# the latent layer prints what it is and what the clipping cost

    Code
      print(latent)
    Output
      <effect_latent_geometry>
        measurements: 4
        component:    total
        projection:   psd_projection (eigenvalue truncation at zero)
        moved_mass:   7 moved, 3 of 4 measurements clipped, max share 1,
                      1 moved none
        n_eff:        median 1.98 (range 1.6 to 2.65), 1 measurement masked
        reading:      latent descriptive layer; not for inference
       measurement n_eff moved_mass moved_share root1 root2 root3 root4
               psd 2.647          0      0.0000     4     2   1.0   0.5
        indefinite 1.600          2      0.3333     3     1   0.0   0.0
          negative    NA          4      1.0000     0     0   0.0   0.0
             block 1.976          1      0.1818     3     1   0.5   0.0
        next:         x$cumulative, x$projection

---

    Code
      format(latent)
    Output
      [1] "<effect_latent_geometry: 4 measurements, 4 nonnegative roots>"

---

    Code
      print(latent$projection)
    Output
      <effect_latent_projection_receipt>
        method:     psd_projection
        operator:   eigenvalue truncation at zero
        component:  total
        clipped:    3 of 4 measurements
        moved_mass: 7 total, max share 1
        masked:     1 measurement
        source:     plan-sha256:abc

---

    Code
      format(latent$projection)
    Output
      [1] "<effect_latent_projection_receipt: psd_projection, 3 of 4 measurements clipped>"

---

    Code
      print(clean)
    Output
      <effect_latent_geometry>
        measurements: 2
        component:    total
        projection:   psd_projection (eigenvalue truncation at zero)
        moved_mass:   none, no source root was negative
        n_eff:        median 1.8 (range 1.6 to 2)
        reading:      latent descriptive layer; not for inference
       measurement n_eff moved_mass moved_share root1 root2
                m1   1.6          0           0     3     1
                m2   2.0          0           0     2     2
        next:         x$cumulative, x$projection

