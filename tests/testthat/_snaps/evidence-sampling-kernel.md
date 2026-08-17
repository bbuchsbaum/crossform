# a sampling covariance in the general evidence basis renders

    Code
      print(fixture$covariance)
    Output
      <effect_sampling_covariance>
        basis:          evidence
        coordinates:    7
        labels:         d1, d2, d3 (+4 more)
        partitions:     4
        target:         point_estimate / nonnegative_plugin
        factors:        signal 7 x 7, xi 7 x 7
        residual noise: unrecorded
        storage:        exact factorized covariance
        signature:      sha256:<digest>...

---

    Code
      format(fixture$covariance)
    Output
      [1] "<effect_sampling_covariance: 7 coordinates, nonnegative_plugin>"

