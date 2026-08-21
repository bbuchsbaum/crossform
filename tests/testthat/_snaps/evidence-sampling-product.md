# an RDM sampling covariance and its batch render for a reader

    Code
      print(single)
    Output
      <effect_sampling_covariance>
        basis:        rdm
        distances:    6
        measurement:  1
        partitions:   4 (dependent pair products)
        target:       null / fixed_zero
        metric:       fixed
        residual:     240 residual df, 2.997 effective dimensions (Wishart-corrected tr(Sigma^2))
        storage:      exact factorized covariance
        spatial law:  local marginal only

---

    Code
      format(single)
    Output
      [1] "<effect_sampling_covariance: 6 distances, fixed_zero, factorized>"

---

    Code
      print(batched)
    Output
      <effect_sampling_covariance_batch>
        basis:        rdm
        measurements: 3
        distances:    6
        execution:    batched / shared_pair_statistics
        shared residual statistics: yes

