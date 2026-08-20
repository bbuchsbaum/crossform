# the worked example and its domain print exactly

    Code
      print(fixture$example)
    Output
      <effect_example_effects>
        fit:       <effect_relation_fit: 4 measurements, 4 effects, 280 features>
        domain:    <effect_domain: crossform:generated-fmri-example, 280 volume...
        frame:     <effect_frame: 280 nodes, additive_diagonal, searchlights>
        contrast:  4 weights (0.5, 0.5, -0.5, -0.5)
        model_rdm: 4 x 4 model dissimilarities
        planted:   19 pattern voxels + 19 mean-shift voxels
        truth:     104 signal measurements, seed 20260814
        next:      plan_geometry(example$fit$relation, example$frame, pairing)

---

    Code
      print(fixture$example$domain)
    Output
      <effect_domain>
        id:          crossform:generated-fmri-example
        kind:        volume
        features:    280
        feature ids: 1, 2, 3, 4 (+276 more)
        coordinates: 280 x 3 (mm, mm, mm)
        geometry:    sha256:be038e4744bd...
        next:        compile_frame(whole_brain(), domain) to place a frame on it

---

    Code
      print(abstract_domain(4, id = "snapshot:abstract:v1"))
    Output
      <effect_domain>
        id:          snapshot:abstract:v1
        kind:        abstract
        features:    4
        feature ids: 1, 2, 3, 4
        coordinates: none (abstract)
        geometry:    sha256:10d104c2adb5...
        next:        compile_frame(whole_brain(), domain) to place a frame on it

# frame specifications and a compiled frame print exactly

    Code
      print(whole_brain())
    Output
      <effect_frame_spec>
        kind:          whole_brain
        normalization: local
        state:         unplaced; call compile_frame(spec, domain) to bind a domain

---

    Code
      print(voxelwise())
    Output
      <effect_frame_spec>
        kind:          voxels
        normalization: conservative
        state:         unplaced; call compile_frame(spec, domain) to bind a domain

---

    Code
      print(searchlights(radius = 4))
    Output
      <effect_frame_spec>
        kind:          searchlights
        radius:        4
        normalization: local
        state:         unplaced; call compile_frame(spec, domain) to bind a domain

---

    Code
      print(regions(rep(c("a", "b"), length.out = 4)))
    Output
      <effect_frame_spec>
        kind:          regions
        normalization: local
        state:         unplaced; call compile_frame(spec, domain) to bind a domain

---

    Code
      print(fixture$example$frame)
    Output
      <effect_frame>
        representation: additive_diagonal
        nodes:          280
        features:       280
        specification:  searchlights, radius 4 mm
        normalization:  local
        weights:        280 x 280, 1698 stored
        domain:         crossform:generated-fmri-example (volume)
        fixed:          yes; locally estimated no
        next:           plan_geometry(relation, frame, pairing)

# the effect space prints exactly

    Code
      print(effect_space(c("face", "body", "house", "tool"), basis_id = "snapshot:conditions:v1"))
    Output
      <effect_space>
        coordinates: 4 (face, body, house, tool)
        basis:       snapshot:conditions:v1
        units:       arbitrary
        scale:       1
        signature:   sha256:02ed17e017d0...

# a relation, its fit, and a pairing print exactly

    Code
      print(fixture$relation)
    Output
      <effect_relation>
        effects:    4 (face, body, house, tool)
        partitions: 4 (run1, run2, run3, run4)
        features:   280
        domain:     crossform:generated-fmri-example
        basis:      generated-condition-means:v1
        sources:    4 matrix (unread)
        extraction: 4 extractors, linear_model
        state:      lazy; sources are read only by neural feature block
        next:       plan_geometry(relation, frame, pairing)

---

    Code
      print(fixture$example$fit)
    Output
      <effect_relation_fit>
        effects:      4
        features:     280
        partitions:   4
        residuals:    4/4 partitions
        covariance:   4/4 partitions

---

    Code
      print(fixture$pairing)
    Output
      <effect_pairing>
        pairs:        6
        left:         run1, run2, run3
        right:        run2, run3, run4
        weights:      equal (0.1667)
        directed:     no
        self pairs:   forbid
        independence: independent
        generalizes:  axis undeclared
        estimate:     cross_generalized

# a geometry plan and a materialized geometry print exactly

    Code
      print(fixture$plan)
    Output
      <effect_geometry_plan>
        effects:      4 x 4
        measurements: 280
        features:     280
        metric:       implicit identity
        generalizes:  6 partition pairs (axis undeclared), endpoints independent
        execution:    query-first, in memory
        state:        nothing computed yet
        next:         contrast_energy(plan, weights), rdm(plan), rsa(plan)

---

    Code
      print(fixture$geometry)
    Output
      <effect_geometry>
        effects:      4 x 4
        measurements: 280
        components:   total + coherent + configuration
        storage:      memory
        estimate:     signed cross-generalized form; PSD not assumed

# the compute policy prints exactly

    Code
      print(compute_policy())
    Output
      <effect_compute_policy>
        workers:       1
        backend:       sequential
        feature block: chosen by the kernel
        workspace:     unbounded

# sampling capabilities print exactly, granted and withheld

    Code
      print(sampling_capabilities(fixture$plan, fixture$example$fit))
    Output
      <effect_sampling_capabilities>
        analytic sampling law: available 
        metric: fixed | partitions: equal | error channel: relation_fit 

---

    Code
      print(sampling_capabilities(fixture$plan))
    Output
      <effect_sampling_capabilities>
        analytic sampling law: unavailable 
        metric: fixed | partitions: equal | error channel: absent 
        unmet requirements:
        * missing_error_channel - this evidence plan has only a precomputed relation and no error channel. Refit raw observations with `lm_relation_fit()` or supply a validated, identity-bound external error channel; beta matrices alone cannot recover residual uncertainty 
            remedy: Refit raw observations with `lm_relation_fit()`. 
        * sampling_axis_missing_or_inconsistent - no single sampling axis is declared, or the declared axis conflicts with the error channel's recorded sampling unit 
            remedy: Declare one `sampling_axis` that matches the error channel's sampling unit. 
        note: requirements that describe the error channel itself cannot be
              evaluated until one exists, and are not listed.

# a capability refusal prints exactly

    Code
      print(refusal)
    Output
      <effect_capability_refusal>
        capability:  sampling_covariance
        namespace:   evidence_sampling
        reasons:
          - missing_error_channel
          - sampling_axis_missing_or_inconsistent
        remedies:
          - Refit raw observations with `lm_relation_fit()`.
          - Declare one `sampling_axis` that matches the error channel's sampling
            unit.
        state:       refused; no partial result was produced

# the view classes summarize themselves exactly

    Code
      cat(format(contrast_energy(fixture$geometry, fixture$example$contrast)))
    Output
      <effect_contrast_view: 280 measurements, signed + energy decomposition>

---

    Code
      cat(format(rdm(fixture$geometry)))
    Output
      <effect_rdm_view: 280 measurements, 6 distances>

---

    Code
      cat(format(rsa(fixture$plan, models = list(category = fixture$example$model_rdm))))
    Output
      <effect_rsa_view: 280 measurements, 2 coefficients>

---

    Code
      cat(format(rdm_sampling_covariance(fixture$plan, fixture$example$fit, target = "null",
      at = 1L)))
    Output
      <effect_sampling_covariance: 6 distances, fixed_zero, factorized>

