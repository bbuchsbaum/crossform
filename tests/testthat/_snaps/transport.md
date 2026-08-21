# a transport prints its shape, semantics, sink and provenance

    Code
      print(transport)
    Output
      <effect_location_transport>
        nodes:      5 native -> 2 group + sink
        semantics:  budget
        sink:       mass 1.3 of 5 rows, 26.0% of territory
        provenance: anatomical (cross-fit: none)
        built:      fixed warp, atlas v1
        signature:  sha256:<digest>...

---

    Code
      format(transport)
    Output
      [1] "<effect_location_transport: 5 native -> 2 group + sink, budget, anatomical>"

---

    Code
      print(transport_fixture(semantics = "density", row_mass = c(1, 2, 3, 4, 5)))
    Output
      <effect_location_transport>
        nodes:      5 native -> 2 group + sink
        semantics:  density
        row mass:   declared, total 15
        sink:       mass 1.3 of 5 rows, 41.3% of territory
        provenance: anatomical (cross-fit: none)
        built:      fixed warp, atlas v1
        signature:  sha256:<digest>...

