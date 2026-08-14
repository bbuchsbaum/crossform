# Mission of `effectagram`

Status: working mission
Date: 2026-08-13

> Our mission is to build a scientifically explicit and maintainable successor to `rMVPA` by organizing task-fMRI analysis around partitioned brain–experiment relations and effect geometry. We will replace method-specific pipelines with a small, law-tested core; extend it through disciplined learning, inference, domain, and workflow modules; and make every reported value traceable to its spatial measurement, generalization requirement, units, assumptions, and provenance.

We will pursue that mission in six ways:

1. **Prove the foundation.** Every optimized computation will be checked against an independent reference implementation and the algebraic laws that define its meaning.
2. **Improve interpretation.** The system will retain geometry before scalarization, distinguish coherent regional effects from spatial configuration, and make cross-partition generalization explicit. It will keep effect, generalization, and calibration distinct: an estimated evidence value and its sampling law must remain separately identifiable, and calibration unsupported by the available uncertainty information must be refused rather than improvised.
3. **Recover the important work, not the old architecture.** We will replace `rMVPA` user journeys through clearer common contracts. We will not reproduce every class, engine, callback, default, or historical feature.
4. **Grow without another monolith.** Adaptive methods, population inference, neuroimaging domains, calibration, and user-facing workflows will connect through narrow interfaces rather than enter the numerical core as parallel method families.
5. **Earn the transition with evidence.** Matched-estimand tests, simulations, public-data analyses, migration guides, and honest scope statements will determine when an effectagram workflow is ready to replace its `rMVPA` counterpart.
6. **Make the principled path practical.** We will use algebraic elimination,
   factorization, sparsity, query fusion, streaming, and matrix-free action to
   remove avoidable work. Every fast lowering must preserve the declared
   estimand under an explicit numerical contract; approximations are separate,
   named scientific choices rather than hidden performance modes.

Progress will be measured by the number and importance of research tasks that become easier to specify, verify, interpret, and extend—not by the number of exported functions. The immediate work is to freeze and prove the core contracts. The continuing work is to build the better system those contracts make possible.
