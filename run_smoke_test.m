function report = run_smoke_test()
%RUN_SMOKE_TEST Verify model conversion, constraint handling, and reductions.

root = fileparts(mfilename("fullpath"));
addpath(genpath(root));
cfg = experiment_config();

benchmark = benchmark12_model();
projected = project_constraints(benchmark);
system12 = second_order_to_first_order(projected.M, projected.D, projected.K, projected.B, projected.Cq, projected.Cv);
bt12 = sobt_reduce(projected, cfg.reduction.order12);
flbt12 = sobt_reduce(projected, cfg.reduction.order12, ...
    FrequencyBand=cfg.reduction.frequencyBand, QuadratureNodes=cfg.reduction.quadratureNodes);

iss = load_iss_model(fullfile(root, "iss.mat"));
systemISS = second_order_to_first_order(iss.M, iss.D, iss.K, iss.B, iss.Cq, iss.Cv);
btISS = sobt_reduce(iss, cfg.reduction.orderISS);
flbtISS = sobt_reduce(iss, cfg.reduction.orderISS, ...
    FrequencyBand=cfg.reduction.frequencyBand, QuadratureNodes=cfg.reduction.quadratureNodes);

factorUpper = ds_gauss("C", system12.E, system12.A, system12.B, ...
    cfg.reduction.frequencyBand(1), cfg.reduction.frequencyBand(2), cfg.reduction.quadratureNodes);
factorLower = ds_gauss("c", system12.E, system12.A, system12.B, ...
    cfg.reduction.frequencyBand(1), cfg.reduction.frequencyBand(2), cfg.reduction.quadratureNodes);
observationUpper = ds_gauss("O", system12.E, system12.A, system12.C, ...
    cfg.reduction.frequencyBand(1), cfg.reduction.frequencyBand(2), cfg.reduction.quadratureNodes);
observationLower = ds_gauss("o", system12.E, system12.A, system12.C, ...
    cfg.reduction.frequencyBand(1), cfg.reduction.frequencyBand(2), cfg.reduction.quadratureNodes);
scalarDecay = 2;
scalarFactor = ds_gauss("C", 1, -scalarDecay, 1, ...
    cfg.reduction.frequencyBand(1), cfg.reduction.frequencyBand(2), cfg.reduction.quadratureNodes);
scalarExactGramian = (atan(cfg.reduction.frequencyBand(2) / scalarDecay) ...
    - atan(cfg.reduction.frequencyBand(1) / scalarDecay)) / (pi * scalarDecay);

frequencies = logspace(log10(cfg.reduction.frequencyBand(1)), log10(cfg.reduction.frequencyBand(2)), 200);
report.constraintResidual = projected.constraintResidual;
report.projectorIdempotenceResidual = projected.projectorIdempotenceResidual;
report.rankFullProjector = rank(projected.theta);
report.rankFirstOrderProjectedDescriptor = 2 * rank(projected.theta);
report.orderAfterConstraintElimination = system12.order;
report.fullMassRank = rank(benchmark.M);
report.eliminatedMassRank = rank(projected.M);
report.eliminatedMassCondition = projected.massCondition;
report.controllabilitySelectorResidual = norm(factorUpper - factorLower, "fro");
report.observabilitySelectorResidual = norm(observationUpper - observationLower, "fro");
report.gramianScalingResidual = abs(scalarFactor * scalarFactor.' - scalarExactGramian);
report.bt12MassBiorthogonalityResidual = bt12.massBiorthogonalityResidual;
report.flbt12MassBiorthogonalityResidual = flbt12.massBiorthogonalityResidual;
report.btISSMassBiorthogonalityResidual = btISS.massBiorthogonalityResidual;
report.flbtISSMassBiorthogonalityResidual = flbtISS.massBiorthogonalityResidual;
report.bt12Stable = all(real(eig(bt12.A)) < 0);
report.flbt12Stable = all(real(eig(flbt12.A)) < 0);
report.btISSStable = all(real(eig(btISS.A)) < 0);
report.flbtISSStable = all(real(eig(flbtISS.A)) < 0);
report.bt12Frequency = frequency_metrics(system12, bt12, frequencies);
report.flbt12Frequency = frequency_metrics(system12, flbt12, frequencies);
report.btISSFrequency = frequency_metrics(systemISS, btISS, frequencies);
report.flbtISSFrequency = frequency_metrics(systemISS, flbtISS, frequencies);

assert(report.constraintResidual < 1e-10, "Constraint projector residual is too large.");
assert(report.projectorIdempotenceResidual < 1e-10, "Constraint projector is not idempotent.");
assert(report.rankFullProjector == 11, "Unexpected constrained-position rank.");
assert(report.orderAfterConstraintElimination == 22, "Unexpected eliminated first-order model size.");
assert(report.eliminatedMassRank == 11, "Constraint-eliminated mass matrix is singular.");
assert(report.controllabilitySelectorResidual == 0, "C/c selectors produced different factors.");
assert(report.observabilitySelectorResidual == 0, "O/o selectors produced different factors.");
assert(report.gramianScalingResidual < 1e-12, "Frequency-limited Gramian 1/pi scaling is incorrect.");
assert(report.bt12MassBiorthogonalityResidual < 1e-8, "12-coordinate SOBT projection is not mass-biorthogonal.");
assert(report.flbt12MassBiorthogonalityResidual < 1e-8, "12-coordinate FL-SOBT projection is not mass-biorthogonal.");
assert(report.btISSMassBiorthogonalityResidual < 1e-8, "ISS SOBT projection is not mass-biorthogonal.");
assert(report.flbtISSMassBiorthogonalityResidual < 1e-8, "ISS FL-SOBT projection is not mass-biorthogonal.");
assert(report.bt12Stable && report.flbt12Stable && report.btISSStable && report.flbtISSStable, ...
    "One or more reduced models is unstable.");

disp(report)
end
