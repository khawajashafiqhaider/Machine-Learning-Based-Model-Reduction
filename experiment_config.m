function cfg = experiment_config()
%EXPERIMENT_CONFIG Canonical settings shared by all MLCR experiments.

cfg.randomSeeds = 1:30;
cfg.time = (0:0.02:20).';
cfg.dt = cfg.time(2) - cfg.time(1);

cfg.reduction.order12 = 5;
cfg.reduction.orderISS = 10;
cfg.reduction.frequencyBand = [0.1, 10];
cfg.reduction.quadratureNodes = 80;

cfg.data.trainTrajectories = 60;
cfg.data.validationTrajectories = 20;
cfg.data.prbsTestTrajectories = 10;
cfg.data.windTestTrajectories = 10;

cfg.dynamic.delayLength = 5;
cfg.dynamic.hiddenUnits = [32, 32];
cfg.dynamic.maxEpochs = 500;
cfg.dynamic.validationPatience = 10;
cfg.dynamic.miniBatchSize = 32;
cfg.dynamic.initialLearnRate = 1e-3;

cfg.static.hiddenUnits = 48;
cfg.static.trainFunction = "trainlm";

cfg.autoencoder.encoderUnits = [64, 32, 16];
cfg.autoencoder.maxEpochs = 1000;
cfg.autoencoder.miniBatchSize = 64;
cfg.autoencoder.initialLearnRate = 1e-3;

cfg.outputRoot = fullfile(fileparts(mfilename("fullpath")), "results");
cfg.rawOutput = fullfile(cfg.outputRoot, "raw");
cfg.summaryOutput = fullfile(cfg.outputRoot, "summary");
end
