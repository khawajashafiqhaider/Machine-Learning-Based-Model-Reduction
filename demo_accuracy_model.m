function metric=demo_accuracy_model(benchmark)
%DEMO_ACCURACY_MODEL Run a supplied improved LTI model; no neural toolbox needed.
if nargin<1,benchmark="official";end
assert(any(string(benchmark)==["official","nominal","constrained12"]));
root=fileparts(mfilename('fullpath'));addpath(genpath(root));cfg=experiment_config();
file=fullfile(root,'results','accuracy_confirmation','deployment_models',sprintf('%s_learned_seed_60001.mat',benchmark));
saved=load(file,'model');model=saved.model;
if string(benchmark)=="constrained12",m=project_constraints(benchmark12_model());else,m=load_iss_model();end
if string(benchmark)=="nominal",m=set_modal_damping_ratio(m,0.05);end
f=second_order_to_first_order(m.M,m.D,m.K,m.B,m.Cq,m.Cv);
u=generate_excitation('prbs',cfg.time,size(f.B,2),900001,1);
y=simulate_lti(f,cfg.time,u);yp=simulate_lti(model,cfg.time,u);
metric=trajectory_metrics(y,yp);disp(metric);
fprintf('Model order %d; maximum real pole %.6g\n',size(model.A,1),max(real(eig(model.A))));
end
