%% MAIN FILE:
%
% Use this to run all the tests in the paper.
% Different tests can be run selectively.
% For convenience all results are stored as .mat files.
%
% Author: Daniel Machado, 2013

addpath('tests')
addpath('utils')
addpath('wrappers')
addpath('plotting')

%% PROBLEM SETUP

% select dataset and organism (see load_dataset.m for details)
%DATASET = 'ishii'; ORGANISM = 'ecoli';
DATASET = 'holm'; ORGANISM = 'ecoli';
%DATASET = 'rintala'; ORGANISM = 'yeast';

% select experiment type
experiment_type = 'sim_all'; % simulate all fluxes 
%experiment_type = 'sim_intra'; % simulate intracellular fluxes 
%experiment_type = 'sim_secr'; % simulate secretion and growth rates

COBRA_SOLVER = 'gurobi5';

%initialize cobra
changeCobraSolver(COBRA_SOLVER, 'all');

%initialize cmpi (requires MADE package; see README)
cmpi.init()

%load model
model = load_model(ORGANISM);

%load data set
dataset = load_dataset(DATASET);

options.experiment_type = experiment_type;
options.reestimate_data = true; % fit experimental data to model


%% BENCHMARKING

methods = {'pFBA', 'GIMME', 'iMAT', 'E-Flux', 'Lee-12', 'RELATCH', 'MADE', 'GX-FBA'};

for i = 1:length(methods)
	benchmark_method(methods{i}, model, dataset, options);
end


%% SENSITIVITY ANALYSIS

conf = cell(1,6);
conf{1} = {'GIMME', 'OBJ_FRAC', 0, 1, 'lin'};
conf{2} = {'GIMME', 'GIMME_LOWER_QUANTILE', 0, 1, 'lin'}; 
conf{3} = {'iMAT', 'IMAT_LOWER_QUANTILE', 0, 0.75, 'lin'}; 
conf{4} = {'iMAT', 'IMAT_UPPER_QUANTILE', 0.25, 1, 'lin'}; 
conf{5} = {'iMAT', 'IMAT_EPS', -2, 2, 'log'}; 
conf{6} = {'MADE', 'OBJ_FRAC', 0, 1, 'lin'}; 
 
points = 100; 
 
for i = 1:6
    [method, parameter, min_val, max_val, scale] = deal(conf{i}{:}); 
    sensitivity_analysis(model, dataset, method, parameter, min_val, max_val, points, scale, options);
end


%% ROBUSTNESS ANALYSIS

methods = {'GIMME', 'iMAT', 'E-Flux', 'Lee-12', 'RELATCH', 'MADE', 'GX-FBA'};
steps = 2; points = 2;

for i = 1:length(methods)
     robustness_analysis(model, dataset, methods{i}, dataset.conditions{2}, dataset.conditions{1}, steps, points, options);
end

%% PLOTTING

build_figures()
