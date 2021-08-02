function build_figures()
% Build all figures from the paper.
%
% Author: Daniel Machado, 2013

    DPI = '-r300';
    methods = {'pFBA', 'GIMME', 'iMAT', 'MADE', 'E-Flux', 'Lee-12', 'RELATCH', 'GX-FBA'};
    datasets = {'ishii', 'holm', 'rintala'};

    dataset_labels = {'{\it E. coli} (Ishii 2007)', '{\it E. coli} (Holm 2010)', 'Yeast (Rintala 2009)'};
    ymaxs = [3 3 3 3 3 3];
    build_error_boxplots_together(datasets, methods, dataset_labels, ymaxs, DPI)
 
    datasets2 = {'ishii','ishii-protein','rintala-red','rintala-protein'};
    dataset_labels2 = {'Ishii (transcript)', 'Ishii (protein)', 'Rintala (transcript)', 'Rintala (protein)'};
    ymaxs = [4 4 4 4 4 4 4 4];
    build_gene_vs_protein_plots(datasets2, methods(2:end), dataset_labels2, ymaxs, DPI)
  
    build_heatmaps('ishii', methods, 'sim_intra', [-5 5], DPI);
    build_heatmaps('holm', methods, 'sim_intra', [-10 10], DPI);
    build_heatmaps('rintala', methods, 'sim_intra', [-1 1], DPI);
 
    build_secretion_plots_ishii('WT_0.7h-1', methods, DPI);
    build_secretion_plots_holm(methods, DPI);
    build_secretion_plots_rintala(methods, DPI);
 
    plot_sensitivity_analysis('holm', [3 3 3 3 3 5], DPI)
    plot_sensitivity_analysis('rintala', [3 3 3 3 3 5], DPI)
  
    methods = { 'GIMME', 'iMAT', 'MADE', 'E-Flux', 'Lee-12',  'GX-FBA', 'RELATCH', 'RELATCH'};
    plot_robustness_analysis(methods, 'rintala', 'O2_0.0', [0.4 0.4 0.4 1 0.4 0.4 0.4 0], DPI)
end


function build_error_boxplots_together(datasets, methods, dataset_labels, ymaxs, dpi)
	
    experiment_types = {'sim_all', 'sim_intra' };    
    sublabels = {'a)', 'b)','c)','d)','e)','f)'};
    
    k = 0;
    for i = 1:length(experiment_types)
        for j = 1:length(datasets)
            data = load_error_results(datasets{j}, experiment_types{i}, methods);
            
            k = k+1; subplot(length(experiment_types),length(datasets),k);
            box_plotting(data, methods, ymaxs(k));
            ylabel('normalized error')
            title([sublabels{k} '        ' dataset_labels{j} '          '])
        end
    end
    
	set(gcf,'PaperUnits', 'points')
	set(gcf,'PaperPosition', [0 0 600 400])
    print('-dtiff', dpi, 'images/error_boxplots_all.tiff');
    close
end

function build_gene_vs_protein_plots(datasets, methods, dataset_labels, ymaxs, dpi)
	
    experiment_types = {'sim_all', 'sim_intra' };
    sublabels = {'a)', 'b)','c)','d)','e)','f)', 'g)', 'h)'};

    k = 0;
    for i = 1:length(experiment_types)
        for j = 1:length(datasets)
            data = load_error_results(datasets{j}, experiment_types{i}, methods);
            k = k+1; subplot(length(experiment_types),length(datasets),k);
            box_plotting(data, methods, ymaxs(k));
            ylabel('normalized error')
            title([sublabels{k} '        ' dataset_labels{j} '          '])
        end
    end
    
	set(gcf,'PaperUnits', 'points')
	set(gcf,'PaperPosition', [0 0 800 400])
    print('-dtiff', dpi, 'images/gene_vs_protein.tiff');
    close
end

function build_heatmaps(dataset, methods, experiment_type, clim, dpi)
	
    n = length(methods);
    order = [];
    for i = 1:n
        experiment = load_experiment(dataset, experiment_type, methods{i});
        
        if i == 1
            [experiment, order] = re_sort(experiment, []);
        else
            experiment = re_sort(experiment, order);
        end
        subplot(5, n, i);
        plot_error(experiment.error_all);
        title(methods{i}, 'FontWeight', 'bold')
        
        subplot(5, n, n*(1:4)+i);
        show_labels = i == 1;
        heatmap(experiment, show_labels, clim);
    end
    
    set(gcf,'PaperUnits', 'points')
    set(gcf,'PaperPosition', [0 0 1000 600])
    filename = sprintf('images/heatmaps_%s.tiff', dataset);
    print('-dtiff', dpi, filename);
    close
end

function build_secretion_plots_ishii(condition, methods, dpi)
    
    reaction_labels = {'Acetate', 'CO_2', 'Ethanol', 'Formate', ...
        'Lactate', 'Pyruvate', 'Succinate', 'Growth'};
    
    reaction_ids = {'EX_ac(e)';'EX_co2(e)';'EX_etoh(e)';'EX_for(e)';
        'EX_lac_D(e)';'EX_pyr(e)';'EX_succ(e)';'Ec_biomass_iAF1260_core_59p81M'};
    
    nrows = 2; ncols = 4;
    img_size = [0 0 800 500];
    
    TOL = 1e-6;
    
    expdata = load_dataset('ishii');
    cond_idx = strcmp(condition, expdata.conditions);
    [~, reaction_idx_data] = ismember(reaction_ids, expdata.reactions);
    reference = expdata.fluxomics(reaction_idx_data, cond_idx);
    
    data = load_results_for_condition('ishii', 'sim_secr', condition, methods);
    data(abs(data) < TOL) = 0;
    
    for i = 1:length(reaction_ids)
        subplot(nrows, ncols, i)
        barplot(data(i,:), reference(i), reaction_labels{i}, methods)
        set(gca,'OuterPosition', get(gca,'OuterPosition'))
    end
    
	set(gcf,'PaperUnits', 'points')
    set(gcf,'PaperPosition', img_size)
    filename = sprintf('images/secretion_%s.tiff', 'ishii');
    print('-dtiff', dpi, filename);
    close
    
end


function build_secretion_plots_holm(methods, dpi)
    
    reaction_ids = {'EX_ac(e)';'Ec_biomass_iAF1260_core_59p81M'};
    TOL = 1e-6;
    
    dataset = 'holm';
    expdata = load_dataset(dataset);
    
    condition = 'NOX';
    reaction_labels = {'Acetate (NADP oxidase)', 'Growth (NADP oxidase)'};
    
    cond_idx = strcmp(condition, expdata.conditions);
    [~, reaction_idx_data] = ismember(reaction_ids, expdata.reactions);
    reference = expdata.fluxomics(reaction_idx_data, cond_idx);
    data = load_results_for_condition(dataset, 'sim_secr', condition, methods);
    data(abs(data) < TOL) = 0;
    
    for i = 1:length(reaction_ids)
        subplot(1, 4, i)
        barplot(data(i,:), reference(i), reaction_labels{i}, methods)
        set(gca,'OuterPosition', get(gca,'OuterPosition'))
    end
    
    condition = 'ATP';
    reaction_labels = {'Acetate (ATPase)', 'Growth (ATPase)'};
    
    cond_idx = strcmp(condition, expdata.conditions);
    [~, reaction_idx_data] = ismember(reaction_ids, expdata.reactions);
    reference = expdata.fluxomics(reaction_idx_data, cond_idx);
    data = load_results_for_condition(dataset, 'sim_secr', condition, methods);
    data(abs(data) < TOL) = 0;
    
    for i = 1:length(reaction_ids)
        subplot(1, 4, 2+i)
        barplot(data(i,:), reference(i), reaction_labels{i}, methods)
        set(gca,'OuterPosition', get(gca,'OuterPosition'))
    end
    
    set(gcf,'PaperUnits', 'points')
    set(gcf,'PaperPosition', [0 0 800 200])
    filename = sprintf('images/secretion_%s.tiff', dataset);
    print('-dtiff', dpi, filename);
    close
    
end


function build_secretion_plots_rintala(methods, dpi)
    
    reaction_ids = {'ACxtO';'CBIOMASS';'ETHxtO';'GLxtO'};
    TOL = 1e-6;
    
    dataset = 'rintala';
    expdata = load_dataset(dataset);
    
    condition = 'O2_20.9';
    reaction_labels = {'Acetate (aerobic)', 'Growth (aerobic)', 'Ethanol (aerobic)', 'Glycerol (aerobic)'};
    
    cond_idx = strcmp(condition, expdata.conditions);
    [~, reaction_idx_data] = ismember(reaction_ids, expdata.reactions);
    reference = expdata.fluxomics(reaction_idx_data, cond_idx);
    data = load_results_for_condition(dataset, 'sim_secr', condition, methods);
    data(abs(data) < TOL) = 0;
    
    position = [1, 4, 2, 3];
    for i = 1:length(reaction_ids)
        subplot(2, 4, position(i))
        barplot(data(i,:), reference(i), reaction_labels{i}, methods)
        set(gca,'OuterPosition', get(gca,'OuterPosition'))
    end
    
    condition = 'O2_0.0';
    reaction_labels = {'Acetate (anaerobic)', 'Growth (anaerobic)', 'Ethanol (anaerobic)', 'Glycerol (anaerobic)'};
    
    cond_idx = strcmp(condition, expdata.conditions);
    [~, reaction_idx_data] = ismember(reaction_ids, expdata.reactions);
    reference = expdata.fluxomics(reaction_idx_data, cond_idx);
    data = load_results_for_condition(dataset, 'sim_secr', condition, methods);
    data(abs(data) < TOL) = 0;
    
    position = [5, 8, 6, 7];
    for i = 1:length(reaction_ids)
        subplot(2, 4, position(i))
        barplot(data(i,:), reference(i), reaction_labels{i}, methods)
        set(gca,'OuterPosition', get(gca,'OuterPosition'))
    end
    
    set(gcf,'PaperUnits', 'points')
    set(gcf,'PaperPosition', [0 0 800 500])
    filename = sprintf('images/secretion_%s.tiff', dataset);
    print('-dtiff', dpi, filename);
    close
    
end


function plot_sensitivity_analysis(dataset, ymaxs, dpi)
    methods = {'GIMME', 'GIMME', 'MADE', 'iMAT', 'iMAT', 'iMAT'};
    parameters = {'GIMME_LOWER_QUANTILE', 'OBJ_FRAC', 'OBJ_FRAC', 'IMAT_LOWER_QUANTILE', 'IMAT_UPPER_QUANTILE', 'IMAT_EPS'};
    labels = {'expression threshold', 'objective fraction', 'objective fraction', 'low expression threshold', 'high expression threshold', 'flux activation threshold'};
    
    for i = 1:6
        [values, data, scale] = load_sensitivity_analysis(dataset, methods{i}, parameters{i});
        subplot(2,3,i); plot_sensitivity(values, data, scale, ymaxs(i), labels{i}, 'Error')
        title(methods{i})
    end
    
    set(gcf,'PaperUnits', 'points')
    set(gcf,'PaperPosition', [0 0 600 300])
    filename = sprintf('images/sensitivity_%s.tiff', dataset);
    print('-dtiff', dpi, filename);
    close
end

function plot_robustness_analysis(methods, dataset, condition, ymaxs, dpi)
    
    for i = 1:length(methods)
        [alpha, data] = load_robustness_analysis(methods{i}, dataset, condition);
        subplot(2,4,i); plot_sensitivity(alpha, data, 'lin', ymaxs(i), 'noise', 'Error')
        title(methods{i})
    end
    
    set(gcf,'PaperUnits', 'points')
    set(gcf,'PaperPosition', [0 0 600 300])
    filename = sprintf('images/robustness_%s.tiff', dataset);
    print('-dtiff', dpi, filename);
    close
end

function data = load_error_results(dataset, experiment_type, methods)
    data = cell(1, length(methods));
    for i = 1:length(methods)
        filename = sprintf('results/%s_%s_%s.mat', methods{i}, dataset, experiment_type);
        load(filename);
        data{i} = experiment.error_all(experiment.status_all == 1);
    end
end

function experiment = load_experiment(dataset, experiment_type, method)
    filename = sprintf('results/%s_%s_%s.mat', method, dataset, experiment_type);
    load(filename);
end

function [values, data, scale] = load_sensitivity_analysis(dataset, method, parameter)
    filename = sprintf('results/sensitivity_%s_%s_%s.mat', method, parameter, dataset);
    load(filename);
    values = analysis.param_values;
    data = analysis.error_data;
    scale = analysis.scale;
end

function [alpha, data] = load_robustness_analysis(method, dataset, condition)
    filename = sprintf('results/robustness_%s_%s_%s.mat', method, dataset, condition);
    load(filename);
    alpha = analysis.alpha;
    data = analysis.error_data;
end


function data = load_results_for_condition(dataset, experiment_type, condition, methods)
    data = [];
    for i = 1:length(methods)
        filename = sprintf('results/%s_%s_%s.mat', methods{i}, dataset, experiment_type);
        load(filename);
        cond_idx = strcmp(condition, experiment.conditions);
        data = [data experiment.fluxes_sim_all{cond_idx}];
    end
end


function [experiment, order] = re_sort(experiment, order)

    experiment.error_all(experiment.status_all == 0) = nan;

    if isempty(order)
        [~, order] = sort(experiment.error_all);
    end
    
    experiment.error_all = experiment.error_all(order);
    experiment.status_all = experiment.status_all(order);
    experiment.conditions = experiment.conditions(order);
    experiment.fluxes_exp_all = experiment.fluxes_exp_all(order);
    experiment.fluxes_sim_all = experiment.fluxes_sim_all(order);
end
