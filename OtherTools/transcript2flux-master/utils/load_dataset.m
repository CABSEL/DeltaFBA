function data = load_dataset(name)
% Load one of the datasets.
% Options:
%   - 'ishii': Ishii dataset for E. coli
%   - 'ishii-protein': Ishii dataset (replace transcriptome with proteome)
%   - 'holm': Holm dataset for E. coli
%   - 'rintala': Rintala dataset for yeast
%   - 'rintala-red': Rintala dataset (reduced transcriptome)
%   - 'rintala-protein': Rintala dataset (replace transcriptome with proteome)
%
% Author: Daniel Machado, 2013

    ishii.folder = 'datasets/ishii/';
    ishii.transcriptomics = 'transcriptomics_nmol_per_gDW.csv';
    ishii.transcriptomics_std = 'transcriptomics_std_nmol_per_gDW.csv';
    ishii.proteomics = 'protein_nmol_per_gDW.csv';
    ishii.proteomics_std = 'protein_std_nmol_per_gDW.csv';
    ishii.fluxomics = 'fluxomics_mmol_per_gDW_per_h.csv';

    holm.folder = 'datasets/holm/';
    holm.transcriptomics = 'transcriptome.csv';
    holm.transcriptomics_std = 'transcriptomics_std.csv';
    holm.fluxomics = 'fluxes_mmol_per_gDW_per_h.csv';  

    rintala.folder = 'datasets/rintala/';
    rintala.transcriptomics = 'transcriptomics.csv';
    rintala.transcriptomics_std = 'transcriptomics_std.csv';
    rintala.transcriptomics_red = 'transcriptomics_red.csv';
    rintala.transcriptomics_std_red = 'transcriptomics_std_red.csv';
    rintala.proteomics = 'proteomics.csv';
    rintala.proteomics_std = 'proteomics_std.csv';
    rintala.fluxomics = 'fluxomics_mmol_per_gDW_per_h.csv';
    
    switch name
        case 'ishii'
            dataset = ishii;
        case 'ishii-protein'
            dataset = ishii;
            dataset.transcriptomics = ishii.proteomics;
            dataset.transcriptomics_std = ishii.proteomics_std;
        case 'holm'
            dataset = holm;       
        case 'rintala'
            dataset = rintala;
        case 'rintala-red'
            dataset = rintala;
            dataset.transcriptomics = rintala.transcriptomics_red;
            dataset.transcriptomics_std = rintala.transcriptomics_std_red;
        case 'rintala-protein'
            dataset = rintala;
            dataset.transcriptomics = rintala.proteomics;
            dataset.transcriptomics_std = rintala.proteomics_std;
    end
  
    [genes, conditions, transcriptomics] = load_transcript_data([dataset.folder dataset.transcriptomics]);
    [reactions, conditions2, fluxomics] = load_flux_data([dataset.folder dataset.fluxomics]);
    
    if ~all(strcmp(conditions, conditions2))
        disp('Warning: conditions do not match.');
    end
    
    data.name = name;
    data.conditions = conditions;
    data.genes = genes;
    data.transcriptomics = transcriptomics;
    data.reactions = reactions;
    data.fluxomics = fluxomics;
    data.knockouts = cell(size(conditions));

    switch name
        case 'ishii'
            data.knockouts = {'', '', '', '', '', 'b0756', 'b2388', 'b0688', ...
                'b4025', 'b3916', 'b1723', 'b4232', 'b2097', '', 'b0755', ...
                'b3612', 'b1854', 'b1676', 'b1702', 'b1852', 'b0767', 'b2029', ...
                'b3386', 'b2914', 'b4090', 'b2935', 'b2465', 'b2464', 'b0008'};

    end
    
    if isfield(dataset, 'transcriptomics_std')
        [genes2, conditions2, transcriptomics_std] = load_transcript_data([dataset.folder dataset.transcriptomics_std]);

        if ~all(strcmp(genes, genes2))
            disp('Warning: gene names do not match.');
        end

        if ~all(strcmp(conditions, conditions2))
            disp('Warning: conditions do not match.');
        end

        data.transcriptomics_std = transcriptomics_std;
    end
end

function [genes, conditions, transcriptomics] = load_transcript_data(file)
    table = importdata(file);
    conditions = table.textdata(1, 2:end);
    genes = table.textdata(2:end, 1);
    transcriptomics = table.data;
    transcriptomics(transcriptomics < 1e-12) = 0;
end

function [reactions, conditions, fluxomics] = load_flux_data(file)
    table = importdata(file);
    conditions = table.textdata(1, 2:end);
    reactions = table.textdata(2:end, 1);
    fluxomics = table.data;
    
    for i = 1:length(reactions)
        reactions{i} = strrep(reactions{i},'R_','');
        reactions{i} = strrep(reactions{i}, '_e_', '(e)');
    end
end