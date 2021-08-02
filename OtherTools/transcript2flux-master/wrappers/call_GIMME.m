function fluxes = call_GIMME(model, gene_names, gene_exp, obj_frac, threshold)
% Implements GIMME as defined in [Becker and Palsson, PLoS Comp Bio, 2008].
% Adapted from the implementation provided in the cobra toolbox.
%
% INPUTS
%       model - cobra model
%       gene_names - gene ids
%       gene_exp - gene expression
%       obj_frac - minimum fraction of the biological objective
%       threshold - gene expression threshold
%
% OUTPUTS
%       fluxes - flux distribution
%
% Author: Daniel Machado, 2013 


    expressionCol = gene_to_reaction_levels(model, gene_names, gene_exp, @min, @max);
    expressionCol(isnan(expressionCol)) = -1;

    objectiveCol = [find(model.c) obj_frac];

    fluxes = solveGimme(model, objectiveCol, expressionCol, threshold);

end


function fluxes = solveGimme(model,objectiveCol,expressionCol,cutoff)
% Implementation from the cobra toolbox (createTissueSpecificModel.m)

    nRxns = size(model.S,2);

    %first make model irreversible
    [modelIrrev,matchRev,rev2irrev,irrev2rev] = convertToIrreversible(model);

    nExpressionCol = size(expressionCol,1);
    if (nExpressionCol < nRxns)
        display('Warning: Fewer expression data inputs than reactions');
        expressionCol(nExpressionCol+1:nRxns,:) = zeros(nRxns-nExpressionCol, size(expressionCol,2));
    end

    nIrrevRxns = size(irrev2rev,1);
    expressionColIrrev = zeros(nIrrevRxns,1);
    for i=1:nIrrevRxns
    %     objectiveColIrrev(i,:) = objectiveCol(irrev2rev(i,1),:);
        expressionColIrrev(i,1) = expressionCol(irrev2rev(i,1),1);
    end

    nObjectives = size(objectiveCol,1);
    for i=1:nObjectives
        objectiveColIrrev(i,:) = [rev2irrev{objectiveCol(i,1),1}(1,1) objectiveCol(i,2)];
    end

    %Solve initially to get max for each objective
    for i=1:size(objectiveCol)
        %define parameters for initial solution
        modelIrrev.c=zeros(nIrrevRxns,1);
        modelIrrev.c(objectiveColIrrev(i,1),1)=1;

        %find max objective
        FBAsolution = optimizeCbModel(modelIrrev);
        if (FBAsolution.stat ~= 1)
            not_solved=1;
            display('Failed to solve initial FBA problem');
            return
        end
        maxObjective(i)=FBAsolution.f;
    end

    model2gimme = modelIrrev;
    model2gimme.c = zeros(nIrrevRxns,1);


    for i=1:nIrrevRxns
        if (expressionColIrrev(i,1) > -1)   %if not absent reaction
            if (expressionColIrrev(i,1) < cutoff)
                model2gimme.c(i,1) = cutoff-expressionColIrrev(i,1);
            end
        end
    end

    for i=1:size(objectiveColIrrev,1)
        model2gimme.lb(objectiveColIrrev(i,1),1) = objectiveColIrrev(i,2) * maxObjective(i);
    end

    gimmeSolution = optimizeCbModel(model2gimme,'min');

    if (gimmeSolution.stat ~= 1)
    %%        gimme_not_solved=1;
    %        display('Failed to solve GIMME problem'); 
    %        return
    gimmeSolution.x = zeros(nIrrevRxns,1);
    end

    reactionActivityIrrev = zeros(nIrrevRxns,1);
    for i=1:nIrrevRxns
        if ((expressionColIrrev(i,1) > cutoff) | (expressionColIrrev(i,1) == -1))
            reactionActivityIrrev(i,1)=1;
        elseif (gimmeSolution.x(i,1) > 0)
            reactionActivityIrrev(i,1)=2;
        end
    end

    %Translate reactionActivity to reversible model
    reactionActivity = zeros(nRxns,1);
    for i=1:nRxns
        for j=1:size(rev2irrev{i,1},2)
            if (reactionActivityIrrev(rev2irrev{i,1}(1,j)) > reactionActivity(i,1))
                reactionActivity(i,1) = reactionActivityIrrev(rev2irrev{i,1}(1,j));
            end
        end
    end

    %Calculate flux distribution for original model
    fluxes = convertIrrevFluxDistribution(gimmeSolution.x, matchRev);
end

function [modelIrrev,matchRev,rev2irrev,irrev2rev] = convertToIrreversible(model)
%convertToIrreversible Convert model to irreversible format
%
% Copied from the cobra toolbox with some fixes.
%
% [modelIrrev,matchRev,rev2irrev,irrev2rev] = convertToIrreversible(model)
%
%INPUT
% model         COBRA model structure
%
%OUTPUTS
% modelIrrev    Model in irreversible format
% matchRev      Matching of forward and backward reactions of a reversible
%               reaction
% rev2irrev     Matching from reversible to irreversible reactions
% irrev2rev     Matching from irreversible to reversible reactions
%
% Uses the reversible list to construct a new model with reversible
% reactions separated into forward and backward reactions.  Separated
% reactions are appended with '_f' and '_b' and the reversible list tracks
% these changes with a '1' corresponding to separated forward reactions.
% Reactions entirely in the negative direction will be reversed and
% appended with '_r'.
%
% written by Gregory Hannum 7/9/05
%
% Modified by Markus Herrgard 7/25/05
% Modified by Jan Schellenberger 9/9/09 for speed.
%
% Some fixes by Daniel Machado, 2013.

    %declare variables
    modelIrrev.S = spalloc(size(model.S,1),0,2*nnz(model.S));
    modelIrrev.rxns = [];
    modelIrrev.rev = zeros(2*length(model.rxns),1);
    modelIrrev.lb = zeros(2*length(model.rxns),1);
    modelIrrev.ub = zeros(2*length(model.rxns),1);
    modelIrrev.c = zeros(2*length(model.rxns),1);
    matchRev = zeros(2*length(model.rxns),1);

    nRxns = size(model.S,2);
    irrev2rev = zeros(2*length(model.rxns),1);

    %loop through each column/rxn in the S matrix building the irreversible
    %model
    cnt = 0;
    for i = 1:nRxns
        cnt = cnt + 1;

        %expand the new model (same for both irrev & rev rxns  
        modelIrrev.rev(cnt) = model.rev(i);
        irrev2rev(cnt) = i;

        % Check if reaction is declared as irreversible, but bounds suggest
        % reversible (i.e., having both positive and negative bounds
        if (model.ub(i) > 0 && model.lb(i) < 0) && model.rev(i) == false
            model.rev(i) = true;
            warning(cat(2,'Reaction: ',model.rxns{i},' is classified as irreversible, but bounds are positive and negative!'))

        end


    % FIX Daniel M. 2013-01-11 - Temporary fix
    % This causes problems when comparing two models under different environmental
    % conditions, because they can end up with flux vectors of different sizes.  

    % Reaction entirely in the negative direction
    %    if (model.ub(i) <= 0 && model.lb(i) < 0)
    %         % Retain original bounds but reversed
    %         modelIrrev.ub(cnt) = -model.lb(i);
    %         modelIrrev.lb(cnt) = -model.ub(i);
    %         % Reverse sign
    %         modelIrrev.S(:,cnt) = -model.S(:,i);
    %         modelIrrev.c(cnt) = -model.c(i);
    %         modelIrrev.rxns{cnt} = [model.rxns{i} '_r'];
    %         model.rev(i) = false;
    %         modelIrrev.rev(cnt) = false;
    %     else
            % Keep positive upper bound
            modelIrrev.ub(cnt) = model.ub(i);
            %if the lb is less than zero, set the forward rxn lb to zero 

    %         if model.lb(i) < 0
    %            modelIrrev.lb(cnt) = 0;
    %         else
    %            modelIrrev.lb(cnt) = model.lb(i);
    %         end
            modelIrrev.lb(cnt) = max(0, model.lb(i));
            modelIrrev.ub(cnt) = max(0, model.ub(i));

            modelIrrev.S(:,cnt) = model.S(:,i);
            modelIrrev.c(cnt) = model.c(i);
            modelIrrev.rxns{cnt} = model.rxns{i};
    %    end


        %if the reaction is reversible, add a new rxn to the irrev model and
        %update the names of the reactions with '_f' and '_b'
        if model.rev(i) == true
            cnt = cnt + 1;
            matchRev(cnt) = cnt - 1;
            matchRev(cnt-1) = cnt;
            modelIrrev.rxns{cnt-1} = [model.rxns{i} '_f'];
            modelIrrev.S(:,cnt) = -model.S(:,i);
            modelIrrev.rxns{cnt} = [model.rxns{i} '_b'];
            modelIrrev.rev(cnt) = true;

    % FIX Daniel M. 2013-01-09 - if original reaction has a positive lb,
    % backwards reaction should have nonnegative upper bound.
    %        modelIrrev.lb(cnt) = 0;
    %        modelIrrev.ub(cnt) = -model.lb(i);    
            modelIrrev.lb(cnt) = max(0, -model.ub(i));
            modelIrrev.ub(cnt) = max(0, -model.lb(i));

            modelIrrev.c(cnt) = 0;
            rev2irrev{i} = [cnt-1 cnt];
            irrev2rev(cnt) = i;
        else
            matchRev(cnt) = 0;
            rev2irrev{i} = cnt;
        end
    end

    rev2irrev = columnVector(rev2irrev);
    irrev2rev = irrev2rev(1:cnt);
    irrev2rev = columnVector(irrev2rev);

    % Build final structure
    modelIrrev.S = modelIrrev.S(:,1:cnt);
    modelIrrev.ub = columnVector(modelIrrev.ub(1:cnt));
    modelIrrev.lb = columnVector(modelIrrev.lb(1:cnt));
    modelIrrev.c = columnVector(modelIrrev.c(1:cnt));
    modelIrrev.rev = modelIrrev.rev(1:cnt);
    modelIrrev.rev = columnVector(modelIrrev.rev == 1);
    modelIrrev.rxns = columnVector(modelIrrev.rxns); 
    modelIrrev.mets = model.mets;
    matchRev = columnVector(matchRev(1:cnt));
    modelIrrev.match = matchRev;
    if (isfield(model,'b'))
        modelIrrev.b = model.b;
    end
    if isfield(model,'description')
        modelIrrev.description = [model.description ' irreversible'];
    end
    if isfield(model,'subSystems')
        modelIrrev.subSystems = model.subSystems(irrev2rev);
    end
    if isfield(model,'genes')
        modelIrrev.genes = model.genes;
        genemtxtranspose = model.rxnGeneMat';
        modelIrrev.rxnGeneMat = genemtxtranspose(:,irrev2rev)';
        modelIrrev.rules = model.rules(irrev2rev);
        modelIrrev.grRules = model.grRules(irrev2rev);
    end
    modelIrrev.reversibleModel = false;

end