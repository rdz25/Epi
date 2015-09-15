% function [comb, ep] = epiComb(positive, indexposi, episolo, model)
    %
    % epistasis= e of rxn1&2 - e of rxn1 *e of rxn2
    %
    % INPUTS
    %       positive - array of rxn that produce positive
    %       indexposi - index of rxn in the model
    %       episolo - vector of individual positive effects values when intake is incremented by -1 
    %       model - model used for FBA
    % OUTPUTS
    %       comb - vector of name of pair with positive epistasis value
    %       ep - vector of epistasis values for each combination
    %
    % I know it's kinda messy, maybe I should make an object/class
    % Richard Zhang
    
    count=0;
    ep=[];
    comb=[];
    load('ecoli_core_model.mat');
    for i=1:length(positive)
        loweri=model.lb(indexposi(i)); %positive index lb
        for j=i+1:length(positive)
            lowerj=model.lb(indexposi(j));
            model = changeRxnBounds(model,positive{i},loweri-1,'l');    %model = changeRxnBounds(model,'EX_ac(e)',-1,'l');
            model = changeRxnBounds(model,positive{j},lowerj-1,'l');
            FBAsolution = optimizeCbModel(model,'max');
            if episolo(i)*episolo(j)<FBAsolution.f  %positive epistasis
                count=count+1;
                ep(count)=FBAsolution.f-episolo(i)*episolo(j);
                comb{count}=[positive{i},'&',positive{j}];
            end
            model = changeRxnBounds(model,positive{i},loweri,'l');  %return lower bound to normal
            model = changeRxnBounds(model,positive{j},lowerj,'l');
        end
    end
% end