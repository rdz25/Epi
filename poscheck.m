rxnnames=[{'EX_ac(e)';'EX_acald(e)';'EX_akg(e)';'EX_co2(e)';'EX_etoh(e)';'EX_for(e)';'EX_fru(e)';'EX_fum(e)';'EX_glc(e)';'EX_gln_L(e)';'EX_glu_L(e)';'EX_h(e)';'EX_h2o(e)';'EX_lac_D(e)';'EX_mal_L(e)';'EX_nh4(e)';'EX_o2(e)';'EX_pi(e)';'EX_pyr(e)';'EX_succ(e)'}];
initCobraToolbox;
load('ecoli_core_model.mat');

model = changeObjective(model,'Biomass_Ecoli_core_w_GAM');
FBAsolution = optimizeCbModel(model,'max');

standard=FBAsolution.f;
positive=[];
indexposi=[];
episolo=[];
count=0;
for i=1:length(rxnnames)
    %upper=model.ub(i+19);  %effect of greater output not yet checked for
    lower=model.lb(i+19);   %currently checks for positive effect from greater intake/lower bound
    model = changeRxnBounds(model,rxnnames{i},lower-1,'l');   %always decrease lb? or also check pos effects for increased ub?
    FBAsolution = optimizeCbModel(model,'max');
    if FBAsolution.f>standard
        count= count+1;
        positive{count}= rxnnames{i};   %rxns that have positive effect on biomass
        indexposi(count)=i+19;
        episolo(count)=FBAsolution.f-standard;  %effect of solo mutation
    end
    model = changeRxnBounds(model,rxnnames{i},lower,'l');   %return value to normal
end