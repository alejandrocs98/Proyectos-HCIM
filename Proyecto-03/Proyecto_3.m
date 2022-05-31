initCobraToolbox
%% Load Models
model = readCbModel('../Codes-and-models/iMAC868_GSM.xlsx')
model_native = readCbModel('../Codes-and-models/iMAC868_GSM_native.xlsx')
model_methane = readCbModel('../Codes-and-models/iMAC868_GSM_methane.xlsx')
%% Principal Component Analysis
[U,S,V] = svds(model.S)
Vt = V.'
%% FBA for Biomass in native substrate and in substrate with Methane
FBA_overall_native = optimizeCbModel(model_native, 'max') 
FBA_overall_methane = optimizeCbModel(model_methane, 'max')
%% FBA for Acetate in native substrate and in substrate with Methane
model_native = changeObjective(model_native, 'EX_ac[e]')
model_methane = changeObjective(model_methane, 'EX_ac[e]')
FBA_acetate_native = optimizeCbModel(model_native, 'max') 
FBA_acetate_methane = optimizeCbModel(model_methane, 'max')
%% FBA for Methanol in native substrate and in substrate with Methane
model_native = changeObjective(model_native, 'EX_meoh[e]')
model_methane = changeObjective(model_methane, 'EX_meoh[e]')
FBA_methanol_native = optimizeCbModel(model_native, 'max') 
FBA_methanol_methane = optimizeCbModel(model_methane, 'max')
%% PhPP in native substrate: function of CO2 and CH4 = No Phenotypic Phase Plane results
% % For Biomass
% model_native = changeObjective(model_native, 'overall');
% [growthRates, shadowPrices1, shadowPrices2] = phenotypePhasePlane(model_native, 'EX_co2[e]', 'EX_ch4[e]'); 
% % For acetate
% model_native = changeObjective(model_native, 'EX_ac[e]'); 
% [growthRates, shadowPrices1, shadowPrices2] = phenotypePhasePlane(model_native, 'EX_co2[e]', 'EX_ch4[e]'); 
% % For methanol 
% model_native = changeObjective(model_native, 'EX_meoh[e]'); 
% [growthRates, shadowPrices1, shadowPrices2] = phenotypePhasePlane(model_native, 'EX_co2[e]', 'EX_ch4[e]'); 
%% PhPP for Biomass in substrate with Methane: function of CO2 and CH4
model_methane = changeObjective(model_methane, 'overall')
[growthRates_overall, shadowPrices1_overall, shadowPrices2_overall] = phenotypePhasePlane(model_methane, 'EX_co2[e]', 'EX_ch4[e]')
%% PhPP for Acetate in substrate with Methane: function of CO2 and CH4
model_methane = changeObjective(model_methane, 'EX_ac[e]')
[growthRates_ac, shadowPrices1_ac, shadowPrices2_ac] = phenotypePhasePlane(model_methane, 'EX_co2[e]', 'EX_ch4[e]')
%% PhPP for Methanol in substrate with Methane: function of CO2 and CH4
model_methane = changeObjective(model_methane, 'EX_meoh[e]')
[growthRates_meoh, shadowPrices1_meoh, shadowPrices2_meoh] = phenotypePhasePlane(model_methane, 'EX_co2[e]', 'EX_ch4[e]')
%% FVA for Biomass in native substrate and in substrate with methane
model_native = changeObjective(model_native, 'overall')
model_methane = changeObjective(model_methane, 'overall')
[minflux_ns_overall, maxflux_ns_overall] = fluxVariability(model_native)
[minflux_m_overall, maxflux_m_overall] = fluxVariability(model_methane)
%% FVA for Acetate in native substrate and in substrate with methane
model_native = changeObjective(model_native, 'EX_ac[e]')
model_methane = changeObjective(model_methane, 'EX_ac[e]')
[minflux_ns_ac, maxflux_ns_ac] = fluxVariability(model_native)
[minflux_m_ac, maxflux_m_ac] = fluxVariability(model_methane)
%% FVA for Methanol in native substrate and in substrate with methane
model_native = changeObjective(model_native, 'EX_meoh[e]')
model_methane = changeObjective(model_methane, 'EX_meoh[e]')
[minflux_ns_meoh, maxflux_ns_meoh] = fluxVariability(model_native)
[minflux_m_meoh, maxflux_m_meoh] = fluxVariability(model_methane)