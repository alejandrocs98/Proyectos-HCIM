initCobraToolbox
%% Principal Component Analysis
model = readCbModel('../Codes-and-models/iMAC868_GSM.xlsx');
[U,S,V] = svds(model.S);
Vt = V.';
%% FBA for overall growth in native substrate and in substrate with methane
model_native = readCbModel('../Codes-and-models/iMAC868_GSM_native.xlsx');
model_methane = readCbModel('../Codes-and-models/iMAC868_GSM_methane.xlsx');
FBA_native = optimizeCbModel(model_native, 'max') ;
FBA_methane = optimizeCbModel(model_methane, 'max');
%% PhPP in substarte with methane
[growthRates_m, shadowPrices1_m, shadowPrices2_m] = phenotypePhasePlane(model_methane, 'EX_ac[e]', 'EX_ch4[e]');
% [growthRates_m, shadowPrices1_m, shadowPrices2_m] = phenotypePhasePlane(model_methane, 'EX_ac[e]', 'EX_co2[e]');
% [growthRates_m, shadowPrices1_m, shadowPrices2_m] = phenotypePhasePlane(model_methane, 'EX_ch4[e]', 'EX_co2[e]');
% [growthRates_m, shadowPrices1_m, shadowPrices2_m] = phenotypePhasePlane(model_methane, 'EX_ch4[e]', 'EX_ch4[e]');
%% PhPP in native substrate
% [growthRates_ns, shadowPrices1_ns, shadowPrices2_ns] = phenotypePhasePlane(model_native, 'EX_ac[e]', 'EX_ch4[e]');
[growthRates_ns, shadowPrices1_ns, shadowPrices2_ns] = phenotypePhasePlane(model_native, 'EX_ac[e]', 'EX_co2[e]');
% [growthRates_ns, shadowPrices1_ns, shadowPrices2_ns] = phenotypePhasePlane(model_native, 'EX_ch4[e]', 'EX_co2[e]');
% [growthRates_ns, shadowPrices1_ns, shadowPrices2_ns] = phenotypePhasePlane(model_native, 'EX_ch4[e]', 'EX_ch4[e]');
%% Thermodynamic FBA
model_methane_ = readCbModel('../Codes-and-models/iMAC868_GSM_methane.xlsx');
model_methane_ = assignQualDir(model_methane_);
model_methane_ = configureSetupThermoModelInputs(model_methane_);
model_methane_ = setupThermoModel(model_methane_, 0.95)

% model_methane = estimateDG_temp(model_methane)
% model_methane = addThermoToModel(model_methane)