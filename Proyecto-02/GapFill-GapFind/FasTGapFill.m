% nitialize the cobra toolbox
initCobraToolbox(false);
changeCobraSolver ('gurobi', 'all', 1);
% Load model
pathModel = '/home/acs98/Documents/Semestre-XII/Herramientas-computacionales-para-ingenieria-metabolica/Matlab/cobratoolbox/'
modelFileName = 'iAF1260.mat'
load([pathModel, modelFileName])
model = iAF1260;
% model = readCbModel(modelFileName)

% Detect deadend metabolites
outputMets = detectDeadEnds(model)
% print the corresponding metabolite names
DeadEnds = model.mets(outputMets)
% Identify associated reactions
[rxnList, rxnFormulaList] = findRxnsFromMets(model, DeadEnds)
% look at the lower and upper bounds of these reactions
model.lb(find(ismember(model.rxns,rxnList)))
model.ub(find(ismember(model.rxns,rxnList)))
% GapFind
[allGaps, rootGaps, downstreamGaps] = gapFind(model, 'true');
% Identification of blocked reactions
BlockedReactions = findBlockedReaction(model)
% FastGapFill allows to set different priorities for reaction types (MetabolicRxns = internal metabolic reactions in the reference database (here, KEGG), ExchangeRxns = exchange reactions, and TransportRxns = intra- and extracellular transport reactions) using the weights. The lower the weight for a reaction type, the higher is its priority. Generally, a metabolic reaction should be prioritised in a solution over transport and exchange reactions.
weights.MetabolicRxns = 0.1; % Kegg metabolic reactions
weights.ExchangeRxns = 0.5; % Exchange reactions
weights.TransportRxns = 10; % Transport reactions
% Prepare the output table with statistics
cnt = 1;
Stats{cnt,1} = 'Model name';cnt = cnt+1;
Stats{cnt,1} = 'Size S (original model)';cnt = cnt+1;
Stats{cnt,1} = 'Number of compartments';cnt = cnt+1;
Stats{cnt,1} = 'List of compartments';cnt = cnt+1;
Stats{cnt,1} = 'Number of blocked reactions';cnt = cnt+1;
Stats{cnt,1} = 'Number of solvable blocked reactions';cnt = cnt+1;
Stats{cnt,1} = 'Size S (flux consistent)';cnt = cnt+1;
Stats{cnt,1} = 'Size SUX (including solvable blocked reactions)';cnt = cnt+1;
Stats{cnt,1} = 'Number of added reactions (all)';cnt = cnt+1;
Stats{cnt,1} = 'Number of added metabolic reactions ';cnt = cnt+1;
Stats{cnt,1} = 'Number of added transport reactions ';cnt = cnt+1;
Stats{cnt,1} = 'Number of added exchange reactions ';cnt = cnt+1;
Stats{cnt,1} = 'Time preprocessing';cnt = cnt+1;
Stats{cnt,1} = 'Time fastGapFill';cnt = cnt+1;
% Initiate parameters
col = 1;
RxnList={};
cnt = 1;
% Remove constraints from exchange reactions
EX = strmatch('EX_',model.rxns);
model.lb(EX)=-100;
model.ub(EX)=100;
clear EX
% Get basic model statistics
Stats{cnt,i+1} = modelFileName;cnt = cnt+1;
[a,b] = size(model.S);
Stats{cnt,i+1} = strcat(num2str(a),'x',num2str(b));cnt = cnt+1;
% List of compartments in the model that will be considered during the gap filling.
[tok,rem] = strtok(model.mets,'\[');
rem = unique(rem);
Stats{cnt,i+1} = num2str(length(rem));cnt = cnt+1;
Rem = rem{1};
for j = 2:length(rem)
Rem = strcat(Rem,',',rem{j});
end
Stats{cnt,i+1} = Rem;cnt = cnt+1;
clear Rem tok rem;
% Prepare fastGapFill
tic; [consistModel,consistMatricesSUX,BlockedRxns] = prepareFastGapFill(model);
tpre=toc;
% Add on more statistics to the table. Here, we will report how many gaps cannot be filled at all in the starting reconstruction.
Stats{cnt,i+1} = num2str(length(BlockedRxns.allRxns));cnt = cnt+1;
Stats{cnt,i+1} = num2str(length(BlockedRxns.solvableRxns));cnt = cnt+1;
[a,b] = size(consistModel.S);
Stats{cnt,i+1} = strcat(num2str(a),'x',num2str(b));cnt = cnt+1;
[a,b] = size(consistMatricesSUX.S);
Stats{cnt,i+1} = strcat(num2str(a),'x',num2str(b));cnt = cnt+1;
% Perform fastGapsFill
epsilon = 1e-4;
tic; [AddedRxns] = fastGapFill(consistMatricesSUX,epsilon, weights);
tgap=toc;
Stats{cnt,i+1} = num2str(length(AddedRxns.rxns));cnt = cnt+1;
% Postprocessing of results
IdentifyPW = 0;
[AddedRxnsExtended] = postProcessGapFillSolutions(AddedRxns,model,BlockedRxns);
clear AddedRxns;
% Add further statistics about computed solutions to the Stats table.
Stats{cnt,i+1} = num2str(AddedRxnsExtended.Stats.metabolicSol);cnt = cnt+1;
Stats{cnt,i+1} = num2str(AddedRxnsExtended.Stats.transportSol);cnt = cnt+1;
Stats{cnt,i+1} = num2str(AddedRxnsExtended.Stats.exchangeSol);cnt = cnt+1;
Stats{cnt,i+1} = num2str(tpre);cnt = cnt+1;
Stats{cnt,i+1} = num2str(tgap);cnt = cnt+1;
clear a b
% Generate a reaction list
RxnList{1,col}=modelFileName;RxnList(2:length(AddedRxnsExtended.rxns)+1,col) = AddedRxnsExtended.rxns; col = col + 1;
RxnList{1,col}=modelFileName;RxnList(2:length(AddedRxnsExtended.rxns)+1,col) = AddedRxnsExtended.rxnFormula; col = col + 1;
RxnList{1,col}=modelFileName;RxnList(2:length(AddedRxnsExtended.rxns)+1,col) = AddedRxnsExtended.subSystem; col = col + 1;
clear AddedRxnsExtended tgap tpre j BlockedRxns i cnt consistM*
% Analysis of results
RxnList
% Identify the deadend metabolites
metid = detectDeadEnds(model,false);
% Print the metabolite abbreviations and metabolite names on the screen.
horzcat(model.mets(metid),model.metNames(metid))