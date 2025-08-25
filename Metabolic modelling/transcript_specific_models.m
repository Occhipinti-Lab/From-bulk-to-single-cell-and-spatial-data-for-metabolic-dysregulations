%%% Code for setting the model's bounds according to the provided transcriptomic data
% Boundary provided by NE
% gamma = [2, 1.5, 1, 0.5];                            
% threshold = [0,10,25,50];                        
% lb = [-50,-30,-10,-5];   
% type = ["to_one","to_zero","binary"];
%clear, clc
%feasTol =1e-02;

load('human1.mat');
model = human1;


model = changeRxnBounds(model,'CITRtm',0.00368209748124444,'u');%CITRtm
model = changeRxnBounds(model,'CITRtm',0.00351384832876872,'l');
 
model = changeRxnBounds(model,'MAL_Lte',0.000843099,'u');%malate
model = changeRxnBounds(model,'MAL_Lte',0.000732398763287308,'l');

model = changeRxnBounds(model,'HMR_9048',-2,'1'); %o2
model = changeRxnBounds(model,'HMR_9048',0,'u');


model = changeRxnBounds(model,'HMR_9034',-0.436212496491953,'1'); %glucose
model = changeRxnBounds(model,'HMR_9034',-0.433578512991697,'u')

model = changeRxnBounds(model,'HMR_9135',0.602986807828337,'1');%lactate 
model = changeRxnBounds(model,'HMR_9135',0.602986807828337,'u')


model = changeRxnBounds(model,'HMR_9063',-0.110504785,'1'); %glutamine
model = changeRxnBounds(model,'HMR_9063',-0.105713826230368,'u');


model = changeRxnBounds(model,'HMR_9071',0.017333082,'1'); %glutamate
model = changeRxnBounds(model,'HMR_9071',0.019031666,'u');

model = changeRxnBounds(model,'HMR_9070',0.000750958,'l'); %aspartate
model = changeRxnBounds(model,'HMR_9070',0.004253695,'u');


model = changeRxnBounds(model,'HMR_9062',-0.005030074,'l'); %asparagine
model = changeRxnBounds(model,'HMR_9062',0.000125257,'u');


model = changeRxnBounds(model,'HMR_9068',0.00204167687175758,'1'); %proline
model = changeRxnBounds(model,'HMR_9068',0.0032211735344259,'u');


model = changeRxnBounds(model,'HMR_9066',-0.005552074,'1'); %arginine
model = changeRxnBounds(model,'HMR_9066',0.00264078139633385,'u');


model = changeRxnBounds(model,'HMR_9061',0.0125778798412317,'1'); %alanine
model = changeRxnBounds(model,'HMR_9061',0.0155448644221056,'u');



model = changeRxnBounds(model,'HMR_9069',-0.0235542728719781,'l'); %serine
model = changeRxnBounds(model,'HMR_9069',-0.0216456965485817,'u');


model = changeRxnBounds(model,'HMR_9067',0.000822982323203474,'l'); %glycine
model = changeRxnBounds(model,'HMR_9067',0.00191853,'u');


model = changeRxnBounds(model,'HMR_9041',-0.0103457563037912,'l'); %lysine
model = changeRxnBounds(model,'HMR_9041',-0.00994377466392736,'u');



model = changeRxnBounds(model,'HMR_9045',-0.001321123,'l'); %tryptophan
model = changeRxnBounds(model,'HMR_9045',-0.00100677090358275,'u');


model = changeRxnBounds(model,'HMR_9040',-0.015474964,'l'); %leucine
model = changeRxnBounds(model,'HMR_9040',-0.007186224,'u');


model = changeRxnBounds(model,'HMR_9064',-0.0081131265759043,'l'); %tyrosine
model = changeRxnBounds(model,'HMR_9064',-0.004467569,'u');


model = changeRxnBounds(model,'HMR_9043',-0.004821991,'l'); %phenylalanine
model = changeRxnBounds(model,'HMR_9043',-0.003758165,'u');


model = changeRxnBounds(model,'HMR_9039',-0.011731877,'l'); %isoleucine
model = changeRxnBounds(model,'HMR_9039',-0.00645537,'u');

model = changeRxnBounds(model,'HMR_9046',-0.009963059,'l'); %Valine
model = changeRxnBounds(model,'HMR_9046',-0.00728358,'u');

model = changeRxnBounds(model,'HMR_9044',-0.009634843,'l'); %threonine
model = changeRxnBounds(model,'HMR_9044',-0.00731823,'u');

model = changeRxnBounds(model,'HMR_9087',-0.001561292,'l'); %ornithine
model = changeRxnBounds(model,'HMR_9087',0.000728721893945584,'u');



% model.c(13350) = 0; %  remove biomass from being the objective if you want to optimise another function

% define parameters over which to iterate

gamma = [2];
threshold = [25];

% addpath(genpath('PATH_TO_COBRATOOLBOX'));     % add here the path to your cobratoolbox folder
% initCobraToolbox
% by using the human1 model(in our case)calculate the 
[reaction_expression, pos_genes_in_react_expr, ixs_geni_sorted_by_length] = compute_reaction_expression(model);

% load gene expression and gene ids
load('gene_exp.mat');             % patients on the rows, gene expression along the columns (it is a matrix of only numbers, no ids)
load('gene_ids.mat');                   % gene ids of the genes present in the gene expression dataset, sorted by the same order in the gene_exp matrix, it is a string array
load('patient_ids.mat');              % patient id, sorted by the same order in the gene_exp matrix
%gene_exp=load('gene-exp.mat');
genes = model.genes;
%your data gene
genes_in_dataset = gene_ids; 
C{size(gene_exp, 1), length(model.rxns)} = [];

GeneExpressionArray = ones(numel(genes),1); 

k = 1;
tic
for g = 1:numel(gamma)                
    changeCobraSolver('gurobi', 'all');
	gam = gamma(g);
    new_k = main(threshold, gam, gene_exp,genes_in_dataset,patient_ids,model,genes,GeneExpressionArray,g,reaction_expression,pos_genes_in_react_expr,ixs_geni_sorted_by_length,k);
    k = new_k;
end
toc

function new_k = main(threshold, gam, gene_exp,genes_in_dataset,patient_ids,model,genes,GeneExpressionArray,g,reaction_expression,pos_genes_in_react_expr,ixs_geni_sorted_by_length,k)
    gamma = gam;
    for tr = 1:numel(threshold)
        prc_cut = threshold(tr);                                                            
        % cut data with respect to a threshold  
         data = gene_exp;
        %data = gene_exp./mean(gene_exp);                  % get the fold change (If A is a matrix, then mean(A) returns a row vector containing the mean of each column.)
        % data = data.*(data>prctile(data,prc_cut,1));    % if you want to binaruse data according to a percentile                                        
                
        % applying the bounds
        fprintf("Iteration (k): %d, Gamma: %d, Threshold: %d\n",k,gamma(g),threshold(tr));
        for t=1:size(data,1)          % in here we select a unique profile(one patient at the time)
           	expr_profile = data(t,:);
            pos_genes_in_dataset = zeros(numel(genes),1);% gene in the model human 1
            for i=1:length(pos_genes_in_dataset)
                position = find(strcmp(genes{i},genes_in_dataset),1); 
                if ~isempty(position)                                   
                    pos_genes_in_dataset(i) = position(1);              
                    GeneExpressionArray(i) = expr_profile(pos_genes_in_dataset(i));         
                end
            end
            if or(sum(isinf(GeneExpressionArray)) >= 1, sum(isnan(GeneExpressionArray)) >= 1)
                fprintf("\nError in the gene expression data!");
            end
            [fluxes] = transcriptomic_bounds(gamma(g), GeneExpressionArray, model, genes, reaction_expression, pos_genes_in_react_expr, ixs_geni_sorted_by_length);
            fprintf("\nNon-zero reactions: %d, sample number: %d, Solution status: %d\n", length(find(fluxes.v)), t, fluxes.stat);
            %fprintf("\nNon-zero reactions: %d, sample number: %d, Solution status: %d\n", length(find(fluxes.v)), t, fluxes.stat);
            C(k,:) = num2cell([patient_ids(k), transpose(fluxes.v)]);
            %C(k,:) = num2cell([patient_ids(k), transpose(fluxes)]);
            k = k +1;
        end
    end
    new_k = k; 
    T = cell2table(C);
    T.Properties.VariableNames{1} = 'patient_id';  
    writetable(T,"D:\Second-year\Code\FBA\flux-test.csv");     % set the folder where you want the data to be saved
end


%Find the index of a reaction, e.g., ‘ATPM', in the model
%rxnList = 'ATPM';
%rxnID = findRxnIDs(modelEcore, rxnList)
%Set ATP demand at the real rate now
%valATPM = 1.07;
%modelCore = changeRxnBounds(modelCore,'DM_atp_c_',valATPM,'l');
%modelCore = changeRxnBounds(modelCore,'DM_atp_c_',1000,'u');
%
%
