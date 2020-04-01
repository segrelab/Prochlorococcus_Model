%base model used: iJC568a.mat - look in model log file

%1
%adding R00188 to resolve fake Adenosine 3,5 bisphophate transport
model=addReaction(model,'R00188','metaboliteList',{'Adenosine_3_5_bisphosphate[c]',...
    'H2O[c]','AMP[c]','Orthophosphate[c]'},...
    'stoichCoeffList',[-1;-1;1;1], 'reversible', true);
%model saved as iJC568a in the model_files folder

%2
%resolving 7 aminomethyl 7 carbaguanine. Adding three reactions sinking 
%it into tRNA
model=addReaction(model,'R10209m','metaboliteList',{'7_Aminomethyl_7_carbaguanine[c]',...
    'preqtRNA[c]',},...
    'stoichCoeffList',[-1;1], 'reversible', false);
model=addReaction(model,'SAMTRIm','metaboliteList',{'S_Adenosyl_L_methionine[c]',...
    'preqtRNA[c]','H[c]','Adenine[c]','L_Methionine[c]','epxytRNA[c]'},...
    'stoichCoeffList',[-1;-1;1;1;1;1], 'reversible', false);
model=addReaction(model,'EPXQRm','metaboliteList',{'H[c]','epxytRNA[c]','NADPH[c]','H2O[c]',...
    'NADP[c]','tRNA[c]'},...
    'stoichCoeffList',[-1;-1;-1;1;1;1], 'reversible', false); 
%since 7 aminomethyl 7 carbaguanine still prefes to go out of the fake
%exchange if possible even with the new reactions in place I am removing
%the fake exchange.
model=removeRxns(model,'7NMeth7carbEX');

%resolving nicotinate, adding R01724
model=addReaction(model,'R01724','metaboliteList',{'Nicotinate_D_ribonucleotide[c]','Diphosphate[c]',...
    'ADP[c]','Orthophosphate[c]','Nicotinate[c]','5_Phospho_alpha_D_ribose_1_diphosphate[c]',...
    'ATP[c]','H2O[c]','H[c]'},'stoichCoeffList',[-1;-1;-1;-1;1;1;1;1;1], 'reversible', true);

%resolving S-Adenosyl-4 methylthio-2-oxobutanoate, unresolvable- only one
%reaction known to produce this metabolite and no known reactions
%degrading it. Keep for next round of curation.

%Resolving S-Methyl-5-thio-D-ribose,hopefully. I found some upstream
%missing reactions that might help. I'm updateing them into the model.
%Directionality was check with the SEED DB. - saved as iJC568a.4.mat

%Adding R01777
model=addReaction(model,'R01777','metaboliteList',{'Succinyl_CoA[c]','L_Homoserine[c]',...
    'CoA[c]','O_Succinyl_L_homoserine[c]'},'stoichCoeffList',[-1;-1;1;1], 'reversible', false);
%Adding R01776
model=addReaction(model,'R01776','metaboliteList',{'Acetyl[c]','L_Homoserine[c]',...
    'CoA[c]','O_Acetyl_L_homoserine[c]'},'stoichCoeffList',[-1;-1;1;1], 'reversible', false);
%Adding R03260
model=addReaction(model,'R03260','metaboliteList',{'O_Succinyl_L_homoserine[c]','L_Cysteine[c]',...
    'L_Cystathionine[c]','Succinate[c]'},'stoichCoeffList',[-1;-1;1;1], 'reversible', true);
%adding R08549
model=addReaction(model,'R08549','metaboliteList',{'2_Oxoglutarate[c]','CoA[c]',...
    'NAD[c]','Succinyl_CoA[c]','CO2[c]','NADH[c]','H[c]'},'stoichCoeffList',...
    [-1;-1;-1;1;1;1;1], 'reversible', false);

% Adding citrate, pyruvate and acetate exchanges as found in Biller et.al,
% 2018. Also addint formate exchange as found in Bertillson et.al, 2015

% adding acetate
model=addReaction(model,'Trans_Acetate','metaboliteList',{'Acetate[e]','Acetate[c]'},...
    'stoichCoeffList',[-1;1], 'reversible', true);
model = addReaction(model, 'AcetateEX', 'metaboliteList', {'Acetate[e]'} ,'stoichCoeffList', [-1]);

% adding citrate
model=addReaction(model,'Trans_Citrate','metaboliteList',{'Citrate[e]','Citrate[c]'},...
    'stoichCoeffList',[-1;1], 'reversible', true);
model = addReaction(model, 'CitrateEX', 'metaboliteList', {'Citrate[e]'} ,'stoichCoeffList', [-1]);

% adding formate
model=addReaction(model,'Trans_Formate','metaboliteList',{'Formate[e]','Formate[c]'},...
    'stoichCoeffList',[-1;1], 'reversible', true);
model = addReaction(model, 'FormateEX', 'metaboliteList', {'Formate[e]'} ,'stoichCoeffList', [-1]);

% adding pyruvate
model=addReaction(model,'Trans_Pyruvate','metaboliteList',{'Pyruvate[e]','Pyruvate[c]'},...
    'stoichCoeffList',[-1;1], 'reversible', true);
model = addReaction(model, 'PyruvateEX', 'metaboliteList', {'Pyruvate[e]'} ,'stoichCoeffList', [-1]);


% Adding exchanges after looking at metabolights study MTBLS567 about proc.
% exometabolome. These data were not collected for MED4 but for another
% proc strain. Still valueble. Found 13 exchanges to add to metabolites
% already in the model, 1 reaction plus exchange (only if connected). Also
% included 2 exchanges found interesting from transportDB.
% this will be saved as iJC568a.6

% adding fumarate
model=addReaction(model,'Trans_Fumarate','metaboliteList',{'Fumarate[e]','Fumarate[c]'},...
    'stoichCoeffList',[-1;1], 'reversible', true);
model = addReaction(model, 'FumarateEX', 'metaboliteList', {'Fumarate[e]'} ,'stoichCoeffList', [-1]);

%adding glutathione
model=addReaction(model,'Trans_Glutathione','metaboliteList',{'Glutathione[e]','Glutathione[c]'},...
    'stoichCoeffList',[-1;1], 'reversible', true);
model = addReaction(model, 'GlutathioneEX', 'metaboliteList', {'Glutathione[e]'} ,'stoichCoeffList', [-1]);

% adding 4_Methyl_2_oxopentanoate

model=addReaction(model,'Trans_4_Methyl_2_oxopentanoate','metaboliteList',{'4_Methyl_2_oxopentanoate[e]','4_Methyl_2_oxopentanoate[c]'},...
    'stoichCoeffList',[-1;1], 'reversible', true);
model = addReaction(model, '4_Methyl_2_oxopentanoateEX', 'metaboliteList', {'4_Methyl_2_oxopentanoate[e]'} ,'stoichCoeffList', [-1]);

% adding Pantothenate
model=addReaction(model,'Trans_Pantothenate','metaboliteList',{'Pantothenate[e]','Pantothenate[c]'},...
    'stoichCoeffList',[-1;1], 'reversible', true);
model = addReaction(model, 'PantothenateEX', 'metaboliteList', {'Pantothenate[e]'} ,'stoichCoeffList', [-1]);

%adding succinate
model=addReaction(model,'Trans_Succinate','metaboliteList',{'Succinate[e]','Succinate[c]'},...
    'stoichCoeffList',[-1;1], 'reversible', true);
model = addReaction(model, 'SuccinateEX', 'metaboliteList', {'Succinate[e]'} ,'stoichCoeffList', [-1]);

% adding thymidine
model=addReaction(model,'Trans_Thymidine','metaboliteList',{'Thymidine[e]','Thymidine[c]'},...
    'stoichCoeffList',[-1;1], 'reversible', true);
model = addReaction(model, 'ThymidineEX', 'metaboliteList', {'Thymidine[e]'} ,'stoichCoeffList', [-1]);

% adding Xanthosine
model=addReaction(model,'Trans_Xanthosine','metaboliteList',{'Xanthosine[e]','Xanthosine[c]'},...
    'stoichCoeffList',[-1;1], 'reversible', true);
model = addReaction(model, 'XanthosineEX', 'metaboliteList', {'Xanthosine[e]'} ,'stoichCoeffList', [-1]);

% adding 4-aminobenzoate
model=addReaction(model,'Trans_4_Aminobenzoate','metaboliteList',{'4_Aminobenzoate[e]','4_Aminobenzoate[c]'},...
    'stoichCoeffList',[-1;1], 'reversible', true);
model = addReaction(model, '4_AminobenzoateEX', 'metaboliteList', {'4_Aminobenzoate[e]'} ,'stoichCoeffList', [-1]);

% adding 5-methylthioadenosine
model=addReaction(model,'Trans_5_Methylthioadenosine','metaboliteList',{'5_Methylthioadenosine[e]','5_Methylthioadenosine[c]'},...
    'stoichCoeffList',[-1;1], 'reversible', true);
model = addReaction(model, '5_MethylthioadenosineEX', 'metaboliteList', {'5_Methylthioadenosine[e]'} ,'stoichCoeffList', [-1]);

% adding adenine
model=addReaction(model,'Trans_Adenine','metaboliteList',{'Adenine[e]','Adenine[c]'},...
    'stoichCoeffList',[-1;1], 'reversible', true);
model = addReaction(model, 'AdenineEX', 'metaboliteList', {'Adenine[e]'} ,'stoichCoeffList', [-1]);

% adding dethiobiotin
model=addReaction(model,'Trans_Dethiobiotin','metaboliteList',{'Dethiobiotin[e]','Dethiobiotin[c]'},...
    'stoichCoeffList',[-1;1], 'reversible', true);
model = addReaction(model, 'DethiobiotinEX', 'metaboliteList', {'Dethiobiotin[e]'} ,'stoichCoeffList', [-1]);

% adding guanine
model=addReaction(model,'Trans_Guanine','metaboliteList',{'Guanine[e]','Guanine[c]'},...
    'stoichCoeffList',[-1;1], 'reversible', true);
model = addReaction(model, 'GuanineEX', 'metaboliteList', {'Guanine[e]'} ,'stoichCoeffList', [-1]);

% adding guanosine
model=addReaction(model,'Trans_Guanosine','metaboliteList',{'Guanosine[e]','Guanosine[c]'},...
    'stoichCoeffList',[-1;1], 'reversible', true);
model = addReaction(model, 'GuanosineEX', 'metaboliteList', {'Guanosine[e]'} ,'stoichCoeffList', [-1]);

% adding putrescine
model=addReaction(model,'Trans_Putrescine','metaboliteList',{'Putrescine[e]','Putrescine[c]'},...
    'stoichCoeffList',[-1;1], 'reversible', true);
model = addReaction(model, 'PutrescineEX', 'metaboliteList', {'Putrescine[e]'} ,'stoichCoeffList', [-1]);

%adding spermidine
model=addReaction(model,'Trans_Spermidine','metaboliteList',{'Spermidine[e]','Spermidine[c]'},...
    'stoichCoeffList',[-1;1], 'reversible', true);
model = addReaction(model, 'SpermidineEX', 'metaboliteList', {'Spermidine[e]'} ,'stoichCoeffList', [-1]);

% add exchange for 4-hydroxybenzoate
model=addReaction(model,'Trans_4_Hydroxybenzoate','metaboliteList',{'4_Hydroxybenzoate[e]','4_Hydroxybenzoate[c]'},...
    'stoichCoeffList',[-1;1], 'reversible', true);
model = addReaction(model, '4_HydroxybenzoateEX', 'metaboliteList', {'4_Hydroxybenzoate[e]'} ,'stoichCoeffList', [-1]);

% add ethanol exchange, ethanol can passively diffuse out of a phospholipid
% layer so if it's produced it can go out
model=addReaction(model,'Trans_ethanol','metaboliteList',{'Ethanol[e]','Ethanol[c]'},...
    'stoichCoeffList',[-1;1], 'reversible', true);
model = addReaction(model, 'EthanolEX', 'metaboliteList', {'Ethanol[e]'} ,'stoichCoeffList', [-1]);

% Add S-Malate exchange and transport following the Metabolic evolution and the self-organization of ecosystems by Baakman et al., (2017)
model=addReaction(model,'Trans_S_Malate','metaboliteList',{'S_Malate[e]','S_Malate[c]'},...
    'stoichCoeffList',[-1;1], 'reversible', true);
model = addReaction(model, 'S_MalateEX', 'metaboliteList', {'S_Malate[e]'} ,'stoichCoeffList', [-1]);
%% Changes added to go from iSO595 v4 to v5
% Change reaction R06050 from reversible to forward. This reaction convers
% Amylose to glucose 1-phosphate. Knockout mutant of glgC (PMM0769 -> R00948) does not store glycogen
% (Shinde et al., 2020), and the change introduced here matches this phenotype. 
model = changeRxnBounds(model, 'R06050', 0, 'l');

% Add transport of hydrogen peroxide. Photosynthesis is known to create
% reactive oxygen species and Prochlorococcus does not have catalase. It is
% assumed (but not completely clear) that hydrogen peroxide can diffuse
% across the cell membrane
model = addMetabolite(model, 'Hydrogen_peroxide[e]', 'metName', 'Hydrogen peroxide', 'metFormula', 'H2O2','Charge', 0);
model = addReaction(model, 'Trans_H2O2', 'reactionName', 'H2O2 diffusion', 'metaboliteList', {'Hydrogen_peroxide[e]', 'Hydrogen_peroxide[c]'},'stoichCoeffList', [1, -1], 'lowerBound', -1000);
model = addReaction(model, 'H2O2EX', 'reactionName', 'H2O2 exchange','metaboliteList', {'Hydrogen_peroxide[e]'}, 'stoichCoeffList', [-1]);

% Add missing reaction 6PG-dehydratase/ 4.2.1.12 in Entner-Doudoroff Pathway. 
% Chen et al., 2016 clearly demonstrate that ED is present in Synechocystis, 
% suggesting that the gene slr0452 is encoding for this enzyme. A blast search 
% of this gene shows that PMM0774 is a homolog of this gene in P. marinus
% MED4. A search on equilibrator show a dG'm of -43 kJ/mol, so therefore I
% set this reaction to forward. It is striking that adding this reaction
% increase the growth rate (only limited by rubisco) from 0.0786 to 0.0997
model = addReaction(model, 'R02036', 'reactionName', '6-Phospho-D-gluconate hydro-lyase','metaboliteList', {'6_Phospho_D_gluconate[c]', '2_Dehydro_3_deoxy_6_phospho_D_gluconate[c]', 'H2O[c]'},...
                    'stoichCoeffList', [-1, 1,1], 'lowerBound', 0, 'geneRule', 'PMM0774');
                
%% Changes added to go to iSO595 v6
% By setting these two reactions to the same direction (producing
% Chlorophyllide) we avoid that these reactions are used in a cycle to generate ATP /
% NADPH
model = changeRxnBounds(model, 'R06282', -1000, 'l');
model = changeRxnBounds(model, 'R06282', 0, 'u');
model = changeRxnBounds(model, 'R03845', -1000, 'l');
model = changeRxnBounds(model, 'R03845', 0, 'u');

%%

writeCbModel(model, 'format', 'sbml', 'fileName', 'iSO595v6.xml');
writeCbModel(model, 'format', 'mat', 'fileName', 'iSO595v6.mat');
