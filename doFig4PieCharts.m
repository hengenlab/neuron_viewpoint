%make timeline plot
clear

% load bibliographic spreadsheet
opts = spreadsheetImportOptions("NumVariables", 48);

% Specify sheet and range
opts.Sheet = "Sheet1";
opts.DataRange = "A2:AV320";

% Specify column names and types
opts.VariableNames = ["mouse", "rat", "cat", "monkey", "human", "fish", "turtle", "leech", "crayfish", "culture", "retina", "spikes", "CaImaging", "LFP", "voltageImaging", "wholecellPatch", "EEG", "MEG", "BOLD", "avalanches", "LRTC", "IsingFit", "ARFit", "RG", "fPowerSpectrum", "spatialCorr", "autocorr", "branching", "epilepsy", "Parkinsons", "depression", "autism", "Alzheimers", "schozophrenia", "stroke", "burstSuppression", "ptsd", "timeLineColor", "timeLineTopVBottom", "citeRate", "firstAuthor", "year", "journal", "title", "fullAuthorList", "bcolorR", "bcolorG", "bcolorB"];
opts.VariableTypes = ["string", "categorical", "categorical", "categorical", "categorical", "string", "string", "string", "string", "categorical", "string", "categorical", "string", "categorical", "string", "string", "categorical", "categorical", "categorical", "categorical", "categorical", "string", "string", "string", "categorical", "string", "categorical", "categorical", "categorical", "string", "string", "string", "string", "categorical", "string", "string", "string", "string", "categorical", "double", "string", "string", "categorical", "string", "string", "double", "double", "double"];

% Specify variable properties
opts = setvaropts(opts, ["mouse", "fish", "turtle", "leech", "crayfish", "retina", "CaImaging", "voltageImaging", "wholecellPatch", "IsingFit", "ARFit", "RG", "spatialCorr", "Parkinsons", "depression", "autism", "Alzheimers", "stroke", "burstSuppression", "ptsd", "timeLineColor", "firstAuthor", "year", "title", "fullAuthorList"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["mouse", "rat", "cat", "monkey", "human", "fish", "turtle", "leech", "crayfish", "culture", "retina", "spikes", "CaImaging", "LFP", "voltageImaging", "wholecellPatch", "EEG", "MEG", "BOLD", "avalanches", "LRTC", "IsingFit", "ARFit", "RG", "fPowerSpectrum", "spatialCorr", "autocorr", "branching", "epilepsy", "Parkinsons", "depression", "autism", "Alzheimers", "schozophrenia", "stroke", "burstSuppression", "ptsd", "timeLineColor", "timeLineTopVBottom", "firstAuthor", "year", "journal", "title", "fullAuthorList"], "EmptyFieldRule", "auto");

% Import the data
bibtab = readtable("G:\CriticalityReview2023\biblioTable012025.xlsx", opts, "UseExcel", false);
bibtab(1,:)=[];

speciesNames=opts.VariableNames;
speciesNames(12:end)=[];

measType=opts.VariableNames;
measType=measType(12:19);

analysisType=opts.VariableNames;
analysisType=analysisType(20:28);

disease=opts.VariableNames;
disease=disease(29:37);

%%%%%%%%%%%%%%%%% count species/preps %%%%%%%%%%%%%%%%%%%%%
nSpec=length(speciesNames);
counts=zeros(1,nSpec);
for i=1:nSpec
    counts(i) = sum(table2array(bibtab(:,i))=='x');
end
tbl = table(speciesNames,counts);

figure(10)
subplot(221)
p=piechart(tbl,"counts","speciesNames")
% Change the labels
p.Labels = p.Names + " (" + string(counts) + ")";

%%%%%%%%%%%%%%%%% count measurement modalities %%%%%%%%%%%%%%%%%%%%%
nSpec=length(measType);
counts=zeros(1,nSpec);
for i=1:nSpec
    counts(i) = sum(table2array(bibtab(:,i+11))=='x');
end
tbl = table(measType,counts);

subplot(222)
p=piechart(tbl,"counts","measType")
% Change the labels
p.Labels = p.Names + " (" + string(counts) + ")";

%%%%%%%%%%%%%%%%% count analysis type %%%%%%%%%%%%%%%%%%%%%
nSpec=length(analysisType);
counts=zeros(1,nSpec);
for i=1:nSpec
    counts(i) = sum(table2array(bibtab(:,i+19))=='x');
end
tbl = table(analysisType,counts);

subplot(223)
p=piechart(tbl,"counts","analysisType")
% Change the labels
p.Labels = p.Names + " (" + string(counts) + ")";

%%%%%%%%%%%%%%%%% count disease studies %%%%%%%%%%%%%%%%%%%%%
nSpec=length(disease);
counts=zeros(1,nSpec);
for i=1:nSpec
    counts(i) = sum(table2array(bibtab(:,i+28))=='x');
end
tbl = table(disease,counts);

subplot(224)
p=piechart(tbl,"counts","disease")
% Change the labels
p.Labels = p.Names + " (" + string(counts) + ")";