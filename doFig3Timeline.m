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

citeRate=bibtab.citeRate + 0.0001; %citions/year
bcol=[bibtab.bcolorR bibtab.bcolorG bibtab.bcolorB]; %color code
bcol(isnan(bcol))=255;
bcol(bcol==300)=255;
bcol=bcol/255;

tvb=bibtab.timeLineTopVBottom=='t'; %top v bottom indicator
pyear=str2double(bibtab.year);
n=length(pyear);
pndx=1:n;

strexp=0.4;
vpos=(citeRate.^strexp);

figure(100)
for y=1998:2024
    yx=find(diff([0; pyear]<y)<0,1)-0.5;
    plot([yx yx],[-7 7],'Color',[1 1 1]*0.7)
    hold on;
    if mod(y,5)==0; text(yx,-7,num2str(y)); end
end
for xdiv=[1 10:10:1.5*max(citeRate)]
    plot([1 n],[1 1]*xdiv^strexp,'Color',[1 1 1]*0.7)
    plot([1 n],-[1 1]*xdiv^strexp,'Color',[1 1 1]*0.7)
end

toplist=find(tvb)'; nt=length(toplist);
botlist=find(~tvb)'; nb=length(botlist);

for i=1:nt
plot(pndx(toplist(i)),vpos(toplist(i)),'o', ...
    'markersize', (citeRate(toplist(i))).^0.5*4, ...
    'MarkerFaceColor',bcol(toplist(i),:), ...
    'MarkerEdgeColor',[1 1 1]*0.8)
hold on
end


for i=1:nb
plot(pndx(botlist(i)),-vpos(botlist(i)),'o', ...
    'markersize', (citeRate(botlist(i))).^0.5*4, ...
    'MarkerFaceColor',bcol(botlist(i),:), ...
    'MarkerEdgeColor',[1 1 1]*0.8)
hold on
end
hold off

% scatter(pndx(tvb),vpos(tvb),25*ones(1,sum(tvb)),bcol(tvb,:))