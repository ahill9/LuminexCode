%% Import Data
cd(expDir)
warning('off','MATLAB:xlswrite:AddSheet')
beadCounts=readtable(filename,'Sheet','BeadCounts','ReadRowNames',false);
beadCounts.Properties.VariableNames = regexprep(beadCounts.Properties.VariableNames,'[-_.]','');
sampleMFI = readtable(filename,'Sheet','SampleMFI','ReadRowNames',false);
sampleMFI.Properties.VariableNames = regexprep(sampleMFI.Properties.VariableNames,'[-_.]','');
sampleMFI{1:end,2:end}(beadCounts{1:end,2:end}<50)=NaN; %Remove MFIs with low bead counts

standardLabels = {'S0','S01','S02','S03','S04','S05','S06','S07','S08','S09','S10','S11','S12','S13','S14','S15','S16'}; %Possible sample names for standards
blankLabels = {'Blank','blank','BLANK'};
standardMFI = sampleMFI(ismember(sampleMFI{:,1},standardLabels),:);
[~, idx] = sort(standardMFI{:,1});
standardMFI = standardMFI(idx,:);
blankMFI = sampleMFI(ismember(sampleMFI{:,1},blankLabels),:);
sampleMFIsorted = sampleMFI(~ismember(sampleMFI{:,1},[standardLabels,blankLabels]),:);
[~, idx] = sort(sampleMFIsorted{:,1});
sampleMFIsorted = sampleMFIsorted(idx,:);
[~, sortIdx] = sort(sampleMFIsorted.Properties.VariableNames(2:end));
sampleMFIsorted = sampleMFIsorted(:,[1,sortIdx+1]);
standardMFIsorted = standardMFI(:,[1,sortIdx+1]);
blankMFIsorted = blankMFI(:,[1,sortIdx+1]);
blank0s=zeros(size(blankMFIsorted));
Cytokines = regexprep(sampleMFIsorted.Properties.VariableNames(2:end),'_','-');

%Check for outliers in standard MFIs
%*********************************************************************************
%To Do
%Grubbs test a nice balance of easy+appropriate, but supposed to be for >6
%replicates; see also http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-7-123
%**********************************************************************************

%Read in standard concentrations
standardTable = readtable('D:\OneDrive\Research\Protocols\Luminex\LuminexStandards.xlsx','Sheet',standardLot{1});
standardConc=standardTable(standardTable{:,1}==standardLot{2},:);
standardConc.Properties.VariableNames = regexprep(standardConc.Properties.VariableNames,'_','');
[sorted, sortIdx] = sort(standardConc.Properties.VariableNames(2:end));
standardConcSorted = standardConc(:,[1,sortIdx+1]);
%Calculate concentrations of standards
currDilution = 1; standards = unique(standardMFIsorted{:,1});
dilutionFactor=table(standardMFIsorted{:,1},ones(length(standardMFIsorted{:,1}),1));
for n = 1:length(standards)
    dilutionFactor{strcmp(dilutionFactor{:,1},standards(n)),2}=currDilution;
    currDilution = currDilution*standardDilutions;
end
allStandardsConc = dilutionFactor{:,2}*standardConcSorted{1,2:end};
allStandardsConc = [array2table(dilutionFactor{:,1}),array2table(allStandardsConc)];
allStandardsConc.Properties.VariableNames = ['Standard',standardTable.Properties.VariableNames(2:end)];

%Check that samples and standards match up
assert(isequal(standardMFIsorted.Properties.VariableNames(2:end),sampleMFIsorted.Properties.VariableNames(2:end)),'Analyte names do not match')
assert(isequal(standardMFIsorted.Properties.VariableNames(2:end),allStandardsConc.Properties.VariableNames(2:end)),'Analyte names do not match')
assert(isequal(standardMFIsorted.Properties.VariableNames(2:end),blankMFIsorted.Properties.VariableNames(2:end)),'Analyte names do not match')

%% Write useful stuff to file so you have it saved
sampleMFIsortedSTABLE=sampleMFIsorted;
standardMFIsortedSTABLE=standardMFIsorted;
currDateTime=datestr(now,'yyyy_mm_dd_HH_MM');
mkdir(cd,[analysisName,'_',currDateTime]); cd([analysisName,'_',currDateTime])
outFileName=[analysisName,'_',currDateTime,'.xlsx'];
writetable(sampleMFIsortedSTABLE,outFileName,'Sheet','sampleMFI')
writetable(standardMFIsortedSTABLE,outFileName,'Sheet','standardMFI')
writetable(allStandardsConc,outFileName,'Sheet','trueStandardsConc')


%% Method 1: Fit 5-parameter logistic curves to standards, omitting blanks
sampleConc_1 = sampleMFIsorted; sampleConc_1{1:end,2:end}=0;
A=zeros(1,size(sampleConc_1,2)-1); B=A; C=A; D=A; E=A; %initialize
flag=zeros(size(sampleConc_1,1),size(sampleConc_1,2)-1);
sampleMFIadj=sampleMFIsorted;
for c=1:size(sampleConc_1,2)-1
    [cf{c} gf{c}] = L5P(allStandardsConc{1:end,c+1},standardMFI{1:end,c+1}); 
    A(c)=cf{c}.A; B(c)=cf{c}.B; C(c)=cf{c}.C; D(c)=cf{c}.D; E(c)=cf{c}.E;
    flag(sampleMFIsorted{:,c+1}<A(c),c)=-1; sampleMFIadj{sampleMFIsorted{:,c+1}<A(c),c+1}=A(c); 
    flag(sampleMFIsorted{:,c+1}>D(c),c)=1; sampleMFIadj{sampleMFIsorted{:,c+1}>D(c),c+1}=0.99*D(c); 
    sampleConc_1{:,c+1} = L5Pinv(cf{c},sampleMFIadj{:,c+1});
    rsquares(c) = gf{c}.rsquare;
end
%Average replicates - assumes that sample replicates have identical names,
%omit NaNs
% **********************************************************************
% To Do: check for bad wells (e.g. no sample added? all MFIs lower than
% blank?)
% **********************************************************************
uniqueSampleNames = unique(sampleMFIsorted{:,1});
meanSampleConc_1 = zeros(length(uniqueSampleNames),size(sampleConc_1,2)-1);
for sampleIdx=1:length(uniqueSampleNames)
    meanSampleConc_1(sampleIdx,:) = mean(sampleMFIsorted{strcmp(sampleMFIsorted{:,1},uniqueSampleNames(sampleIdx)),2:end},1,'omitnan');
end
meanSampleConcTab_1 = sampleConc_1(1:length(uniqueSampleNames),:);
meanSampleConcTab_1{:,1} = uniqueSampleNames;
meanSampleConcTab_1{:,2:end} = meanSampleConc_1;

%% Method 2: Fit 5-parameter logistic curves to standards + blanks
sampleMFIsorted=sampleMFIsortedSTABLE;
standardMFIsorted=standardMFIsortedSTABLE;
sampleConc_2 = sampleMFIsorted; sampleConc_2{1:end,2:end}=0;
A=zeros(1,size(sampleConc_2,2)-1); B=A; C=A; D=A; E=A; %initialize
flag=zeros(size(sampleConc_2,1),size(sampleConc_2,2)-1);
sampleMFIadj=sampleMFIsorted;
for c=1:size(sampleConc_2,2)-1
    [cf{c},gf{c}] = L5P([allStandardsConc{1:end,c+1};blank0s(:,c+1)],[standardMFIsorted{1:end,c+1};blankMFIsorted{1:end,c+1}]); 
    A(c)=cf{c}.A; B(c)=cf{c}.B; C(c)=cf{c}.C; D(c)=cf{c}.D; E(c)=cf{c}.E;
    flag(sampleMFIsorted{:,c+1}<A(c),c)=-1; sampleMFIadj{sampleMFIsorted{:,c+1}<A(c),c+1}=A(c); 
    flag(sampleMFIsorted{:,c+1}>D(c),c)=1; sampleMFIadj{sampleMFIsorted{:,c+1}>D(c),c+1}=0.99*D(c); 
    sampleConc_2{:,c+1} = L5Pinv(cf{c},sampleMFIadj{:,c+1});
    rsquares(c) = gf{c}.rsquare;
end
%Average replicates - assumes that sample replicates have identical names,
%omit NaNs
% **********************************************************************
% To Do: check for bad wells (e.g. no sample added? all MFIs lower than
% blank?)
% **********************************************************************
uniqueSampleNames = unique(sampleMFIsorted{:,1});
meanSampleConc_2 = zeros(length(uniqueSampleNames),size(sampleConc_2,2)-1);
for sampleIdx=1:length(uniqueSampleNames)
    meanSampleConc_2(sampleIdx,:) = mean(sampleMFIsorted{strcmp(sampleMFIsorted{:,1},uniqueSampleNames(sampleIdx)),2:end},1,'omitnan');
end
meanSampleConcTab_2 = sampleConc_2(1:length(uniqueSampleNames),:);
meanSampleConcTab_2{:,1} = uniqueSampleNames;
meanSampleConcTab_2{:,2:end} = meanSampleConc_2;

%% Method 3: Subtract blank MFI, then fit 5-parameter logistic curve (ignore negatives)
sampleMFIsorted=sampleMFIsortedSTABLE;
standardMFIsorted=standardMFIsortedSTABLE;
sampleConc_3 = sampleMFIsorted; sampleConc_3{1:end,2:end}=0;
A=zeros(1,size(sampleConc_3,2)-1); B=A; C=A; D=A; E=A; %initialize
flag=zeros(size(sampleConc_3,1),size(sampleConc_3,2)-1);
sampleMFIadj=sampleMFIsorted;
sampleMFIsorted{:,2:end} = gsubtract(sampleMFIsorted{:,2:end},mean(blankMFIsorted{:,2:end}));
standardMFIsorted{:,2:end}=gsubtract(standardMFIsorted{:,2:end},mean(blankMFIsorted{:,2:end}));
for c=1:size(sampleConc_3,2)-1
    [cf{c},gf{c}] = L5P([allStandardsConc{1:end,c+1};blank0s(:,c+1)],[standardMFIsorted{1:end,c+1};blankMFIsorted{1:end,c+1}]); 
    A(c)=cf{c}.A; B(c)=cf{c}.B; C(c)=cf{c}.C; D(c)=cf{c}.D; E(c)=cf{c}.E;
    flag(sampleMFIsorted{:,c+1}<A(c),c)=-1; sampleMFIadj{sampleMFIsorted{:,c+1}<A(c),c+1}=A(c); 
    flag(sampleMFIsorted{:,c+1}>D(c),c)=1; sampleMFIadj{sampleMFIsorted{:,c+1}>D(c),c+1}=0.99*D(c); 
    sampleConc_3{:,c+1} = L5Pinv(cf{c},sampleMFIadj{:,c+1});
    rsquares(c) = gf{c}.rsquare;
end
%Average replicates - assumes that sample replicates have identical names,
%omit NaNs
% **********************************************************************
% To Do: check for bad wells (e.g. no sample added? all MFIs lower than
% blank?)
% **********************************************************************
uniqueSampleNames = unique(sampleMFIsorted{:,1});
meanSampleConc_3 = zeros(length(uniqueSampleNames),size(sampleConc_3,2)-1);
for sampleIdx=1:length(uniqueSampleNames)
    meanSampleConc_3(sampleIdx,:) = mean(sampleMFIsorted{strcmp(sampleMFIsorted{:,1},uniqueSampleNames(sampleIdx)),2:end},1,'omitnan');
end
meanSampleConcTab_3 = sampleConc_3(1:length(uniqueSampleNames),:);
meanSampleConcTab_3{:,1} = uniqueSampleNames;
meanSampleConcTab_3{:,2:end} = meanSampleConc_3;
