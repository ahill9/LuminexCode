%Analysis of Proliferative samples July 2016
clear; clc; close all

%% Input
%Must modify this section for each assay:
analysisName = 'ProlifPrimary_CellLines_Compare';
expDir='D:\OneDrive\Research\0C Triculture Christi\Exp7\';
filename='ProlifPrimary_CellLines_Compare_revised.xlsx';
standardLot = {'27plex',64020782};
standardDilutions = 1/3;

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
warning('off','L5P:IgnorningNansAndInfs')
sampleMFIsorted=sampleMFIsortedSTABLE;
standardMFIsorted=standardMFIsortedSTABLE;
sampleConc_1 = sampleMFIsorted; sampleConc_1{1:end,2:end}=0;
A=zeros(1,size(sampleConc_1,2)-1); B=A; C=A; D=A; E=A; %initialize
flag=zeros(size(sampleConc_1,1),size(sampleConc_1,2)-1);
sampleMFIadj=sampleMFIsorted;
for c=1:size(sampleConc_1,2)-1
    [cf{c} gf{c}] = L5P(allStandardsConc{1:end,c+1},standardMFIsorted{1:end,c+1}); 
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


%% Specific to this experiment
meanSampleConcTab = meanSampleConcTab_1; %Select which method to use for plots
figuresFilename = 'exp7figures_3.ps';
[~, ~, ~, ~, ~, bhpvals0188day7v15] = makeLuminexSuperbarPlot(meanSampleConcTab,'0188 D7 1','0188 D15 1','0188 Day 7 vs. Day 15 No Stim','paired');
print(gcf, '-dpsc', figuresFilename);
[~, ~, ~, ~, ~, bhpvals0190day7v15] =makeLuminexSuperbarPlot(meanSampleConcTab,'0190 D8 1','0190 D16 1','0190 Day 8 vs. Day 16 No Stim','paired');
print(gcf, '-dpsc', figuresFilename,'-append')
[~, ~, ~, ~, ~, bhpvals0191day7v15] =makeLuminexSuperbarPlot(meanSampleConcTab,'0191 D7 1[1-5]','0191 D15 1[1-5]','0191 Day 7 vs. Day 15 No Stim','paired');
print(gcf, '-dpsc', figuresFilename,'-append')

[~, ~, ~, ~, ~, bhpvals0188preVpostStim] =makeLuminexSuperbarPlot(meanSampleConcTab,'0188 D7 2','0188 D15 2','0188 Day 7 Pre-Stim vs. Day 15 Post-Stim','paired');
print(gcf, '-dpsc', figuresFilename,'-append');
[~, ~, ~, ~, ~, bhpvals0190preVpostStim] =makeLuminexSuperbarPlot(meanSampleConcTab,'0190 D8 2','0190 D16 2','0190 Day 8 Pre-Stim vs. Day 16 Post-Stim','paired');
print(gcf, '-dpsc', figuresFilename,'-append');
[~, ~, ~, ~, ~, bhpvals0191preVpostStim] =makeLuminexSuperbarPlot(meanSampleConcTab,'0191 D7 2[1-5]','0191 D15 2[1-5]','0191 Day 7 Pre-Stim vs. Day 15 Post-Stim','paired');
print(gcf, '-dpsc', figuresFilename,'-append');

[~, ~, ~, ~, ~, bhpvals0188stimVunstim] =makeLuminexSuperbarPlot(meanSampleConcTab,'0188 D15 2','0188 D15 1','0188 Day 15 No Stim vs. Stim','unpaired');
print(gcf, '-dpsc', figuresFilename,'-append');
[~, ~, ~, ~, ~, bhpvals0190stimVunstim] =makeLuminexSuperbarPlot(meanSampleConcTab,'0190 D16 2','0190 D16 1','0190 Day 16 No Stim vs. Stim','unpaired');
print(gcf, '-dpsc', figuresFilename,'-append');
[~, ~, ~, ~, ~, bhpvals0191stimVunstim] =makeLuminexSuperbarPlot(meanSampleConcTab,'0191 D15 2[1-5]','0191 D15 1[1-5]','0191 Day 15 No Stim vs. Stim','unpaired');
print(gcf, '-dpsc', figuresFilename,'-append');

makeLuminexSuperbarPlot(meanSampleConcTab,'CL D6 [5-9] A','CL D14 [5-9] A','CL Day 6 Pre-Stim vs. Day 14 Post-Stim','unpaired');
print(gcf, '-dpsc', figuresFilename,'-append');
makeLuminexSuperbarPlot(meanSampleConcTab,'CL D6 5[1-4] A','0188 D7 1','CL Day 6 Pre-Stim vs. 0188 Day 7 Pre-Stim','unpaired');
print(gcf, '-dpsc', figuresFilename,'-append');
makeLuminexSuperbarPlot(meanSampleConcTab,'CL D6 5[1-4] A','0190 D8 1','CL Day 6 Pre-Stim vs. 0190 D8 Pre-Stim','unpaired');
print(gcf, '-dpsc', figuresFilename,'-append');
makeLuminexSuperbarPlot(meanSampleConcTab,'CL D6 5[1-4] A','0191 D7 1','CL Day 6 Pre-Stim vs. 0191 Day 7 Pre-Stim','unpaired');
print(gcf, '-dpsc', figuresFilename,'-append')

% All Day 6-8 samples, Cell Lines and Primaries (unstim and pre-stim)
    sample1temp = regexp(uniqueSampleNames,'CL D6 [5-6][1-6]');
    sample1 = meanSampleConcTab(cellfun(@length,sample1temp(:))>0,:);
    sample2temp = regexp(uniqueSampleNames,'CL D6 [1-9] A');
    sample2 = meanSampleConcTab(cellfun(@length,sample2temp(:))>0,:);
    sample3temp = regexp(uniqueSampleNames,'0188 D7 [1-2]');
    sample3 = meanSampleConcTab(cellfun(@length,sample3temp(:))>0,:);
    sample4temp = regexp(uniqueSampleNames,'0190 D8 [1-2]');
    sample4 = meanSampleConcTab(cellfun(@length,sample4temp(:))>0,:);
    sample5temp = regexp(uniqueSampleNames,'0191 D7 [1-2]');
    sample5 = meanSampleConcTab(cellfun(@length,sample5temp(:))>0,:);
    [h,p]=ttest2([mean(sample1{:,2:end},1);mean(sample2{:,2:end},1)],[mean(sample3{:,2:end},1);mean(sample4{:,2:end},1);mean(sample5{:,2:end},1)]);
    bhpvals = mafdr(p,'BHFDR',true);
    pvals4plot = zeros(2*length(bhpvals),2*length(bhpvals))/0;
    for n = 1:size(sample1,2)-1
        pvals4plot(n,27+n) = bhpvals(n); pvals4plot(27+n,n) = bhpvals(n);
    end
    figure
    [hb, ~, ~, ~, ~] = superbar(1:27,[mean(sample1{:,2:end});mean(sample2{:,2:end});mean(sample3{:,2:end});mean(sample4{:,2:end});mean(sample5{:,2:end})]',...
        'E',[std(sample1{:,2:end});std(sample2{:,2:end});std(sample3{:,2:end});std(sample4{:,2:end});std(sample5{:,2:end})]',...
        'BarFaceColor',permute([0.2, 0.8, 0.8; 0.4,0.6,0.6; .8, 0.2, 0.2; 0.7,0.3,0.3; 0.6,0.4,0.4],[3,1,2]),...
        'ErrorbarColor','k');
    xlabel('Cytokine'); ylabel('Concentration (pg/mL)'); title('Cell Lines and Primaries All Days 6-8')
    set(gca,'XTick',1:27,'XTickLabel',Cytokines,'XTickLabelRotation',90);
    legend(hb(1,1:5),{'Cell Lines 1','Cell Lines 2','Primary 0188','Primary 0190','Primary 0191'})
print(gcf, '-dpsc', figuresFilename,'-append')


%% Specific cytokines

%IL-10
figure
sample1temp = regexp(uniqueSampleNames,'0188 D7 2');
earlyNoStim0188 = meanSampleConcTab.IL10(cellfun(@length,sample1temp(:))>0,:);
sample1temp = regexp(uniqueSampleNames,'0188 D7 1');
earlyUnStim0188 = meanSampleConcTab.IL10(cellfun(@length,sample1temp(:))>0,:);
sample1temp = regexp(uniqueSampleNames,'0188 D15 2');
lateNoStim0188 = meanSampleConcTab.IL10(cellfun(@length,sample1temp(:))>0,:);
sample1temp = regexp(uniqueSampleNames,'0188 D15 1');
lateStim0188 = meanSampleConcTab.IL10(cellfun(@length,sample1temp(:))>0,:);

sample2temp = regexp(uniqueSampleNames,'0190 D8 2');
earlyNoStim0190 = meanSampleConcTab.IL10(cellfun(@length,sample2temp(:))>0,:);
sample2temp = regexp(uniqueSampleNames,'0190 D8 1');
earlyUnStim0190 = meanSampleConcTab.IL10(cellfun(@length,sample2temp(:))>0,:);
sample2temp = regexp(uniqueSampleNames,'0190 D16 2');
lateNoStim0190 = meanSampleConcTab.IL10(cellfun(@length,sample2temp(:))>0,:);
sample2temp = regexp(uniqueSampleNames,'0190 D16 1');
lateStim0190 = meanSampleConcTab.IL10(cellfun(@length,sample2temp(:))>0,:);

sample3temp = regexp(uniqueSampleNames,'0191 D7 2');
earlyNoStim0191= meanSampleConcTab.IL10(cellfun(@length,sample3temp(:))>0,:);
sample3temp = regexp(uniqueSampleNames,'0191 D7 1');
earlyUnStim0191= meanSampleConcTab.IL10(cellfun(@length,sample3temp(:))>0,:);
sample3temp = regexp(uniqueSampleNames,'0191 D15 2');
lateNoStim0191= meanSampleConcTab.IL10(cellfun(@length,sample3temp(:))>0,:);
sample3temp = regexp(uniqueSampleNames,'0191 D15 1');
lateStim0191= meanSampleConcTab.IL10(cellfun(@length,sample3temp(:))>0,:);

Labels = {'0188','0190','0191'};

IL10pvals = NaN*ones(12,12);
IL10pvals(4,10)=bhpvals0191preVpostStim(strcmp(Cytokines,'IL10'));
IL10pvals(10,4)=IL10pvals(4,10);
IL10pvals(7,10)=bhpvals0188stimVunstim(strcmp(Cytokines,'IL10'));
IL10pvals(10,7)=IL10pvals(7,10);
IL10pvals(5,11)=bhpvals0191preVpostStim(strcmp(Cytokines,'IL10'));
IL10pvals(11,5)=IL10pvals(5,11);
IL10pvals(8,11)=bhpvals0190stimVunstim(strcmp(Cytokines,'IL10'));
IL10pvals(11,8)=IL10pvals(8,11);
IL10pvals(6,12)=bhpvals0191preVpostStim(strcmp(Cytokines,'IL10'));
IL10pvals(12,6)=IL10pvals(6,12);
IL10pvals(9,12)=bhpvals0191stimVunstim(strcmp(Cytokines,'IL10'));
IL10pvals(12,9)=IL10pvals(9,12);
IL10pvals(IL10pvals>=0.05)=NaN; %Prevent non-significant ones from showing up

[hb, ~, ~, ~, ~] = superbar(1:3,[mean(earlyNoStim0188),mean(earlyUnStim0188),mean(lateNoStim0188),mean(lateStim0188);
    mean(earlyNoStim0190),mean(earlyUnStim0190),mean(lateNoStim0190),mean(lateStim0190);
    mean(earlyNoStim0191),mean(earlyUnStim0191),mean(lateNoStim0191),mean(lateStim0191)], ...
    'E',[std(earlyNoStim0188),std(earlyUnStim0188),std(lateNoStim0188),std(lateStim0188);
    std(earlyNoStim0190),std(earlyUnStim0190),std(lateNoStim0190),std(lateStim0190);
    std(earlyNoStim0191),std(earlyUnStim0191),std(lateNoStim0191),std(lateStim0191)],...
    'P',IL10pvals,...
    'BarFaceColor',permute([0.2, 0.8, 0.8; 0.4,0.6,0.6; .6,.6,.6; .8, 0.2, 0.2] ,[3,1,2]),...
    'ErrorbarColor','k');
xlabel('IL-10'); ylabel('Concentration (pg/mL)'); title('IL-10 Primaries')
set(gca,'XTick',1:3,'XTickLabel',Labels,'XTickLabelRotation',90);
axis([0,4,0,3200])
legend(hb(1,1:4),{'Early, No Stim','Early, Pre-Stim','Late, No Stim','Late, Post-Stim'})

%IL-15
figure
sample1temp = regexp(uniqueSampleNames,'0188 D7 2');
earlyNoStim0188 = meanSampleConcTab.IL15(cellfun(@length,sample1temp(:))>0,:);
sample1temp = regexp(uniqueSampleNames,'0188 D7 1');
earlyUnStim0188 = meanSampleConcTab.IL15(cellfun(@length,sample1temp(:))>0,:);
sample1temp = regexp(uniqueSampleNames,'0188 D15 2');
lateNoStim0188 = meanSampleConcTab.IL15(cellfun(@length,sample1temp(:))>0,:);
sample1temp = regexp(uniqueSampleNames,'0188 D15 1');
lateStim0188 = meanSampleConcTab.IL15(cellfun(@length,sample1temp(:))>0,:);

sample2temp = regexp(uniqueSampleNames,'0190 D8 2');
earlyNoStim0190 = meanSampleConcTab.IL15(cellfun(@length,sample2temp(:))>0,:);
sample2temp = regexp(uniqueSampleNames,'0190 D8 1');
earlyUnStim0190 = meanSampleConcTab.IL15(cellfun(@length,sample2temp(:))>0,:);
sample2temp = regexp(uniqueSampleNames,'0190 D16 2');
lateNoStim0190 = meanSampleConcTab.IL15(cellfun(@length,sample2temp(:))>0,:);
sample2temp = regexp(uniqueSampleNames,'0190 D16 1');
lateStim0190 = meanSampleConcTab.IL15(cellfun(@length,sample2temp(:))>0,:);

sample3temp = regexp(uniqueSampleNames,'0191 D7 2');
earlyNoStim0191= meanSampleConcTab.IL15(cellfun(@length,sample3temp(:))>0,:);
sample3temp = regexp(uniqueSampleNames,'0191 D7 1');
earlyUnStim0191= meanSampleConcTab.IL15(cellfun(@length,sample3temp(:))>0,:);
sample3temp = regexp(uniqueSampleNames,'0191 D15 2');
lateNoStim0191= meanSampleConcTab.IL15(cellfun(@length,sample3temp(:))>0,:);
sample3temp = regexp(uniqueSampleNames,'0191 D15 1');
lateStim0191= meanSampleConcTab.IL15(cellfun(@length,sample3temp(:))>0,:);

Labels = {'0188','0190','0191'};


IL15pvals = NaN*ones(12,12);
IL15pvals(4,10)=bhpvals0191preVpostStim(strcmp(Cytokines,'IL15'));
IL15pvals(10,4)=IL15pvals(4,10);
IL15pvals(7,10)=bhpvals0188stimVunstim(strcmp(Cytokines,'IL15'));
IL15pvals(10,7)=IL15pvals(7,10);
IL15pvals(5,11)=bhpvals0191preVpostStim(strcmp(Cytokines,'IL15'));
IL15pvals(11,5)=IL15pvals(5,11);
IL15pvals(8,11)=bhpvals0190stimVunstim(strcmp(Cytokines,'IL15'));
IL15pvals(11,8)=IL15pvals(8,11);
IL15pvals(6,12)=bhpvals0191preVpostStim(strcmp(Cytokines,'IL15'));
IL15pvals(12,6)=IL15pvals(6,12);
IL15pvals(9,12)=bhpvals0191stimVunstim(strcmp(Cytokines,'IL15'));
IL15pvals(12,9)=IL15pvals(9,12);
IL15pvals(IL15pvals>=0.05)=NaN; %Prevent non-significant ones from showing up

[hb, ~, ~, ~, ~] = superbar(1:3,[mean(earlyNoStim0188),mean(earlyUnStim0188),mean(lateNoStim0188),mean(lateStim0188);
    mean(earlyNoStim0190),mean(earlyUnStim0190),mean(lateNoStim0190),mean(lateStim0190);
    mean(earlyNoStim0191),mean(earlyUnStim0191),mean(lateNoStim0191),mean(lateStim0191)], ...
    'E',[std(earlyNoStim0188),std(earlyUnStim0188),std(lateNoStim0188),std(lateStim0188);
    std(earlyNoStim0190),std(earlyUnStim0190),std(lateNoStim0190),std(lateStim0190);
    std(earlyNoStim0191),std(earlyUnStim0191),std(lateNoStim0191),std(lateStim0191)],...
    'P',IL15pvals,...
    'BarFaceColor',permute([0.2, 0.8, 0.8; 0.4,0.6,0.6; .6,.6,.6; .8, 0.2, 0.2] ,[3,1,2]),...
    'ErrorbarColor','k');
xlabel('IL-15'); ylabel('Concentration (pg/mL)'); title('IL-15 Primaries')
set(gca,'XTick',1:3,'XTickLabel',Labels,'XTickLabelRotation',90);
axis([0,4,0,3600])
%legend(hb(1,1:4),{'Early, No Stim','Early, Pre-Stim','Late, No Stim','Late, Post-Stim'})


%IL-8
figure
sample1temp = regexp(uniqueSampleNames,'0188 D7 2');
earlyNoStim0188 = meanSampleConcTab.IL8(cellfun(@length,sample1temp(:))>0,:);
sample1temp = regexp(uniqueSampleNames,'0188 D7 1');
earlyUnStim0188 = meanSampleConcTab.IL8(cellfun(@length,sample1temp(:))>0,:);
sample1temp = regexp(uniqueSampleNames,'0188 D15 2');
lateNoStim0188 = meanSampleConcTab.IL8(cellfun(@length,sample1temp(:))>0,:);
sample1temp = regexp(uniqueSampleNames,'0188 D15 1');
lateStim0188 = meanSampleConcTab.IL8(cellfun(@length,sample1temp(:))>0,:);

sample2temp = regexp(uniqueSampleNames,'0190 D8 2');
earlyNoStim0190 = meanSampleConcTab.IL8(cellfun(@length,sample2temp(:))>0,:);
sample2temp = regexp(uniqueSampleNames,'0190 D8 1');
earlyUnStim0190 = meanSampleConcTab.IL8(cellfun(@length,sample2temp(:))>0,:);
sample2temp = regexp(uniqueSampleNames,'0190 D16 2');
lateNoStim0190 = meanSampleConcTab.IL8(cellfun(@length,sample2temp(:))>0,:);
sample2temp = regexp(uniqueSampleNames,'0190 D16 1');
lateStim0190 = meanSampleConcTab.IL8(cellfun(@length,sample2temp(:))>0,:);

sample3temp = regexp(uniqueSampleNames,'0191 D7 2');
earlyNoStim0191= meanSampleConcTab.IL8(cellfun(@length,sample3temp(:))>0,:);
sample3temp = regexp(uniqueSampleNames,'0191 D7 1');
earlyUnStim0191= meanSampleConcTab.IL8(cellfun(@length,sample3temp(:))>0,:);
sample3temp = regexp(uniqueSampleNames,'0191 D15 2');
lateNoStim0191= meanSampleConcTab.IL8(cellfun(@length,sample3temp(:))>0,:);
sample3temp = regexp(uniqueSampleNames,'0191 D15 1');
lateStim0191= meanSampleConcTab.IL8(cellfun(@length,sample3temp(:))>0,:);

Labels = {'0188','0190','0191'};


IL8pvals = NaN*ones(12,12);
IL8pvals(4,10)=bhpvals0191preVpostStim(strcmp(Cytokines,'IL8'));
IL8pvals(10,4)=IL8pvals(4,10);
IL8pvals(7,10)=bhpvals0188stimVunstim(strcmp(Cytokines,'IL8'));
IL8pvals(10,7)=IL8pvals(7,10);
IL8pvals(5,11)=bhpvals0191preVpostStim(strcmp(Cytokines,'IL8'));
IL8pvals(11,5)=IL8pvals(5,11);
IL8pvals(8,11)=bhpvals0190stimVunstim(strcmp(Cytokines,'IL8'));
IL8pvals(11,8)=IL8pvals(8,11);
IL8pvals(6,12)=bhpvals0191preVpostStim(strcmp(Cytokines,'IL8'));
IL8pvals(12,6)=IL8pvals(6,12);
IL8pvals(9,12)=bhpvals0191stimVunstim(strcmp(Cytokines,'IL8'));
IL8pvals(12,9)=IL8pvals(9,12);
IL8pvals(IL8pvals>=0.05)=NaN; %Prevent non-significant ones from showing up

[hb, ~, ~, ~, ~] = superbar(1:3,[mean(earlyNoStim0188),mean(earlyUnStim0188),mean(lateNoStim0188),mean(lateStim0188);
    mean(earlyNoStim0190),mean(earlyUnStim0190),mean(lateNoStim0190),mean(lateStim0190);
    mean(earlyNoStim0191),mean(earlyUnStim0191),mean(lateNoStim0191),mean(lateStim0191)], ...
    'E',[std(earlyNoStim0188),std(earlyUnStim0188),std(lateNoStim0188),std(lateStim0188);
    std(earlyNoStim0190),std(earlyUnStim0190),std(lateNoStim0190),std(lateStim0190);
    std(earlyNoStim0191),std(earlyUnStim0191),std(lateNoStim0191),std(lateStim0191)],...
    'P',IL8pvals,...
    'BarFaceColor',permute([0.2, 0.8, 0.8; 0.4,0.6,0.6; .6,.6,.6; .8, 0.2, 0.2] ,[3,1,2]),...
    'ErrorbarColor','k');
xlabel('IL-10'); ylabel('Concentration (pg/mL)'); title('IL-8 Primaries')
set(gca,'XTick',1:3,'XTickLabel',Labels,'XTickLabelRotation',90);
legend(hb(1,1:4),{'Early, No Stim','Early, Pre-Stim','Late, No Stim','Late, Post-Stim'})

% %IL-13
% figure
% sample1temp = regexp(uniqueSampleNames,'0188 D7 2');
% earlyNoStim0188 = meanSampleConcTab.IL13(cellfun(@length,sample1temp(:))>0,:);
% sample1temp = regexp(uniqueSampleNames,'0188 D7 1');
% earlyUnStim0188 = meanSampleConcTab.IL13(cellfun(@length,sample1temp(:))>0,:);
% sample1temp = regexp(uniqueSampleNames,'0188 D15 2');
% lateNoStim0188 = meanSampleConcTab.IL13(cellfun(@length,sample1temp(:))>0,:);
% sample1temp = regexp(uniqueSampleNames,'0188 D15 1');
% lateStim0188 = meanSampleConcTab.IL13(cellfun(@length,sample1temp(:))>0,:);
% 
% sample2temp = regexp(uniqueSampleNames,'0190 D8 2');
% earlyNoStim0190 = meanSampleConcTab.IL13(cellfun(@length,sample2temp(:))>0,:);
% sample2temp = regexp(uniqueSampleNames,'0190 D8 1');
% earlyUnStim0190 = meanSampleConcTab.IL13(cellfun(@length,sample2temp(:))>0,:);
% sample2temp = regexp(uniqueSampleNames,'0190 D16 2');
% lateNoStim0190 = meanSampleConcTab.IL13(cellfun(@length,sample2temp(:))>0,:);
% sample2temp = regexp(uniqueSampleNames,'0190 D16 1');
% lateStim0190 = meanSampleConcTab.IL13(cellfun(@length,sample2temp(:))>0,:);
% 
% sample3temp = regexp(uniqueSampleNames,'0191 D7 2');
% earlyNoStim0191= meanSampleConcTab.IL13(cellfun(@length,sample3temp(:))>0,:);
% sample3temp = regexp(uniqueSampleNames,'0191 D7 1');
% earlyUnStim0191= meanSampleConcTab.IL13(cellfun(@length,sample3temp(:))>0,:);
% sample3temp = regexp(uniqueSampleNames,'0191 D15 2');
% lateNoStim0191= meanSampleConcTab.IL13(cellfun(@length,sample3temp(:))>0,:);
% sample3temp = regexp(uniqueSampleNames,'0191 D15 1');
% lateStim0191= meanSampleConcTab.IL13(cellfun(@length,sample3temp(:))>0,:);
% 
% Labels = {'0188','0190','0191'};
% 
% [hb, ~, ~, ~, ~] = superbar(1:3,[mean(earlyNoStim0188),mean(earlyUnStim0188),mean(lateNoStim0188),mean(lateStim0188);
%     mean(earlyNoStim0190),mean(earlyUnStim0190),mean(lateNoStim0190),mean(lateStim0190);
%     mean(earlyNoStim0191),mean(earlyUnStim0191),mean(lateNoStim0191),mean(lateStim0191)], ...
%     'E',[std(earlyNoStim0188),std(earlyUnStim0188),std(lateNoStim0188),std(lateStim0188);
%     std(earlyNoStim0190),std(earlyUnStim0190),std(lateNoStim0190),std(lateStim0190);
%     std(earlyNoStim0191),std(earlyUnStim0191),std(lateNoStim0191),std(lateStim0191)],...
%     'BarFaceColor',permute([0.2, 0.8, 0.8; 0.4,0.6,0.6; .6,.6,.6; .8, 0.2, 0.2] ,[3,1,2]),...
%     'ErrorbarColor','k');
% xlabel('Subjects'); ylabel('Concentration (pg/mL)'); title('IL-13 Primaries')
% set(gca,'XTick',1:3,'XTickLabel',Labels,'XTickLabelRotation',90);
% legend(hb(1,1:4),{'Early, No Stim','Early, Pre-Stim','Late, No Stim','Late, Post-Stim'})
% 
% 
