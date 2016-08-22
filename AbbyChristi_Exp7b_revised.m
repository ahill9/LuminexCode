%Analysis of Proliferative samples July 2016
clear; clc; close all

%% Input
%Must modify this section for each assay:
analysisName = 'ProlifPrimary_CellLines_Compare';
expDir='D:\OneDrive\Research\0C Triculture Christi\Exp7\';
filename='ProlifPrimary_CellLines_Compare_revised.xlsx';
standardLot = {'27plex',64020782};
standardDilutions = 1/3;

LuminexAnalysisGeneral %Calculates concentrations

%% Specific to this experiment
meanSampleConcTab = meanSampleConcTab_1; %Select which fitting method to use for plots
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
    pvals4plot = zeros(5*length(bhpvals),5*length(bhpvals))/0;
    for n = 1:size(sample1,2)-1
        pvals4plot(n,3*27+n) = bhpvals(n); pvals4plot(3*27+n,n) = bhpvals(n);
    end
    figure
    [hb, ~, ~, ~, ~] = superbar(1:27,[mean(sample1{:,2:end});mean(sample2{:,2:end});mean(sample3{:,2:end});mean(sample4{:,2:end});mean(sample5{:,2:end})]',...
        'E',[std(sample1{:,2:end});std(sample2{:,2:end});std(sample3{:,2:end});std(sample4{:,2:end});std(sample5{:,2:end})]',...
        'P',pvals4plot,...
        'BarFaceColor',permute([0.2, 0.8, 0.8; 0.4,0.6,0.6; .8, 0.2, 0.2; 0.7,0.3,0.3; 0.6,0.4,0.4],[3,1,2]),...
        'ErrorbarColor','k');
    xlabel('Cytokine'); ylabel('Concentration (pg/mL)'); title('Cell Lines and Primaries All Days 6-8')
    set(gca,'XTick',1:27,'XTickLabel',Cytokines,'XTickLabelRotation',90);
    legend(hb(1,1:5),{'Cell Lines 1','Cell Lines 2','Primary 0188','Primary 0190','Primary 0191'})
print(gcf, '-dpsc', figuresFilename,'-append')

%% Comparing only equivalent cell line experiment, with all replicates of each donor
    sample1temp = regexp(uniqueSampleNames,'CL D6 [5-6][1-6]');
    sample1 = meanSampleConcTab(cellfun(@length,sample1temp(:))>0,:);
    
    sample3temp = regexp(uniqueSampleNames,'0188 D7 [1-2]');
    sample3 = meanSampleConcTab(cellfun(@length,sample3temp(:))>0,:);
    sample4temp = regexp(uniqueSampleNames,'0190 D8 [1-2]');
    sample4 = meanSampleConcTab(cellfun(@length,sample4temp(:))>0,:);
    sample5temp = regexp(uniqueSampleNames,'0191 D7 [1-2]');
    sample5 = meanSampleConcTab(cellfun(@length,sample5temp(:))>0,:);
    [h,p1]=ttest2(sample1{:,2:end},sample3{:,2:end},'vartype','unequal');
    bhpvals1 = mafdr(p1,'BHFDR',true);
    [h,p2]=ttest2(sample1{:,2:end},sample4{:,2:end},'vartype','unequal');
    bhpvals2 = mafdr(p2,'BHFDR',true);
    [h,p3]=ttest2(sample1{:,2:end},sample5{:,2:end},'vartype','unequal');
    bhpvals3 = mafdr(p3,'BHFDR',true);
    pvals4plot = zeros(4*length(bhpvals),4*length(bhpvals))/0;
    for n = 1:27
        pvals4plot(n,27+n)=bhpvals1(n); pvals4plot(27+n,n)=bhpvals1(n);
        pvals4plot(n,2*27+n)=bhpvals2(n); pvals4plot(2*27+n,n)=bhpvals2(n);
        pvals4plot(n,3*27+n)=bhpvals3(n); pvals4plot(3*27+n,n)=bhpvals3(n);
    end 
    pvals4plot(pvals4plot>0.05)=NaN;

    figure
    [hb, ~, ~, ~, ~] = superbar(1:27,[mean(sample1{:,2:end});mean(sample3{:,2:end});mean(sample4{:,2:end});mean(sample5{:,2:end})]',...
        'E',[std(sample1{:,2:end});std(sample3{:,2:end});std(sample4{:,2:end});std(sample5{:,2:end})]',...
        'P',pvals4plot,'PLineColor',[1,1,1],'PLineWidth',0.01,...
        'BarFaceColor',permute([0.2, 0.8, 0.8; .8, 0.2, 0.2; 0.7,0.3,0.3; 0.6,0.4,0.4],[3,1,2]),...
        'ErrorbarColor','k','PStarShowNS',false,'PStarShowGT',false,'PStarColor',[.2,.8,.8],'PStarOffset',0);
    xlabel('Cytokine'); ylabel('Concentration (pg/mL)'); title('Cell Lines and Primaries All Days 6-8')
    set(gca,'XTick',1:27,'XTickLabel',Cytokines,'XTickLabelRotation',90);
    legend(hb(1,1:4),{'Cell Lines 1','Primary 0188','Primary 0190','Primary 0191'})
    print(gcf, '-dpsc', figuresFilename,'-append')
    
% No pvals
    figure
    [hb, ~, ~, ~, ~] = superbar(1:27,[mean(sample1{:,2:end});mean(sample3{:,2:end});mean(sample4{:,2:end});mean(sample5{:,2:end})]',...
        'E',[std(sample1{:,2:end});std(sample3{:,2:end});std(sample4{:,2:end});std(sample5{:,2:end})]',...
        'BarFaceColor',permute([0.2, 0.8, 0.8; .8, 0.2, 0.2; 0.7,0.3,0.3; 0.6,0.4,0.4],[3,1,2]),...
        'ErrorbarColor','k');
    xlabel('Cytokine'); ylabel('Concentration (pg/mL)'); title('Cell Lines and Primaries All Days 6-8')
    set(gca,'XTick',1:27,'XTickLabel',Cytokines,'XTickLabelRotation',90);
    legend(hb(1,1:4),{'Cell Lines 1','Primary 0188','Primary 0190','Primary 0191'})
    print(gcf, '-dpsc', figuresFilename,'-append')

 %Split into 3 panels based on magnitude
    [~,max_data_idx] = sort(max([mean(sample1{:,2:end});mean(sample3{:,2:end});mean(sample4{:,2:end});mean(sample5{:,2:end})]));
    data_A_idx = max_data_idx(1:9);
    data_B_idx = max_data_idx(10:18);
    data_C_idx = max_data_idx(19:27);
    figure
    subplot(1,3,1)
    [hb, ~, ~, ~, ~] = superbar(1:9,[mean(sample1{:,data_A_idx+1});mean(sample3{:,data_A_idx+1});mean(sample4{:,data_A_idx+1});mean(sample5{:,data_A_idx+1})]',...
        'E',[std(sample1{:,data_A_idx+1});std(sample3{:,data_A_idx+1});std(sample4{:,data_A_idx+1});std(sample5{:,data_A_idx+1})]',...
        'BarFaceColor',permute([0.2, 0.8, 0.8; .8, 0.2, 0.2; 0.7,0.3,0.3; 0.6,0.4,0.4],[3,1,2]),...
        'ErrorbarColor','k');
    xlabel('Cytokine'); ylabel('Concentration (pg/mL)'); title('Cell Lines and Primaries All Days 6-8')
    set(gca,'XTick',1:9,'XTickLabel',Cytokines(data_A_idx),'XTickLabelRotation',90);
    legend(hb(1,1:4),{'Cell Lines 1','Primary 0188','Primary 0190','Primary 0191'})
    print(gcf, '-dpsc', figuresFilename,'-append')
    
        subplot(1,3,2)
    [hb, ~, ~, ~, ~] = superbar(1:9,[mean(sample1{:,data_B_idx+1});mean(sample3{:,data_B_idx+1});mean(sample4{:,data_B_idx+1});mean(sample5{:,data_B_idx+1})]',...
        'E',[std(sample1{:,data_B_idx+1});std(sample3{:,data_B_idx+1});std(sample4{:,data_B_idx+1});std(sample5{:,data_B_idx+1})]',...
        'BarFaceColor',permute([0.2, 0.8, 0.8; .8, 0.2, 0.2; 0.7,0.3,0.3; 0.6,0.4,0.4],[3,1,2]),...
        'ErrorbarColor','k');
    xlabel('Cytokine'); ylabel('Concentration (pg/mL)'); title('Cell Lines and Primaries All Days 6-8')
    set(gca,'XTick',1:9,'XTickLabel',Cytokines(data_B_idx),'XTickLabelRotation',90);
    legend(hb(1,1:4),{'Cell Lines 1','Primary 0188','Primary 0190','Primary 0191'})
    print(gcf, '-dpsc', figuresFilename,'-append')
    
        subplot(1,3,3)
    [hb, ~, ~, ~, ~] = superbar(1:9,[mean(sample1{:,data_C_idx+1});mean(sample3{:,data_C_idx+1});mean(sample4{:,data_C_idx+1});mean(sample5{:,data_C_idx+1})]',...
        'E',[std(sample1{:,data_C_idx+1});std(sample3{:,data_C_idx+1});std(sample4{:,data_C_idx+1});std(sample5{:,data_C_idx+1})]',...
        'BarFaceColor',permute([0.2, 0.8, 0.8; .8, 0.2, 0.2; 0.7,0.3,0.3; 0.6,0.4,0.4],[3,1,2]),...
        'ErrorbarColor','k');
    xlabel('Cytokine'); ylabel('Concentration (pg/mL)'); title('Cell Lines and Primaries All Days 6-8')
    set(gca,'XTick',1:9,'XTickLabel',Cytokines(data_C_idx),'XTickLabelRotation',90);
    legend(hb(1,1:4),{'Cell Lines 1','Primary 0188','Primary 0190','Primary 0191'})
    print(gcf, '-dpsc', figuresFilename,'-append')

    % Sorted all on one axis
    figure
    [hb, ~, ~, ~, ~] = superbar(1:27,[mean(sample1{:,max_data_idx+1});mean(sample3{:,max_data_idx+1});mean(sample4{:,max_data_idx+1});mean(sample5{:,max_data_idx+1})]',...
        'E',[std(sample1{:,max_data_idx+1});std(sample3{:,max_data_idx+1});std(sample4{:,max_data_idx+1});std(sample5{:,max_data_idx+1})]',...
        'BarFaceColor',permute([0.2, 0.8, 0.8; .8, 0.2, 0.2; 0.7,0.3,0.3; 0.6,0.4,0.4],[3,1,2]),...
        'ErrorbarColor','k');
    xlabel('Cytokine'); ylabel('Concentration (pg/mL)'); title('Cell Lines and Primaries All Days 6-8')
    set(gca,'XTick',1:27,'XTickLabel',Cytokines(max_data_idx),'XTickLabelRotation',90);
    legend(hb(1,1:4),{'Cell Lines 1','Primary 0188','Primary 0190','Primary 0191'})
    print(gcf, '-dpsc', figuresFilename,'-append')
    
    % Log Scale
    figure
    logData = log10([mean(sample1{:,max_data_idx+1});mean(sample3{:,max_data_idx+1});mean(sample4{:,max_data_idx+1});mean(sample5{:,max_data_idx+1})]);
    logPosStdevs = log10([mean(sample1{:,max_data_idx+1});mean(sample3{:,max_data_idx+1});mean(sample4{:,max_data_idx+1});mean(sample5{:,max_data_idx+1})]+[std(sample1{:,max_data_idx+1});std(sample3{:,max_data_idx+1});std(sample4{:,max_data_idx+1});std(sample5{:,max_data_idx+1})])-logData;
    logNegStdevs = log10([mean(sample1{:,max_data_idx+1});mean(sample3{:,max_data_idx+1});mean(sample4{:,max_data_idx+1});mean(sample5{:,max_data_idx+1})]-[std(sample1{:,max_data_idx+1});std(sample3{:,max_data_idx+1});std(sample4{:,max_data_idx+1});std(sample5{:,max_data_idx+1})])-logData;
    logBothStdevs(:,:,1) = logNegStdevs;
    logBothStdevs(:,:,2) = logPosStdevs;
    logBothStdevs(:,:,1)=0; %Only show positive error bars
    [hb, ~, ~, ~, ~] = superbar(1:27,logData',...
        'E',logBothStdevs,...
        'BarFaceColor',permute([0.2, 0.8, 0.8; .8, 0.2, 0.2; 0.7,0.3,0.3; 0.6,0.4,0.4],[3,1,2]),...
        'ErrorbarColor','k');
     xlabel('Cytokine'); ylabel('Concentration (pg/mL)'); title('Cell Lines and Primaries All Days 6-8')
     set(gca,'XTick',1:27,'XTickLabel',Cytokines(max_data_idx),'XTickLabelRotation',90);
     legend(hb(1,1:4),{'Cell Lines 1','Primary 0188','Primary 0190','Primary 0191'})
     print(gcf, '-dpsc', figuresFilename,'-append')
     
tableS2 = table(mean(sample1{:,2:end})');
tableS2.Properties.VariableNames={'CellLines'};
tableS2.Subject0188 = mean(sample3{:,2:end})';
%tableS2.Subject0188vCellLines = bhpvals1';
tableS2.Subject0190 = mean(sample4{:,2:end})';
%tableS2.Subject0190vCellLines = bhpvals2';
tableS2.Subject0191 = mean(sample5{:,2:end})';
%tableS2.Subject0191vCellLines = bhpvals3';
tableS2.Properties.RowNames = {sample1.Properties.VariableNames{2:end}};
%Comparison among subjects
[h,p4]=ttest2(sample3{:,2:end},sample4{:,2:end},'vartype','unequal');
[h,p5]=ttest2(sample3{:,2:end},sample5{:,2:end},'vartype','unequal');
[h,p6]=ttest2(sample4{:,2:end},sample5{:,2:end},'vartype','unequal');
combinedPvals = [p1,p2,p3,p4,p5,p6];
bhCombinedPvals = mafdr(combinedPvals, 'BHFDR',1);
bhPvals = [bhCombinedPvals(1:27)', bhCombinedPvals(27+1:2*27)', bhCombinedPvals(2*27+1:3*27)', bhCombinedPvals(3*27+1:4*27)', bhCombinedPvals(4*27+1:5*27)', bhCombinedPvals(5*27+1:6*27)']
tableS2.CellLinesVs0188 = bhPvals(:,1);
tableS2.CellLinesVs0190 = bhPvals(:,2);
tableS2.CellLinesVs0191 = bhPvals(:,3);
tableS2.s0188Vs0190 = bhPvals(:,4);
tableS2.s0188Vs0191 = bhPvals(:,5);
tableS2.s0190Vs0191 = bhPvals(:,6);

%% Stim vs. Unstim Table
% Primary Day 15-16 samples, stim and unstim
    sample1temp = regexp(uniqueSampleNames,'0188 D15 2');
    sample1 = meanSampleConcTab(cellfun(@length,sample1temp(:))>0,:);
    sample2temp = regexp(uniqueSampleNames,'0188 D15 1');
    sample2 = meanSampleConcTab(cellfun(@length,sample2temp(:))>0,:);
    [~,pvals0188] = ttest2(sample1{:,2:end},sample2{:,2:end},'vartype','unequal');
    bhpvals0188 = mafdr(pvals0188,'BHFDR',1);
    
    sample3temp = regexp(uniqueSampleNames,'0190 D16 2');
    sample3 = meanSampleConcTab(cellfun(@length,sample3temp(:))>0,:);
    sample4temp = regexp(uniqueSampleNames,'0190 D16 1');
    sample4 = meanSampleConcTab(cellfun(@length,sample4temp(:))>0,:);
    [~,pvals0190] = ttest2(sample3{:,2:end},sample4{:,2:end},'vartype','unequal');
    bhpvals0190 = mafdr(pvals0190,'BHFDR',1);
    
    sample5temp = regexp(uniqueSampleNames,'0191 D15 2');
    sample5 = meanSampleConcTab(cellfun(@length,sample5temp(:))>0,:);
    sample6temp = regexp(uniqueSampleNames,'0191 D15 1');
    sample6 = meanSampleConcTab(cellfun(@length,sample6temp(:))>0,:);
    [~,pvals0191] = ttest2(sample5{:,2:end},sample6{:,2:end},'vartype','unequal');
    bhpvals0191 = mafdr(pvals0191,'BHFDR',1);
    
    tableS3 = table(mean(sample1{:,2:end},'omitnan')');
    tableS3.Properties.VariableNames = {'s0188unstim'};
    tableS3.s0188stim = mean(sample2{:,2:end},'omitnan')';
    tableS3.pVals0188 = bhpvals0188';
    tableS3.s0190unstim = mean(sample3{:,2:end},'omitnan')';
    tableS3.s0190stim = mean(sample4{:,2:end},'omitnan')';
    tableS3.pVals0190 = bhpvals0190';
    tableS3.s0191unstim = mean(sample5{:,2:end},'omitnan')';
    tableS3.s0191stim = mean(sample6{:,2:end},'omitnan')';
    tableS3.pVals0191 = bhpvals0191';
    tableS3.Properties.RowNames = {sample1.Properties.VariableNames{2:end}};
    
    
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
