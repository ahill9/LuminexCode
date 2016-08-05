function [hb, he, hpt, hpl, hpb, bhpvals] = makeLuminexSuperbarPlot(meanSampleConcTab,grep1, grep2, name,ttestType)
    Cytokines = meanSampleConcTab.Properties.VariableNames(2:end);
    uniqueSampleNames = meanSampleConcTab{:,1};
    sample1temp = regexp(uniqueSampleNames,grep1);
    sample1 = meanSampleConcTab(cellfun(@length,sample1temp(:))>0,:);
    sample2temp = regexp(uniqueSampleNames,grep2);
    sample2 = meanSampleConcTab(cellfun(@length,sample2temp(:))>0,:);
    if strcmp(ttestType,'paired')
        [~,p]=ttest(sample1{:,2:end},sample2{:,2:end});
    else
        [~,p]=ttest2(sample1{:,2:end},sample2{:,2:end});
    end
    bhpvals = mafdr(p,'BHFDR',true);
    pvals4plot = zeros(2*length(bhpvals),2*length(bhpvals))/0;
    for n = 1:size(sample1,2)-1
        pvals4plot(n,27+n) = bhpvals(n); pvals4plot(27+n,n) = bhpvals(n);
    end
    figure
    [hb, he, hpt, hpl, hpb] = superbar(1:27,[mean(sample1{:,2:end});mean(sample2{:,2:end})]',...
        'E',[std(sample1{:,2:end});std(sample2{:,2:end})]',...
        'P',pvals4plot,'PStarShowNS',false,'PLineColor',[1,1,1],'PLineWidth',0.01,...
        'BarFaceColor',permute([0.2, 0.8, 0.8; .8, 0.2, 0.2],[3,1,2]),...
        'ErrorbarColor','k');
    xlabel('Cytokine'); ylabel('Concentration (pg/mL)'); title(name)
    set(gca,'XTick',1:27,'XTickLabel',Cytokines,'XTickLabelRotation',90);
    legend([hb(1,1),hb(1,2)],{grep1,grep2})
end