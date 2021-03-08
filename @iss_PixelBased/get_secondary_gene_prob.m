function o = get_secondary_gene_prob(o,LookupTable)
%% o = o.get_secondary_gene_prob;
%Take all spots with particular BestGene, but with different
%o.pxSpotCodeNo. Find prob after removing a fraction of primary. 
%Fraction of primary removed depends on distance to nearest peak which is
%BestGene. 
nCodes = length(o.CharCodes);
o.pxLogProbOverBackground2 = nan(size(o.pxLogProbOverBackground));
o.pxSpotScore2 = nan(size(o.pxSpotScore));

for BestGene=1:nCodes  
    fprintf('Getting LogProb for spots near gene %.0f\n',BestGene);
    %Find distance to nearest peak that is BestGesne
    Set = o.pxSpotBestGene==BestGene & o.pxSpotCodeNo~=BestGene;
    PeakSpots = o.pxSpotCodeNo==BestGene&o.pxSpotScore>0;
    tree = KDTreeSearcher(o.pxSpotGlobalYX(PeakSpots,:));
    [~,Dist] = tree.knnsearch(o.pxSpotGlobalYX(Set,:));
    DistScore = exp(-Dist.^2/o.pLogProb2DistParam^2); 
    %SpotColors = update_spot_background(o,o.pxSpotColors(Set,:,:),ones(sum(Set),1)*BestGene);
    SpotColors = double(o.pxSpotColors(Set,:));
    %Remove gene code only where gene should be present or significant
    %bleed through
    RemoveThresh = o.pBledCodes(BestGene,:).*o.UnbledCodes(BestGene,:);
    RemoveThresh = repelem(RemoveThresh(RemoveThresh>0),1,o.nRounds)*o.pLogProb2RemoveGeneThresh;
    %Remove fraction of Best gene based on how far to nearest Best gene. 
    SpotColors(:,o.pBledCodes(BestGene,:)>RemoveThresh)=...
        (1-DistScore).*SpotColors(:,o.pBledCodes(BestGene,:)>RemoveThresh);
    LogProbOverBackground = o.get_LogProbOverBackground(...
        reshape(SpotColors,[sum(Set),o.nBP,o.nRounds]),LookupTable);
    LogProbOverBackgroundSort = LogProbOverBackground;
    LogProbOverBackgroundSort(:,BestGene)=-inf;     %So don't assign to gene removed.
    [LogProbOverBackgroundSort,~] = sort(LogProbOverBackgroundSort,2,'descend');
    Ind = sub2ind(size(LogProbOverBackground),(1:sum(Set))',o.pxSpotCodeNo(Set));
    o.pxLogProbOverBackground2(Set) = LogProbOverBackground(Ind);
    o.pxSpotScore2(Set) = o.pxLogProbOverBackground2(Set)-LogProbOverBackgroundSort(:,2);
    if o.Graphics==2
        SetSpotCodeNo = o.pxSpotCodeNo(Set);
        %Plot histogram of spots with high log prob
        figure(462);
        subplot(2,4,1:4);
        NearHist = histogram(SetSpotCodeNo(o.pxSpotScore2(Set)>0&...
            o.pxLogProbOverBackground2(Set)>0),1:74);
        [ModeGenesCount,ModeGenes] = sort(NearHist.BinCounts,'descend');
        xticks(1.5:73.5);
        set(gca,'xticklabel',o.GeneNames);
        xtickangle(90);
        title(sprintf('Distribution of Spots near %s',o.GeneNames{BestGene}));
        
        subplot(2,4,5);
        GeneBledCode = reshape(o.pBledCodes(BestGene,:),[7,7]);
        imagesc(GeneBledCode);
        title(sprintf('%s: %.0f in total',o.GeneNames{BestGene}, sum(ModeGenesCount)));
        xticks(1:7);
        xlabel('Round');
        yticks(1:7);
        yticklabels(o.bpLabels);
        ylabel('Channel');
        OverlapThreshLow = 200;
        OverlapThreshHigh = 1000;
        for i=1:3
            subplot(2,4,i+5);
            NearGeneBledCode = reshape(o.pBledCodes(ModeGenes(i),:),[7,7]);
            imagesc(NearGeneBledCode);
            hold on
            for r=1:o.nRounds
                for b=1:o.nBP
                    if NearGeneBledCode(b,r)>OverlapThreshHigh && GeneBledCode(b,r)>OverlapThreshHigh
                        rectangle('Position',[r-0.5,b-0.5,1,1],'EdgeColor','r','LineWidth',2.5,'LineStyle','-');
                    elseif NearGeneBledCode(b,r)>OverlapThreshLow && GeneBledCode(b,r)>OverlapThreshLow
                        rectangle('Position',[r-0.5,b-0.5,1,1],'EdgeColor','g','LineWidth',1,'LineStyle',':');
                    end
                end
            end
            hold off
            title(sprintf('%s: %.0f spots',o.GeneNames{ModeGenes(i)},...
                ModeGenesCount(i)));
            xticks(1:7);
            xlabel('Round');
            yticks(1:7);
            yticklabels(o.bpLabels);
            ylabel('Channel');
        end
    end
end
fprintf('\n');
end

