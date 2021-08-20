function [AllSpotScore,BestGeneNo] = get_weight_gene_dot_product2(o,...
z_scoredSpotColors,BledCodeWeight,AlreadyAddedGenes)
%% [AllSpotScore,BestGeneNo] = get_weight_gene_dot_product2(o,...
% z_scoredSpotColors,BledCodeWeight,AlreadyAddedGenes)
% This finds a modified dot product taking into account contribution from
% each round more equally and accounting for genes which fail in some
% rounds as given by o.GeneEfficiency. 
%
% Input
%   z_scoredSpotColors: normalised spot colors [o.nBP*o.nRounds, nSpots]
%   BledCodeWeight = sqrt(1./(1+exp(-o.ompNormBledCodeScale*(o.GeneEfficiency-o.ompNormBledCodeShift)))); 
%       It accounts for genes which fail in some rounds.
%   AlreadyAddedGenes: previously added genes for each spot
%       to ensure don't add it again [nSpots, nGenesAdded]
% Output
%   AllSpotScore(s,g): modified dot product score for spot s, gene g.
%   BestGeneNo(s): Gene g such that 
%       max(abs(AllSpotScore(s,:))) = abs(AllSpotScore(s,g))
%%
% DP = ΣW_gr^2 * DP_r where max(DP_r)=1 and min(DP_r)=0
% DP_R = ΣW_gb * W_sg_Rb * g_bR * s_bR
% W_gr weighs each individual round dot product based on gene efficiency
% of gene g. It is such that max(DP) = o.nRounds.
% W_sg_Rb weighs each round of the spot code so they all have norm = 1 (or at
% least approximately) and also weights the unbled channel by a factor. 
% W_gb weighs each round of bled code so they all have norm = 1 and also
% weights the unbled channel by a factor. 
% If s = µg (for any scalar µ) then DP = o.nRounds.

% So max total dot product is o.nRounds
W_gr = BledCodeWeight./vecnorm(BledCodeWeight,2,2) * sqrt(o.nRounds);
W_gr = repelem(W_gr,1,o.nBP);
%% DotProduct for each round, DP_R weights
% DP_R = (Σµ_b * g_bR * s_bR) / [(Σµ_b * g_bR * g_bR)^0.5 * (Σs_bR * s_bR)^0.5]
% Where µ_b = o.ompNormBledCodeUnbledBoost if b is unbled channel of gene
% g, else µ_b = 1. 
% Can write this as:
% DP_R = ΣW_gb * W_sg_Rb * g_bR * s_bR
% where:
% W_gb = [µ_b / (Σµ_b * g_bR * g_bR)]^0.5
% W_sg_Rb = [µ_b / (Σs_bR * s_bR)]^0.5
% Note that the latter is only true if o.ompWeightShift=0 and o.ompWeightPower=1

% W_gb
nCodes = length(o.CharCodes);
BledCodes = o.ompBledCodes(1:nCodes,:);
BledCodes(o.UnbledCodes==1) = BledCodes(o.UnbledCodes==1)*sqrt(o.ompNormBledCodeUnbledBoost);
BledCodes = reshape(BledCodes,[nCodes,o.nBP,o.nRounds]);
W_gb = zeros(nCodes,o.nRounds);
for g=1:nCodes
    for r=1:o.nRounds
        W_gb(g,r) = sqrt(sum(BledCodes(g,:,r).^2));
    end
end
W_gb = repelem(1./W_gb,1,o.nBP);
W_gb(o.UnbledCodes==1) = W_gb(o.UnbledCodes==1)*sqrt(o.ompNormBledCodeUnbledBoost);

%W_sR
%Boost weaker rounds - need contribution from all rounds to be good match
%to gene, not dominated by one intense round.
nSpots = size(z_scoredSpotColors,2);
if o.ompNormBledCodeUnbledBoost == 1
    SpotColorsReshape = reshape(z_scoredSpotColors,[o.nBP,o.nRounds,nSpots]);
    W_sg_rb = zeros(o.nRounds,nSpots);
    for r=1:o.nRounds
        W_sg_rb(r,:) = vecnorm(squeeze(SpotColorsReshape(:,r,:)),2,1);
    end
    W_sg_rb = 1./(W_sg_rb+o.ompWeightShift).^o.ompWeightPower;
    W_sg_rb = repelem(W_sg_rb,o.nBP,1);
end

%% Get dot product for each gene
AllSpotScore = zeros(nSpots,nCodes);
for g=1:nCodes
    gBledCode = repmat(o.ompBledCodes(g,:)',[1,nSpots]);
    gBledCode = gBledCode.*W_gb(g,:)'.*W_gr(g,:)';
    
    % Get W_sg_Rb
    if o.ompNormBledCodeUnbledBoost~=1
        SpotColorForWeight = z_scoredSpotColors;
        InUnbledCode = repmat(o.UnbledCodes(g,:)',1,nSpots)==1;
        SpotColorForWeight(InUnbledCode) = SpotColorForWeight(InUnbledCode)*...
            sqrt(o.ompNormBledCodeUnbledBoost);
        SpotColorForWeight = reshape(SpotColorForWeight,[o.nBP,o.nRounds,nSpots]);
        W_sg_rb = zeros(o.nRounds,nSpots);
        for r=1:o.nRounds
            W_sg_rb(r,:) = vecnorm(squeeze(SpotColorForWeight(:,r,:)),2,1);
        end
        W_sg_rb = 1./(W_sg_rb+o.ompWeightShift);
        W_sg_rb = repelem(W_sg_rb,o.nBP,1);
        W_sg_rb(InUnbledCode) = W_sg_rb(InUnbledCode)*sqrt(o.ompNormBledCodeUnbledBoost);
        W_sg_rb = W_sg_rb.^o.ompWeightPower;
    end
    %
    SpotColorsWeight = z_scoredSpotColors.*W_sg_rb.*W_gr(g,:)';
    %If z_scoredSpotColors(:,s) is µBledCodes(g,:) then
    %AllSpotScore(s,g) = µ*o.nRounds
    AllSpotScore(:,g) = sum(gBledCode.*SpotColorsWeight);
    %Have to norm by sum(WeightFactor.^2) to get actual coefficient value
    %for BledCodes.
    %AllSpotScore(:,g) = sum(gBledCode.*SpotColorsWeight)./...
    %    sum(WeightFactor.^2);
end
if nargin>=4 && size(AlreadyAddedGenes,1)==nSpots
    %Avoid choosing same gene again
    nAddedGenes = size(AlreadyAddedGenes,2);
    SpotInd = repmat((1:nSpots)',nAddedGenes,1);
    PrevAddedInd = sub2ind(size(AllSpotScore),SpotInd,AlreadyAddedGenes(:));
    AllSpotScore(PrevAddedInd)=0;
end
[~,BestGeneNo] = max(abs(AllSpotScore),[],2);

ShowDotProductFitting = false;
if ShowDotProductFitting && nSpots==1
    %Graphically show how dot product is found for a particular round.
    gPlot=1;
    rPlot=7;
    ImageIndex = (rPlot-1)*o.nBP+1:(rPlot-1)*o.nBP+o.nBP;
    figure;
    
    subplot(5,2,1);
    x=0:0.01:round(max(o.GeneEfficiency(:)),2)+0.1;
    BledCodeWeightPlot = ...
        1./(1+exp(-o.ompNormBledCodeScale*(x-o.ompNormBledCodeShift)));
    BledCodeWeightPlot = sqrt(BledCodeWeightPlot);
    plot(x,BledCodeWeightPlot);
    hold on;
    gPlotEff = round(o.GeneEfficiency(gPlot,rPlot),2);
    gRoundWeight = BledCodeWeight(gPlot,rPlot);
    PlotRoundColor = 'r';
    scatter(gPlotEff,gRoundWeight,'x','MarkerEdgeColor',PlotRoundColor);
    NonPlotRoundColor = [0.7,0.7,0.7];
    for r=setdiff(o.UseRounds,rPlot)
        scatter(o.GeneEfficiency(gPlot,r),BledCodeWeight(gPlot,r),...
            '+','MarkerEdgeColor',NonPlotRoundColor);
    end
    text(gPlotEff,gRoundWeight-0.1,...
        [o.GeneNames{gPlot},' Round ', num2str(rPlot)],'Color',PlotRoundColor);
    xlabel('Gene Efficiency');
    ylabel('Round Weight');
    title('How Round Weight is Derived');
    set(gca,'TickLength',[0.001,0.001]);
    
    subplot(5,2,2);
    x=0:0.01:1;
    RoundWeightsNorm = vecnorm(BledCodeWeight(gPlot,:),2,2);
    NormRoundWeight = (x./RoundWeightsNorm).^2*o.nRounds;
    plot(x,NormRoundWeight);
    hold on
    gNormRoundWeight = W_gr(gPlot,rPlot*o.nBP).^2;
    scatter(gRoundWeight,gNormRoundWeight,'rx');
    for r=setdiff(o.UseRounds,rPlot)
        scatter(BledCodeWeight(gPlot,r),W_gr(gPlot,r*o.nBP).^2,...
            '+','MarkerEdgeColor',[0.7,0.7,0.7]);
    end
    xlabel('Round Weight');
    ylabel('Norm Round Weight');
    title(sprintf('Norm Round Weight = %.2f',gNormRoundWeight));
    set(gca,'TickLength',[0.0001,0.0001]);
    
    subplot(5,2,3);
    bImage = zeros(o.nBP,o.nRounds);
    bImage(:,rPlot) = z_scoredSpotColors(ImageIndex,1);
    imagesc(bImage);
    yticks(1:o.nBP);
    yticklabels([]);
    xticks(1:o.nRounds);
    xticklabels([]);
    caxis([min(z_scoredSpotColors(:)),max(z_scoredSpotColors(:))]);
    colormap(gca,bluewhitered);
    colorbar;
    title('Spot Residual');
    
    subplot(5,2,4);
    A_Image = zeros(o.nBP,o.nRounds);
    A_Image(:,rPlot) = o.ompBledCodes(gPlot,ImageIndex);
    imagesc(A_Image);
    yticks(1:o.nBP);
    yticklabels([]);
    xticks(1:o.nRounds);
    xticklabels([]);
    colorbar;
    title([o.GeneNames{gPlot},' Code']);
    
    subplot(5,2,5);
    SpotWeightImage = zeros(o.nBP,o.nRounds);
    SpotWeightImage(:,rPlot) = W_sg_rb(ImageIndex,1);
    imagesc(SpotWeightImage);
    yticks(1:o.nBP);
    yticklabels([]);
    xticks(1:o.nRounds);
    xticklabels([]);
    colorbar;
    title('Spot Weight');
    
    subplot(5,2,6);
    GeneWeightImage = zeros(o.nBP,o.nRounds);
    GeneWeightImage(:,rPlot) = W_gb(gPlot,ImageIndex);
    imagesc(GeneWeightImage);
    yticks(1:o.nBP);
    yticklabels([]);
    xticks(1:o.nRounds);
    xticklabels([]);
    colorbar;
    title([o.GeneNames{gPlot},' Weight']);
    
    subplot(5,2,7);
    bWeightImage = bImage.*SpotWeightImage;
    imagesc(bWeightImage);
    yticks(1:o.nBP);
    yticklabels([]);
    xticks(1:o.nRounds);
    xticklabels([]);
    colormap(gca,bluewhitered);
    colorbar;
    title('Weighted Spot Residual');
    
    subplot(5,2,8);
    A_WeightImage = A_Image.*GeneWeightImage;
    imagesc(A_WeightImage);
    yticks(1:o.nBP);
    yticklabels([]);
    xticks(1:o.nRounds);
    xticklabels([]);
    colorbar;
    title(['Weighted ',o.GeneNames{gPlot},' Code']);
    
    subplot(5,2,9);
    rDotProductImage = bWeightImage.*A_WeightImage;
    imagesc(rDotProductImage);
    yticks(1:o.nBP);
    yticklabels(o.bpLabels);
    ylabel('Channel');
    xticks(1:o.nRounds);
    xlabel('Round');
    colormap(gca,bluewhitered);
    colorbar;
    title(sprintf('Round %.0f Dot Product, Sum = %.3f',...
        rPlot,sum(rDotProductImage(:))));
    
    subplot(5,2,10);
    gPlotBledCode = repmat(o.ompBledCodes(gPlot,:)',[1,nSpots]);
    gPlotBledCode = gPlotBledCode.*W_gb(gPlot,:)'.*W_gr(gPlot,:)';
    gPlotSpotColorsWeight = z_scoredSpotColors.*W_sg_rb.*W_gr(gPlot,:)';
    hold on;
    rDotProduct = zeros(o.nRounds,1);
    rRoundWeight = zeros(o.nRounds,1);
    for r=1:o.nRounds
        rIndex = (r-1)*o.nBP+1:(r-1)*o.nBP+o.nBP;
        rRoundWeight(r) = W_gr(gPlot,r*o.nBP).^2;
        rDotProduct(r) = sum(gPlotBledCode(rIndex).*gPlotSpotColorsWeight(rIndex))./rRoundWeight(r);
    end
    rContribution = rRoundWeight.*rDotProduct;
    NonPlotRounds = setdiff(1:o.nRounds,rPlot);
    scatter(NonPlotRounds, rRoundWeight(NonPlotRounds),...
        'x','MarkerEdgeColor',NonPlotRoundColor,'HandleVisibility','off');
    scatter(NonPlotRounds, rDotProduct(NonPlotRounds),...
        '+','MarkerEdgeColor',NonPlotRoundColor,'HandleVisibility','off');
    scatter(NonPlotRounds, rContribution(NonPlotRounds),...
        'o','MarkerEdgeColor',NonPlotRoundColor,'HandleVisibility','off');
    scatter(rPlot, rRoundWeight(rPlot),...
        'x','MarkerEdgeColor',PlotRoundColor);
    scatter(rPlot, rDotProduct(rPlot),...
        '+','MarkerEdgeColor',PlotRoundColor);
    scatter(rPlot, rContribution(rPlot),...
        'o','MarkerEdgeColor',PlotRoundColor);
    legend('NormRoundWeight','RoundDotProduct','RoundContribution',...
        'FontSize',7,'Location','northwest');
    xticks(1:o.nRounds);
    xlabel('Round');
    xlim([0,o.nRounds+1]);
    title(sprintf('TotalDotProduct = ∑RoundContribution = %.3f',...
        sum(rContribution)));  
    sgtitle(sprintf('Procedure to fit %s with TotalDotProduct = %.3f, highlighting Round %.0f contribution of %.3f',...
        o.GeneNames{gPlot},AllSpotScore(gPlot),rPlot,rContribution(rPlot)));
end

end
