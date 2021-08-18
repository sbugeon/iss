function [AllSpotScore,BestGeneNo] = get_weight_gene_dot_product(o,...
z_scoredSpotColors,BledCodeWeight,AlreadyAddedGenes)
%% OLD VERSION, superseded by get_weight_gene_dot_product2

%% Normalise BledCodes
% Need to use GeneEfficiency information to preferentially weight rounds which 
% are stronger. Also need to ensure no bias to genes with larger average
% GeneEfficiency.
%BledCodeWeight is basically step function i.e. neglect rounds with low
%GeneEfficiency (<0.5) and only keep rounds with high (>0.5).
nCodes = length(o.CharCodes);
BledCodes = reshape(o.ompBledCodes(1:nCodes,:),[nCodes,o.nBP,o.nRounds]);
SqRoundNorm = zeros(nCodes,7);
for g=1:nCodes
    for r=1:o.nRounds
        SqRoundNorm(g,r) = sum(BledCodes(g,:,r).^2);
    end
end

%%Sanity Check, now FinalNorm should equal 7.
%Want ΣWeightedSqRoundNorm = o.nRounds so can compare genes fairly
%SumWeightedSqRoundNorm = sum(BledCodeWeight.*SqRoundNorm,2);
%BledCodeWeightNormFactor = o.nRounds./SumWeightedSqRoundNorm;
%BledCodeWeight = sqrt(BledCodeWeight.*BledCodeWeightNormFactor);
% NormBledCodes = BledCodes;
% for r=1:o.nRounds
%     NormBledCodes(:,:,r) = BledCodes(:,:,r).*BledCodeWeight(:,r);
% end
% SqRoundNorm = zeros(nCodes,7);
% for g=1:nCodes
%     for r=1:o.nRounds
%         SqRoundNorm(g,r) = sum(NormBledCodes(g,:,r).^2);
%     end
% end
% FinalNorm = sum(SqRoundNorm,2);  %Should equal 7 for all genes.

%% SpotColor Weightings
%To find DotProduct use weighting from current residual because wan't to
%see if current residual can be improved by further genes. Not how genes
%already present can be altered as in get_spot_resiudal_weights.

%Boost weaker rounds - need contribution from all rounds to be good match
%to gene, not dominated by one intense round.
nSpots = size(z_scoredSpotColors,2);
SpotColorsReshape = reshape(z_scoredSpotColors,[o.nBP,o.nRounds,nSpots]);
SpotWeight = zeros(o.nRounds,nSpots);
for r=1:o.nRounds
    SpotWeight(r,:) = vecnorm(squeeze(SpotColorsReshape(:,r,:)),2,1);
end
SpotWeight = 1./(SpotWeight+o.ompWeightShift).^o.ompWeightPower;
%Norm so  ΣWeightedSqRoundNorm = o.nRounds for each gene/spot combination
WeightNormFactor = zeros(nSpots,nCodes);
for s=1:nSpots
    WeightNormFactor(s,:) = sum(SpotWeight(:,s)'.^2.*BledCodeWeight.^2.*SqRoundNorm,2);
end
WeightNormFactor = sqrt(o.nRounds./WeightNormFactor);
%SpotWeightFactor(:)=1;
%SpotWeightFactor = repelem(SpotWeightFactor,o.nBP,1);
%SpotColorsWeight = z_scoredSpotColors.*SpotWeightFactor;

%% Get dot product for each gene
AllSpotScore = zeros(nSpots,nCodes);
for g=1:nCodes
    %Overall weight factor is combination of:
    %SpotWeightFactor: Need contribution from every round.
    %BledCodeWeight: For some genes, some rounds failed so ignore these. 
    WeightFactor = SpotWeight.*BledCodeWeight(g,:)';
    WeightFactor = WeightFactor.*WeightNormFactor(:,g)';
    WeightFactor = repelem(WeightFactor,o.nBP,1);
    gBledCode = repmat(BledCodes(g,:)',[1,nSpots]);
    gBledCode = gBledCode.*WeightFactor;
    SpotColorsWeight = z_scoredSpotColors.*WeightFactor;
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

end

