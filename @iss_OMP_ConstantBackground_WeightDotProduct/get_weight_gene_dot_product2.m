function [AllSpotScore,BestGeneNo] = get_weight_gene_dot_product2(o,...
z_scoredSpotColors,BledCodeWeight,AlreadyAddedGenes)
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
    %Overall weight factor is combination of:
    %SpotWeightFactor: Need contribution from every round.
    %BledCodeWeight: For some genes, some rounds failed so ignore these. 
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

end

