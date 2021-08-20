function FinalCoefs = get_omp_coefs2(o,z_scoredSpotColors, SpotYX)
%% coefs = get_omp_coefs2(o,z_scoredSpotColors)
% This performs omp on every pixel in z_scoredSpotColors, the stopping
% criterion is that the residual reduction between iterations falls below
% ResidualThresh or the number of genes added to the pixel exceeds
% o.ompMaxGenes. 
%
% This removes background first and then doesn't update background
%   coefficients.
% Uses a modified OMP with a different weighted dot product.
%   i.e., different method for finding gene to add but once gene known,
%   coefficient found in the same way.
%
% Intense genes found in first iteration have surrounding pixels also added.
% 
% o: iss_OMP object
% z_scoredSpotColors: spot colors that have been z-scored by channel and round.
% coefs(s,g) is the weighting of spot s for gene g found by the omp algorithm.
%   Most are zero. The last o.nBackground genes are background codes and are always non zero.

%% OMP stops when reduction in residual drops below ResidualThresh.
% Prctile bit gets 2nd largest intensity for each spot.
nCodes = length(o.CharCodes);
ResidualThresh = o.ResidualThreshParam*prctile(abs(z_scoredSpotColors(:,:))',47.5*100/49.0)';
%Vecnorm is about double prcntile value hence need to half o.ResidualThreshParam
%ResidualThresh = o.ResidualThreshParam*vecnorm(z_scoredSpotColors(:,:),2,2);  
ResidualThresh(ResidualThresh<o.ResidualThreshMin) = o.ResidualThreshMin;
ResidualThresh(ResidualThresh>o.ResidualThreshMax) = o.ResidualThreshMax;
nSpots = size(z_scoredSpotColors,1);

z_scoredSpotColors = z_scoredSpotColors(:,:)';
FinalCoefs = zeros(nSpots,nCodes+o.nBackground);
%% Remove background and then don't update BackgroundCoefs again
%Subsequent genes will be removed off new z_scoredSpotColors
%Do OMP on post background residual (Not current residual as fitting
%multiple genes all at once). 
[LastResNorm,FinalCoefs(:,nCodes+1:end),z_scoredSpotColors] = ...
    o.get_spot_residual_background(o.ompBledCodes',z_scoredSpotColors);

%% Add one gene to each pixel - do it in vectorized form.
% Weighted dot product to see which gene to add
fprintf('Adding Gene %.0f/%.0f: %.0f Spots\n',1,o.ompMaxGenes,nSpots);
%BledCodeWeight is basically step function i.e. neglect rounds with low
%GeneEfficiency (<0.5) and only keep rounds with high (>0.5).
BledCodeWeight = ...
    1./(1+exp(-o.ompNormBledCodeScale*(o.GeneEfficiency-o.ompNormBledCodeShift))); 
BledCodeWeight = sqrt(BledCodeWeight);

[~,BestGeneNo] = o.get_weight_gene_dot_product2(z_scoredSpotColors,BledCodeWeight);
ResNorm = zeros(size(LastResNorm));
CurrentResidualSpotColors = z_scoredSpotColors;  %Use to get next best gene
IterCoefs = FinalCoefs;  %Will not update all coefs
fprintf('Percentage of spot coefs found:       ');
CurrentSpotNo = 0;
for g=1:nCodes
    Use = BestGeneNo==g;
    nGeneSpots = sum(Use);    
    if nGeneSpots>0
        [ResNorm(Use),IterCoefs(Use,g),CurrentResidualSpotColors(:,Use)] = ...
            o.get_spot_residual(o.ompBledCodes',...
            z_scoredSpotColors(:,Use),g);
        CurrentSpotNo = CurrentSpotNo+nGeneSpots;
        Percent = sprintf('%.6f', round(CurrentSpotNo*100/nSpots));
        fprintf('\b\b\b\b\b%s%%',Percent(1:4));
    end
end
fprintf('\n');

%Ensure coefficient is added near intense genes
NextGeneAdd = zeros(nSpots,o.ompMaxGenes-3);
NextGeneDist = zeros(nSpots,o.ompMaxGenes-3);
gInd = ones(nSpots,1);
for g=1:nCodes
    IntenseSpotInd = find(IterCoefs(:,g)>o.ompIntenseCoefThresh);
    if length(IntenseSpotInd)<=o.ompIntenseClusterNo; continue; end
    ZeroSpotInd = find(IterCoefs(:,g)==0);
    tree = KDTreeSearcher(SpotYX(IntenseSpotInd,:));
    %Only keep pixels in cluster
    [~,NeighbDist] = tree.knnsearch(SpotYX(IntenseSpotInd,:),'K',o.ompIntenseClusterNo);
    IntenseSpotInd = IntenseSpotInd(NeighbDist(:,o.ompIntenseClusterNo)<=...
        o.ompIntenseClusterDist);
    tree = KDTreeSearcher(SpotYX(IntenseSpotInd,:));
    [~,Dist] = tree.knnsearch(SpotYX(ZeroSpotInd,:));
    NearInd = ZeroSpotInd(Dist<o.ompIntenseDistThresh);
    Ind = sub2ind(size(NextGeneAdd),NearInd,gInd(NearInd));
    NextGeneAdd(Ind) = g;
    NextGeneDist(Ind) = Dist(Dist<o.ompIntenseDistThresh);
    gInd(NearInd) = gInd(NearInd)+1;
    %If too many genes near, get rid of furthest away ones.
    MinDistVal = min(NextGeneDist,[],2);
    [~,MaxDistInd] = max(NextGeneDist,[],2);
    gInd(MinDistVal>0) = MaxDistInd(MinDistVal>0);
end

ContinuePixels = find(LastResNorm-ResNorm>ResidualThresh | NextGeneAdd(:,1)~=0);
FinalCoefs(ContinuePixels,:) = IterCoefs(ContinuePixels,:);
LastResNorm = ResNorm;

%% Add multiple genes to each pixel which exceeding residual reduction 

for i=1:o.ompMaxGenes-1
    nSpots = length(ContinuePixels);
    if nSpots==0; break; end
    fprintf('Adding Gene %.0f/%.0f: %.0f Spots\n',1+i,o.ompMaxGenes,nSpots);
    %Find next gene to add - dot product and weights from current 
    %residual to see what to add next (Not InitialSpotWeight). 
    [PrevUseGenes,~] = find(FinalCoefs(ContinuePixels,1:nCodes)'~=0);
    PrevUseGenes = reshape(PrevUseGenes,[i,nSpots])';
    [~,BestGeneNo(ContinuePixels)] = o.get_weight_gene_dot_product2(...
        CurrentResidualSpotColors(:,ContinuePixels),BledCodeWeight,PrevUseGenes);
    %Find updated coefficients of previous and new gene.
    %Updating coefficients of multiple genes so use InitialSpotWeight.
    UseGenes = zeros(nSpots,i);
    UseGenes(:,1:i) = PrevUseGenes;
    UseGenes(:,i+1) = BestGeneNo(ContinuePixels);
    if i<o.ompMaxGenes-3
        NearLargeGene = NextGeneAdd(ContinuePixels,i)~=0;
        UseGenes(NearLargeGene,i+1) =...
            NextGeneAdd(ContinuePixels(NearLargeGene),i);
    end
    fprintf('Percentage of spot coefs found:       ');
    for sIndex=1:nSpots
        s = ContinuePixels(sIndex);
        [ResNorm(s),IterCoefs(s,UseGenes(sIndex,:)),CurrentResidualSpotColors(:,s)] = ...
            o.get_spot_residual(o.ompBledCodes',z_scoredSpotColors(:,s),...
            UseGenes(sIndex,:));
        if mod(sIndex,round(nSpots/100))==0
            Percent = sprintf('%.6f', round(sIndex*100/nSpots));
            fprintf('\b\b\b\b\b%s%%',Percent(1:4));
        end
    end
    fprintf('\n');
    if i<o.ompMaxGenes-3
        ContinuePixels = find(LastResNorm-ResNorm>ResidualThresh |...
            NextGeneAdd(:,i)~=0);
        %ContinuePixels = find(LastResNorm-ResNorm>ResidualThresh);
    else
        ContinuePixels = find(LastResNorm-ResNorm>ResidualThresh);
    end
    FinalCoefs(ContinuePixels,:) = IterCoefs(ContinuePixels,:);
    LastResNorm = ResNorm;    
end

end

