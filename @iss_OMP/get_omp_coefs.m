function coefs = get_omp_coefs(o,z_scoredSpotColors)
%% coefs = get_omp_coefs(o,z_scoredSpotColors)
% coefs(s,g) is the weighting of spot s for gene g. Most are zero.
% The last o.nBackground are background codes and are always non zero.
% o: iss_OMP object
% z_scoredSpotColors: spot colors that have been z-scored by channel and
% round.

% OMP stops when reduction in residual drops below ResidualThresh.
% Prctile bit gets 2nd largest intensity for each spot.
nCodes = length(o.CharCodes);
ResidualThresh = o.ResidualThreshParam*prctile(abs(z_scoredSpotColors(:,:))',47.5*100/49.0)';
NonZeroSpots = ResidualThresh>0;
%Vecnorm is about double prcntile value hence need to half o.ResidualThreshParam
%ResidualThresh = o.ResidualThreshParam*vecnorm(z_scoredSpotColors(:,:),2,2);  
ResidualThresh(ResidualThresh<o.ResidualThreshMin) = o.ResidualThreshMin;
ResidualThresh(ResidualThresh>o.ResidualThreshMax) = o.ResidualThreshMax;
nSpots = size(z_scoredSpotColors,1);

coefs = zeros(nSpots,nCodes+o.nBackground);
BledCodes = o.ompBledCodes;
fprintf('Percentage of spot coefs found:       ');
for s=1:nSpots
    if NonZeroSpots(s)
        coefs(s,:) = omp_free_background(BledCodes(:,:)',z_scoredSpotColors(s,:)',...
            o.ompMaxGenes,ResidualThresh(s),nCodes+1:nCodes+o.nBackground)';
    end
    if mod(s,round(nSpots/100))==0
        Percent = sprintf('%.6f', round(s*100/nSpots));
        fprintf('\b\b\b\b\b%s%%',Percent(1:4));
    end
end
fprintf('\n');
end

