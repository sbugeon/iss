function [BleedMatrix,DiagMeasure] = ...
    get_bleed_matrix_omp(o,SpotColors,SpotIsolated,ScoreThresh,nTries)
%% [BleedMatrix,DiagMeasure] = ...
% o.get_bleed_matrix_omp(o,SpotColors,SpotIsolated,ScoreThresh,nTries)
%
% Gets bleed matrix for SpotColors. Different for omp in that always only
% find one bleedmatrix combining spots from all rounds. Also, never set any
% entries of the bleed matrix to 0.
% Input
%   SpotColors: o.dpSpotColors normalised in some way to equalise channels
%   e.g. z-scoring or dividing by percentile in each channel. 
%   SpotIsolated: which spots are well isolated so used to compute bleed matrix.
%   ScoreThresh: spot round codes have to have DotProduct greater than this
%   with a bleed matrix column to contribute to BleedMatrix. 
%   nTries: current iteration for finding bleed matrix. 
% Output
%   BleedMatrix(:,b) is the expected intensity in any round, r, for a gene
%       which has o.CharCodes{g}(r)+1 = b. 
%   DiagMeasure: should equal nChans if bleed matrix diagonal. 
%%
if nargin<5 || isempty(nTries)
    nTries = 0;
end
nChans = size(o.UseChannels,2);
nRounds = size(o.UseRounds,2);
% now we cluster the intensity vectors to estimate the Bleed Matrix
BleedMatrix = nan(o.nBP,o.nBP); % (Measured, Real, Round)
m = permute(squeeze(squeeze(SpotColors(SpotIsolated,o.UseChannels,o.UseRounds))),[1 3 2]);
m = squeeze(reshape(m,[],size(m,1)*nRounds,nChans));
m = m(~any(isnan(m),2),:);
if strcmpi(o.ompBleedMatrixEigMethod,'Mean')
    [Cluster, v, s2] = ScaledKMeans(m, eye(nChans),ScoreThresh);
elseif strcmpi(o.ompBleedMatrixEigMethod,'Median')
    [Cluster, v, s2] = ScaledKMedians(m, eye(nChans),ScoreThresh);
end
for i=1:nChans
    BleedMatrix(o.UseChannels,o.UseChannels(i)) = v(i,:) * sqrt(s2(i));
end

if o.Graphics
    figure(98043715+nTries); clf
    imagesc(BleedMatrix);
    set(gca, 'xtick', 1:o.nBP);
    set(gca, 'XTickLabel', o.bpLabels);
    set(gca, 'ytick', 1:o.nBP);
    set(gca, 'yTickLabel', o.bpLabels);
    xlabel('Actual')
    ylabel('Measured');
    title('Bleed Matrix showing all bleed through');
    drawnow;
    %     subplot(2,3,6);
    %     caxis([0 1]);
    %     axis off
    %     colormap hot
    % %     colorbar
end

%Find max channel for each column to see if diagonal.  
    %In column i, max square should be in row i if diagonal
[~,BleedMatrixMaxChannel] = max(BleedMatrix);
DiagMeasure = sum(BleedMatrixMaxChannel(o.UseChannels)==o.UseChannels);

end

