function [BleedMatrix,DiagMeasure,BleedMatrixAllBleedThrough] = ...
    get_bleed_matrix(o,SpotColors,SpotIsolated,ScoreThresh,nTries)
%% [BleedMatrix,DiagMeasure] = ...
% o.get_bleed_matrix(SpotColors,SpotIsolated,ScoreThresh,nTries)
% Gets bleed matrix for SpotColors.
% Input
%   SpotColors: o.dpSpotColors normalised in some way to equalise channels
%       e.g. z-scoring or dividing by percentile in each channel. 
%   SpotIsolated: which spots are well isolated so used to compute bleed matrix.
%   ScoreThresh: spot round codes have to have DotProduct greater than this
%       with a bleed matrix column to contribute to BleedMatrix. 
%   nTries: current iteration for finding bleed matrix. 
% Output
%   BleedMatrix: BleedMatrixAllBleedThrough but with entries below a
%       threshold determined by o.BleedThroughThresh set to 0. 
%   DiagMeasure: should equal nChans if bleed matrix diagonal. 
%   BleedMatrixAllBleedThrough(:,b,r) is the expected intensity in round r for a gene
%       which has o.CharCodes{g}(r)+1 = b. 
%%
if nargin<5 || isempty(nTries)
    nTries = 0;
end
nChans = size(o.UseChannels,2);
nRounds = size(o.UseRounds,2);
% now we cluster the intensity vectors to estimate the Bleed Matrix
BleedMatrix = zeros(nChans,nChans,nRounds); % (Measured, Real, Round)
if strcmpi(o.BleedMatrixType,'Separate')
    for r=o.UseRounds
        m = squeeze(SpotColors(SpotIsolated,o.UseChannels,r)); % data: nCodes by nBases
        m = m(~any(isnan(m),2),:);
        if strcmpi(o.BleedMatrixEigMethod,'Mean')
            [Cluster, v, s2] = ScaledKMeans(m, eye(nChans),ScoreThresh);
        elseif strcmpi(o.BleedMatrixEigMethod,'Median')
            [Cluster, v, s2] = ScaledKMedians(m, eye(nChans),ScoreThresh);
        end
        for i=1:nChans
            BleedMatrix(:,i,find(o.UseRounds==r)) = v(i,:) * sqrt(s2(i));
        end
    end
    
elseif strcmpi(o.BleedMatrixType,'Single')
    m = permute(squeeze(squeeze(SpotColors(SpotIsolated,o.UseChannels,o.UseRounds))),[1 3 2]);
    m = squeeze(reshape(m,[],size(m,1)*nRounds,nChans));
    m = m(~any(isnan(m),2),:);
    if strcmpi(o.BleedMatrixEigMethod,'Mean')
        [Cluster, v, s2] = ScaledKMeans(m, eye(nChans),ScoreThresh);
    elseif strcmpi(o.BleedMatrixEigMethod,'Median')
        [Cluster, v, s2] = ScaledKMedians(m, eye(nChans),ScoreThresh);
    end
    for i=1:nChans
        BleedMatrix(:,i,1) = v(i,:) * sqrt(s2(i));
    end
    for r=2:nRounds
        BleedMatrix(:,:,r) = BleedMatrix(:,:,1);
    end
    
else
    warning('Wrong o.BleedMatrixType entry, should be either Separate or Single')
end

if o.Graphics
    figure(98043715+nTries); clf
    for i=1:nRounds
        subplot(ceil(nRounds/3),3,i);
        imagesc(BleedMatrix(:,:,i));
        title(sprintf('Round %d', o.UseRounds(i)));
        set(gca, 'xtick', 1:nChans);
        set(gca, 'XTickLabel', o.bpLabels(o.UseChannels));
        set(gca, 'ytick', 1:nChans);
        set(gca, 'yTickLabel', o.bpLabels(o.UseChannels));
        if i==4
            xlabel('Actual')
            ylabel('Measured');
        end
    end
    sgtitle('Bleed Matrices showing all bleed through');
    drawnow;
    %     subplot(2,3,6);
    %     caxis([0 1]);
    %     axis off
    %     colormap hot
    % %     colorbar
end

%Find max channel for each column to see if diagonal.  
DiagMeasure = zeros(o.nRounds,1);
for r=1:o.nRounds
    %In column i, max square should be in row i if diagonal
    [~,CurrentBleedMatrixMaxChannel] = max(BleedMatrix(:,:,r));
    DiagMeasure(r) = sum(CurrentBleedMatrixMaxChannel==1:nChans); 
end
DiagMeasure = min(DiagMeasure);
BleedMatrixAllBleedThrough = BleedMatrix;
if DiagMeasure==nChans
    %Only keep significant bleedthrough.
    for r=1:nRounds
        for b=1:nChans
            thresh = o.BleedThroughThresh*BleedMatrix(b,b,r);
            BleedMatrix(BleedMatrix(:,b,r)<thresh,b,r)=0;
        end
    end
    
    if o.Graphics
        figure(28043715); clf
        for i=1:nRounds
            subplot(ceil(nRounds/3),3,i);
            imagesc(BleedMatrix(:,:,i));
            title(sprintf('Round %d', o.UseRounds(i)));
            set(gca, 'xtick', 1:nChans);
            set(gca, 'XTickLabel', o.bpLabels(o.UseChannels));
            set(gca, 'ytick', 1:nChans);
            set(gca, 'yTickLabel', o.bpLabels(o.UseChannels));
            if i==4
                xlabel('Actual')
                ylabel('Measured');
            end
        end
        sgtitle('Final Bleed Matrices keeping only significant bleed through');
        drawnow;
    end
    
end

end

