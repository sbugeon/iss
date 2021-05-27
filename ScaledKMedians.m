function [k, v, s2] = ScaledKMedians(x, v0, ScoreThresh)
% [k, v] = ScaledKMedians(x, v0);
%
% does a clustering that minimizes the norm of x_i - g_i * v_{k_i}, where:
% x_i is the i'th row of input matrix x (nPoints by nDims)
% g_i is a gain (not explicitly computed or returned)
% v_k is the k'th row of output v containing cluster medians (nClusters by nDims)
% k_i is the i'th entry of output k giving the cluster of each point (nPoints by 1)
% s2_k is the squared norm of v_k.
% ScoreThresh: x_i will only contribute to v_k if dot product > ScoreThresh
% 
% input v0 (nClusters by nDims) is the starting point. (Required.)
% 
% This assigns each row of x to a row of v based on dot product between
% normalised x and v. Then update v by taking median of all x assigned to
% that row of v. 

if nargin<3 || isempty(ScoreThresh)
    ScoreThresh = 0; % only keep good matches %INCREASE THIS AS A HACK IF NOT DIAGONAL
end

MinClusterSize = 10; % delete clusters with too few points
ConvergenceCriterion = 0; % if this many or less changed, terminate
ScoreDiffThresh = 1e-4;  % if spot equally likely two classes then ignore. 

% normalize starting points and original data
vNorm = bsxfun(@rdivide, v0, sqrt(sum(v0.^2,2)));
v = vNorm;
xNorm = bsxfun(@rdivide, x, sqrt(sum(x.^2,2)));

[nClusters, nDims] = size(v);
s2 = zeros(nClusters,1);


MaxIter = 100;
k = nan; % to make sure it doesn't finish on first iteration

for i=1:MaxIter
    kOld = k;
    vOld = v;
    
    score = xNorm*v'; % project each point onto each cluster. Use normalized so we can interpret score
    [SortScore,ClusterIndex] = sort(score,2,'descend');
    k = ClusterIndex(:,1);      % find best cluster for each point
    k(SortScore(:,1)<ScoreThresh | SortScore(:,1)-SortScore(:,2)<ScoreDiffThresh)=0; % unclusterable points
    
    if sum(k~=kOld)<=ConvergenceCriterion % need a better criterion!
        break;
    end

    % find median for points assigned to each cluster
    for c=1:nClusters
        MyPoints = x(k==c,:); % don't use normalized, to avoid overweighting weak points
        nMyPoints = length(MyPoints);
        if nMyPoints<MinClusterSize
            v(c,:) = 0;
            continue;
        end
        TopEvec = median(MyPoints)';
        %[TopEvec, s2(c)] = eigs(robustcov(MyPoints), 1);   %More Proper but gives too large values.
        s2(c) = vecnorm(TopEvec);
        TopEvec = TopEvec/s2(c);
        v(c,:) = TopEvec*sign(mean(TopEvec)); % make them positive
    end
    

%     figure(9046);
%     hist(k,0:max(k));
%     drawnow

    
end
%Square the norm of v2 so consistent with output of ScaledKMeans. 
s2 = s2.^2; 

