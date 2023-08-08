function [D,xM] = PointCloudRegisterSB(y0, x0,o,shift)

x = [x0,ones(size(x0,1),1)];

%Initialize variables
D = zeros(3,2);
D(1:2,:) = eye(2);
D(3,:) = shift;

fprintf('\nPCR - Finding well isolated points');
% find well isolated points as those whose second neighbor is far
y = y0;

% make kd tree - default options!
k0 = KDTreeSearcher(y0);
[~, d2] = k0.knnsearch(y0, 'k', 2);
if  size(y0,1) > 1
    y = y0(d2(:,2)>2,:);
end

fprintf('\nPCR - Making kd trees');
%Make kd trees out of these well isolated points
k = KDTreeSearcher(y);

%%
Neighbor = [];
IsConverged = false;

for i=1:o.PcIter
    
    LastNeighbor = Neighbor;
    
%     x(:,1:2) = (x0-D(3,:))/D(1:2,:);
    
    x_t = x;
    xM = x_t*D(:,:)+o.TileCentre;
    if IsConverged; continue; end
    Neighbor = k.knnsearch(xM);
    [~,Dist] = k.knnsearch(xM);
    UseMe = Dist<o.PcDist;
    MyNeighb = Neighbor(UseMe>0);
%     nMatches = sum(UseMe);
%     Error = sqrt(mean(Dist(UseMe>0).^2));
    D(:,:) = x_t(UseMe>0,:)\(y(MyNeighb,:)-o.TileCentre);
   IsConverged = isequal(Neighbor,LastNeighbor);
   fprintf('\nPCR - Iteration %d:',i);
   
end
% figure(29387648);
% clf; hold on
% plot(y(:,2), y(:,1), 'g+');
% plot(xM(:,2), xM(:,1), 'r+');
% plot([xM(UseMe>0,2) y(MyNeighb,2)]',...
%     [xM(UseMe>0,1) y(MyNeighb,1)]', 'k-', 'linewidth', 1);
% daspect([1 1 1])
