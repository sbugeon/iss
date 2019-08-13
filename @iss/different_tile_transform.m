function [MyPointCorrectedYXZ, error, nMatches] = different_tile_transform(o, y0, x0,CenteredMyLocalYXZ, t, t2, r, b)
% o = o.PointCloudRegister(y0, x0, A0, Options)
% 
% Using the transformation variables found by the PointCloudRegistration
% algorithm, this tries to match spots which are on tile t on round r but 
% tile t2 on the reference round. Also returns corrected coordinates of
% CenteredScaledMyLocalYXZ
%
% inputs:
% y0 is a cell containig the centered and scaled YXZ location of all spots in all rounds 
% and colour channels for all tiles
%
% x0 is a cell containing the non centered YXZ location of spots in the 
% anchor channel for all tiles
%
% CenteredMyLocalYXZ are the local coordinates of spots on tile t on round r
% and tile t2 on the reference (Anchor) round

y = y0{t,b,r};
x = [(x0{t2} - o.CentreCorrection).*[1,1,o.Zpixelsize/o.XYpixelsize],ones(size(x0{t2},1),1)];

if isempty(o.PcDist)
    o.PcDist = inf;
end

% find well isolated points as those whose second neighbor is far
k0 = KDTreeSearcher(y);
[~, d2] = k0.knnsearch(y, 'k', 2);
if isfinite(o.PcDist) && size(y,1) > 1 
    y = y(d2(:,2)>o.PcDist*2,:);
end

%Make kd trees out of these well isolated points
k = KDTreeSearcher(y);
%Not sure about this but think its correct i.e. when transforming anchor, subtract its origin and add origin of new tile.

xR = x*o.R(:,:,t2,r);        %Shift from anchor to round r
%Then shift from tile t2 to tile t
t2_to_t_shift = (o.TileOrigin(t,:,r) - o.TileOrigin(t2,:,r)).*[1,1,o.Zpixelsize/o.XYpixelsize];
Final_x = xR + t2_to_t_shift;


                
[~,Dist] = k.knnsearch(Final_x);
UseMe = Dist<o.PcDist;               
nMatches = sum(UseMe);
error = sqrt(mean(Dist(UseMe>0).^2));

CenteredMyPointCorrectedYXZ = CenteredMyLocalYXZ*o.R(:,:,t2,r) + t2_to_t_shift;
MyPointCorrectedYXZ = round(CenteredMyPointCorrectedYXZ.*[1,1,o.XYpixelsize/o.Zpixelsize] + o.CentreCorrection);

return




