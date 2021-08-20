function SpotNo = iss_view_spot_omp3(o, FigNo, ImSz, SpotLocation, ScoreMethod, Track, SpotNum, Norm)
%% iss_view_spot_omp3(o, FigNo, ImSz, SpotLocation,ScoreMethod, Track, SpotNum, Norm)
%
% Produces coefficient iss_view_spot_omp2 image as well as spot color
% round/channel iss_view_spot image. Clicking on a gene plot in the coefficient
% image will remove that gene from the spot color plots according to the
% coefficients found. Clicking on any background gene will remove all
% background genes. 
%
% FigNo: o.plot figure number (default, current figure)
% ImSz: radius of image that is plotted for each round and channel.
%   Default value is 7 pixels.
% SpotLocation: logical,  if true, will use location of spot closest to
%   crosshair, otherwise will use actual position of crosshair. Default is false.
% Track: gives plots of residual and gene coefficients at each stage of
%   iteration for central pixel. 
% SpotNum: spot to look at is o.pfSpotGlobalYX(SpotNum,:) where pf
%   corresponds to ScoreMethod. Can also be yx location of interest.
% Norm: true to normalise spot colors by round/channel. May be different
%   for ScoreMethod = 'OMP' or other. false to see raw spot colors.
%
% You can change o.ResidualThreshParam, o.ResidualThreshMin and
% o.ResidualThreshMax, o.ompMaxGenes to produce different coefficients. 


%%
if nargin<3 || isempty(ImSz)
    ImSz = 7;
end
if ImSz>100
    warning('ImSz too large, setting to 7');
    ImSz = 7;
end

if nargin<4 || isempty(SpotLocation)
    SpotLocation = false;
end

if nargin<6 || isempty(Track)
    Track = false;
end

if nargin<5
    ScoreMethod = [];
end

if nargin<7
    SpotNum = [];
end

if nargin<8 || isempty(Norm)
    Norm = true;  %effect of removing genes shows up better with norm hence default
end
[SpotNo,S] = iss_view_spot_omp2(o, FigNo, ImSz, SpotLocation, ScoreMethod, Track, SpotNum);
% Get final residual to get reasonable clims.
NormSpotColorsFinal = (double(S.SpotColors)-o.z_scoreSHIFT)./o.z_scoreSCALE;
for g=1:S.nCodes+o.nBackground
    NormSpotColorsFinal = NormSpotColorsFinal - ...
        permute(repmat(reshape(o.ompBledCodes(g,:),[o.nBP,o.nRounds]),...
        1,1,size(S.coefs,1)),[3,1,2]).*S.coefs(:,g);
end
SpotColorsFinal = NormSpotColorsFinal.*o.z_scoreSCALE + o.z_scoreSHIFT;
S.Clim = zeros(2,o.nBP,o.nRounds);
%climMultiplier = 0.3;
climThresh = 100;
S.Clim(1,:,:) = min(min(min(S.SpotColors),min(SpotColorsFinal)),-climThresh);
S.Clim(2,:,:) = max(max(max(S.SpotColors),max(SpotColorsFinal)),climThresh);
S.FigNo = 27312;
figure(S.FigNo);
set(gcf,'Position',[164,108,1621,805]);
if Norm
    S.Clim = (S.Clim-o.z_scoreSHIFT)./o.z_scoreSCALE;
    S.Clim(1,:,:) = min(S.Clim(:));
    S.Clim(2,:,:) = max(S.Clim(:));
    plot_spot_colors_grid(o, (double(S.SpotColors)-o.z_scoreSHIFT)./o.z_scoreSCALE, S.PointCorrectedLocalYX,...
        S.ImSz, S.Dist, S.SpotCodeNo, S.Clim);
else
    plot_spot_colors_grid(o, S.SpotColors, S.PointCorrectedLocalYX, S.ImSz, S.Dist, S.SpotCodeNo, S.Clim);
end
sgtitle('Original SpotColors');
S.SelectGenes = [];
S.SpotColorsCurrent = S.SpotColors;
S.Norm = Norm;
assignin('base','issViewSpotOMP3Object',S);

end

