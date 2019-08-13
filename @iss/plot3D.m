function plot3D(o, BackgroundImageFile, Roi,name)
% o.plot(BackgroundImage, Roi)
%
% plot the results of in situ sequencing spot detection. 
% SpotYX gives coordinates; Gene is a list of strings; 
% 
% if BackgroundImage is empty or missing, loaded from file o.BigDapiFile
% and subsetted to Roi (if present)
% If it is a numerical array, that is plotted (not subsetted for ROI)
% If zero, nothing is plotted
%
% Roi = [xmin xmax ymin ymax zmin zmax] shows only this part. Whole thing
% shown if empty or missing. Must be integers, xmin and ymin must be 1
%
% sizes can be a vector or a scalar - only used for scatter, which isn't
% called anyway.
% 
% Kenneth D. Harris, 29/3/17
% GPL 3.0 https://www.gnu.org/licenses/gpl-3.0.en.html

if nargin < 4 || isempty(name)
    name = 'Gene Positions';
end

if nargin<3 || isempty(Roi)
    Roi = round([1, max(o.SpotGlobalYXZ(:,2)), ...
    1, max(o.SpotGlobalYXZ(:,1)),...
    min(o.SpotGlobalYXZ(:,3)), max(o.SpotGlobalYXZ(:,3))]);
end

if (nargin<2 || isempty(BackgroundImageFile)) && ~isempty(o.BigDapiFile)
    BackgroundImageFile = o.BigDapiFile;
elseif isempty(o.BigDapiFile)
    warning('not sure what to do with BackgroundImage, setting to off');
end

if Roi(1) ~= 1 || Roi(3) ~= 1
    warning('Set min Roi to 1');
    Roi(1) = 1;
    Roi(3) = 1;
end


%Load in Dapi image
Image3D = zeros(Roi(4),Roi(2),Roi(6)-Roi(5)+1,'uint16');
for z = Roi(5):Roi(6)
    Image3D(:,:,z-Roi(5)+1) = imread(BackgroundImageFile, z,'PixelRegion', {Roi(3:4), Roi(1:2)});
end


S.fh = figure('units','pixels','position',[500 200 800 600],'name',name,'numbertitle','off');  %Left, Bottom, Width, Height
set(gcf, 'color', 'k');
set(gca, 'color', 'k');

S.Image = Image3D;
S.Background = imagesc(S.Image(:,:,1)); hold on; colormap bone;
%set(S.Background, 'XData', [Roi(3), Roi(4)]);
%set(S.Background, 'YData', [Roi(1), Roi(2)]);
xlim([Roi(3) Roi(4)]);
ylim([Roi(1) Roi(2)]);

title(['Z Plane ' num2str(Roi(5))],'Color','w');

hold on;
set(gca, 'YDir', 'normal');
axis on

S.SpotGeneName = o.GeneNames(o.SpotCodeNo);
S.uGenes = unique(S.SpotGeneName);
% which ones pass quality threshold (combi first)
S.QualOK = o.quality_threshold;
S.SpotYXZ = o.SpotGlobalYXZ;
S.Roi = Roi;
InRoi = round(S.SpotYXZ(:,3)) == S.Roi(5) & all(S.SpotYXZ(:,1:2)>=S.Roi([1 3]) & S.SpotYXZ(:,1:2)<=S.Roi([2 4]),2);
PlotSpots = find(InRoi & S.QualOK);
[~, S.GeneNo] = ismember(S.SpotGeneName(PlotSpots), S.uGenes);
S.h = zeros(size(S.uGenes));
for i=1:length(S.uGenes)
    MySpots = PlotSpots(S.GeneNo==i);
    if any(MySpots)
        S.h(i) = plot(S.SpotYXZ(MySpots,2), S.SpotYXZ(MySpots,1), '.');
    end
end 

legend(S.h(S.h~=0), S.uGenes(S.h~=0));
legend off;

set(gca, 'Clipping', 'off');

if ~isempty(PlotSpots)
    change_gene_symbols(0);
else
    set(gcf, 'color', 'k');
    set(gcf, 'InvertHardcopy', 'off');    
end



S.sl = uicontrol('style','slide',...
                 'unit','pix',...
                 'position',[60 8 693 18],...
                 'min',1,'max',Roi(6)-Roi(5)+1,'val',1,...
                 'sliderstep',[1/(Roi(6)-Roi(5)) 1/(Roi(6)-Roi(5))],...
                 'callback',{@sl_call,S});  
set( findall( S.fh, '-property', 'Units' ), 'Units', 'Normalized' )
hold off


function [] = sl_call(varargin)
% Callback for the slider.
[h,S] = varargin{[1,3]};  % calling handle and data structure.
ZPlane = round(get(h,'value'))+S.Roi(5)-1;
S.Background = imagesc(S.Image(:,:,ZPlane-S.Roi(5)+1)); hold on; colormap bone;
%set(S.Background, 'XData', [S.Roi(3), S.Roi(4)]);
%set(S.Background, 'YData', [S.Roi(1), S.Roi(2)]);
xlim([S.Roi(3) S.Roi(4)]);
ylim([S.Roi(1) S.Roi(2)]);

hold on;
set(gca, 'YDir', 'normal');
axis on
title(['Z Plane ' num2str(ZPlane)],'Color','w');
InRoi = round(S.SpotYXZ(:,3)) == ZPlane & all(S.SpotYXZ(:,1:2)>=S.Roi([1 3]) & S.SpotYXZ(:,1:2)<=S.Roi([2 4]),2);
PlotSpots = find(InRoi & S.QualOK);
[~, S.GeneNo] = ismember(S.SpotGeneName(PlotSpots), S.uGenes);
S.h = zeros(size(S.uGenes));
for i=1:length(S.uGenes)
    MySpots = PlotSpots(S.GeneNo==i);
    if any(MySpots)
        S.h(i) = plot(S.SpotYXZ(MySpots,2), S.SpotYXZ(MySpots,1), '.');
    end
end 

legend(S.h(S.h~=0), S.uGenes(S.h~=0));
legend off;


set(gca, 'Clipping', 'off');

if ~isempty(PlotSpots)
    change_gene_symbols(0);
else
    set(gcf, 'color', 'k');
    set(gcf, 'InvertHardcopy', 'off');    
end



hold off

