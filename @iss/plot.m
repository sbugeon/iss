function plot(o)
% plot
%
% plot the results of in situ sequencing spot detection. 
% SpotYX gives coordinates; Gene is a list of strings; 
% 
% BackgroundImage is loaded using options o
%
% sizes can be a vector or a scalar - only used for scatter, which isn't
% called anyway.
% 
% Kenneth D. Harris, 29/3/17
% GPL 3.0 https://www.gnu.org/licenses/gpl-3.0.en.html

BackgroundImageName = 'Background_DAPI';

Scatter=0; % faster
%% plot background image
% first see if it is already in figure to save time loading
cc = get(gca, 'Children');
for i=1:length(cc)
    if strcmp(cc(i).UserData, BackgroundImageName)
        fprintf('reusing background image\n')
        BackgroundImage = cc(i).CData;
        break;
    end
end

if ~exist('BackgroundImage', 'var') && ~isempty(o.BigDapiFile)
    fprintf('loading background image...');
    BackgroundImage = imread(o.BigDapiFile);
    fprintf('done\n');
end

if exist('BackgroundImage', 'var')
    clf; set(gcf, 'color', 'k');
    set(gca, 'color', 'k');
    h = imagesc(BackgroundImage);
    colormap bone;
    set(h, 'UserData', BackgroundImageName);
end

hold on;
set(gca, 'YDir', 'normal');
axis on

%% now plot spots

SpotGeneName = o.GeneNames(o.SpotCodeNo);
uGenes = unique(SpotGeneName);

% which ones pass quality threshold (combi first)
nCombiCodes = sum(~strcmp(o.CharCodes, 'EXTRA'));
PlotMe = (o.SpotCombi & o.SpotScore>o.CombiQualThresh & o.SpotIntensity>o.CombiIntensityThresh);
% now extras - they have their own thresholds, set manually
for i=1:size(o.ExtraCodes,1)
    MySpots = (o.SpotCodeNo == nCombiCodes+i);
    PlotMe(MySpots) = o.SpotIntensity(MySpots)>o.ExtraCodes{i,4};
end

if Scatter
    for i=1:length(uGenes)
        g = uGenes(i);
        MySpots = PlotMe & strcmp(SpotGeneName, g);
        h(i) = scatter(o.SpotGlobalYX(MySpots,2), o.SpotGlobalYX(MySpots,1), Sizes(MySpots));
    end
else
    for i=1:length(uGenes)
        g = uGenes(i);
        MySpots = PlotMe & strcmp(SpotGeneName, g);
        h(i) = plot(o.SpotGlobalYX(MySpots,2), o.SpotGlobalYX(MySpots,1), '.');
    end 
end

legend(uGenes);

% gscatter(SpotsYX(:,2), SpotsYX(:,1), Genes, [], [], Sizes)
change_gene_symbols(0);

end

