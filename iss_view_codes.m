function iss_view_codes(o, FigNo)


    if nargin>=2
        figure(FigNo);
    end
    
    %use normed SpotColors that are actually used to determine spot scores
    cSpotColors = o.cNormSpotColors;
    
    
    set(gca, 'Color', [1 1 1]*.2);
     xy = ginput(1);
     set(gca, 'color', 'k');
%    x = 5832; y = 7936;
    %If specify z plane in workspace as should be from plot3D, find spots closest to (y,x,ZPlane)
    try
        zPlane = evalin('base', 'issPlot3DZPlane');     
        [~,SpotNo] = min(sum(abs(o.SpotGlobalYXZ-[xy(2),xy(1),zPlane]),2));
    catch
        [~,SpotNo] = min(sum(abs(o.SpotGlobalYXZ(:,1:2)-[xy(2),xy(1)]),2));
    end
    CodeNo = o.SpotCodeNo(SpotNo);
    
    MeasuredCode = squeeze(cSpotColors(SpotNo,:,:));
    CodeShape = size(MeasuredCode);
    
    figure(930476530)
    subplot(2,1,1);
    imagesc(MeasuredCode); colorbar
    caxis([0 max(MeasuredCode(:))]);
    title(sprintf('Measured code: match %.3f to %s', o.SpotScore(SpotNo), o.GeneNames{CodeNo}));
    
    set(gca, 'ytick', 1:o.nBP);
    set(gca, 'YTickLabel', o.bpLabels);
    ylabel('Color Channel');
    
    subplot(2,1,2)
    cBledCode = o.BledCodes(CodeNo,:);
    imagesc(reshape(cBledCode, CodeShape)); colorbar
    caxis([0 max(cBledCode(:))]);

    title(sprintf('Predicted Code for %s, code #%d', o.GeneNames{CodeNo}, CodeNo));
    
    
    set(gca, 'ytick', 1:o.nBP);
    set(gca, 'YTickLabel', o.bpLabels);
    ylabel('Color Channel');
    xlabel('Round');

    fprintf('Spot %d at yxz=(%d,%d,%d): code %d, %s\n', ...
        SpotNo, o.SpotGlobalYXZ(SpotNo,1),o.SpotGlobalYXZ(SpotNo,2),round(o.SpotGlobalYXZ(SpotNo,3)), CodeNo, o.GeneNames{CodeNo});

    
end
    