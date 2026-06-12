% Convert the .png output of CellPose to cell boundary coordinates.

function ConvertCellPose2Bound_folder2(o, CellposeF, SliceNb)

MainF = o.OutputDirectory;
%% load CellPose segmentation mask
if nargin>1 && ~isempty(CellposeF) 
     CellMap  = imread(fullfile(CellposeF,[SliceNb,'_cp_masks.png']));
     Suff = '_Cellpose';
else
    CM = load(fullfile(MainF,'CellMap.mat'));
    CellMap  = CM.CellMap;
    Suff = '';
end

%% convert the mask to boundaries
NB = max(max(CellMap));
DapiBound=cell(NB,1);
[S, iSorted] = sort(CellMap(:));
markers = [ find(diff(S)); numel(CellMap)];
for h = 1:numel(markers)-1
    idces = iSorted(markers(h)+1:markers(h+1));
    if length(idces)>10
        idces = idces(round(linspace(1,length(idces),30)));
    end
    [row,col] = ind2sub(size(CellMap),idces);
    B = boundary(row ,col,0.01);
    DapiBound{h,1} = [col(B),row(B)];
end
%% plot the boundaries on the DAPI image
figure
DapiImg = imread(fullfile(MainF,'background_image_fixed.tif'));
imshow(imadjust(DapiImg))
Boundaries = DapiBound;
SliceROI = DapiBound;
for iROI = 1: length(Boundaries)
    if ~isempty(Boundaries{iROI}) && ~iscell(Boundaries{iROI})
        BB = boundary(Boundaries{iROI}(:,1),Boundaries{iROI}(:,2),0.1);
        Boundaries{iROI} = SliceROI{iROI}(BB,:);
        hold on
        plot(Boundaries{iROI}(:,1),Boundaries{iROI}(:,2))
    end
    if mod(iROI,500) == 0
        fprintf('\n ROI %d done over %d',iROI, length(Boundaries))
    elseif iROI == length(Boundaries)
        fprintf('\n All boundaries done \n')
    end
end

% Save DAPI boundaries
save(fullfile(MainF,strcat(['DAPI_Bound',Suff,'.mat'])), 'DapiBound','CellMap');
