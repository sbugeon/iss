% Convert the .png output of CellPose to cell boundary coordinates.

function ConvertCellPose2Bound_folder(MainF,o)
D = dir(MainF);
Dn = {D.name};
SliceL = cellfun(@(x) x(1:end-13),Dn(endsWith(Dn,'_cp_masks.png')),'UniformOutput',false);


for iSlice= 1:length(SliceL)
    SliceNb = SliceL{iSlice};
    %% load CellPose segmentation mask
%     CellMap  = imread(fullfile(MainF,[SliceNb,'_cp_masks.png']));
    CM = load(fullfile(o.OutputDirectory,'CellMap.mat'));
    CellMap  = CM.CellMap;
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
    DapiImg = imread(fullfile(MainF,[SliceNb,'.tif']));
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
    save(fullfile(MainF,strcat('DAPI_Bound_',SliceNb,'.mat')), 'DapiBound');
end