function objNew = switch_class(objOld,DesiredClass,OverrideError)
%% objNew = switch_class(objOld,DesiredClass,OverrideError);
% To convert between classes e.g. o = switch_class(o, iss_PixelBased) will
% change the class of o to iss_PixelBased.
% objOld: iss object
% DesiredClass: iss_Base, iss_GroundTruth, iss_OMP, iss_OMP,
%   iss_OMP_ConstantBackground_WeightDotProduct, iss_PixelBased or iss_Spatial.
% OverrideError: if switch is such that data will be lost e.g. from iss_PixelBased
%   to iss_Base, an error will stop this from happening with OverrideError = false
%   (default). With OverrideError = true, you can carry out the switch.
%
% objNew: iss object of DesiredClass.

if nargin<3 || isempty(OverrideError)
    OverrideError = false;
end
NewClassProp = metaclass(DesiredClass).PropertyList;
if ~OverrideError
    %If will lose data as a result of switching class, raise error.
    OldClassProp = metaclass(objOld).PropertyList;
    NewPropNames = cell(size(NewClassProp));
    for i=1:size(NewClassProp,1)
        NewPropNames{i} = NewClassProp(i).Name;
    end
    OldPropNames = cell(size(OldClassProp));
    for i=1:size(OldClassProp,1)
        OldPropNames{i} = OldClassProp(i).Name;
    end
    DiffProp = setdiff(OldPropNames,NewPropNames);
    if size(DiffProp,1)~=0
        fprintf('Lost Data:\n');
        for i=1:size(DiffProp,1)
            fprintf([DiffProp{i},'\n']);
        end
        error(['Changing class will lose the above data. To overide this error, '...
            'run o=switch_class(o,NEW_CLASS,true);']);
    end
end

objNew = DesiredClass;
for k = 1:length(NewClassProp)
    try
        objNew.(NewClassProp(k).Name) = objOld.(NewClassProp(k).Name);
    catch
        ...
    end
end
end

