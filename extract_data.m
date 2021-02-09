function [datacoord] = extract_data(data,markers,anatLandmark)
%data: input numeric data extracted from the xls file
%markers: Designations given to the anatomical markers, eg, Xyphoid, L.Ischium, etc.
%anatLandmark: Label of the landmark recorded
%RETURNS: An array of position data for the specified anatomical landmark

%New system
label = strfind(markers(3,:),anatLandmark);

% label = strfind(marker(4,:),anatLandmark);

datapos = find(~cellfun(@isempty,label),1); %#ok<STRCLFH>

datax = data(:,datapos);
datay = data(:,datapos + 1);
dataz = data(:,datapos + 2);

datacoord = [datax datay dataz];

end