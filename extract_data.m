function [datacoord] = extract_data(data,markers,anatLandmark)
% data: numeric data extracted from the xls file
% markers: Designations given to the anatomical markers in the xls file
% anatLandmark: Label of the landmark requested
% Returns position data for the requested landmark in x,y,z columns

label = strfind(markers(3,:),anatLandmark);

datapos = find(~cellfun(@isempty,label),1); %#ok<STRCLFH>

datax = data(:,datapos);
datay = data(:,datapos + 1);
dataz = data(:,datapos + 2);

datacoord = [datax datay dataz];

end