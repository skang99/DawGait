function [marker_x,marker_y,marker_z] = extract_XYZ(position_data)
%Extracts the x,y,z coordinates of a given landmark's position data.

marker_x = position_data(:,1);
marker_y = position_data(:,2);
marker_z = position_data(:,3);


end

