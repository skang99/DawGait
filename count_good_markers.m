function [good_markers] = count_good_markers(varargin)
% Returns how many markers are not missing 10 or more consecutive frames

marker_count = length(varargin);
good_markers = marker_count;

for i=1:marker_count 
    if(check_marker(cell2mat(varargin(i))') == 0)
        good_markers = good_markers - 1;
    end
end %for

end

