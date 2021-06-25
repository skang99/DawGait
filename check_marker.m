function [is_good] = check_marker(marker)
% Returns if marker is not missing 10 or more consecutive frames

if(check_missing_frames(marker,10))
    is_good = 0;
else
    is_good = 1;
end

end

