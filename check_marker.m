function [is_good] = check_marker(marker)
% Returns false if marker has >=10 consecutive frames missing, true otherwise

if(check_missing_frames(marker,10))
    is_good = 0;
else
    is_good = 1;
end

end

