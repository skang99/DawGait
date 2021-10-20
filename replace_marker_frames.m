function [return_matrix] = replace_marker_frames(t0,t1)
% Replaces all empty or 0 frames found in t0 with the corresponding
% frames in t1: t0 and t1 should represent the same landmark

empty_frames = find(~t0);
a = zeros(length(t0),3);
a(empty_frames) = t1(empty_frames);

return_matrix = a + t0;

end

