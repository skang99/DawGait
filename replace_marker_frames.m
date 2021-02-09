function [return_matrix] = replace_marker_frames(t0,t1)

empty_frames = find(~t0);
a = zeros(length(t0),3);
a(empty_frames) = t1(empty_frames);

return_matrix = a + t0;

end

