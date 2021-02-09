function [n_data_missing] = count_missing_frames(marker_pos,n)
% Returns whether marker_pos has a consecutive n amount of frames missing
blank_frame = 0;

for i=1:length(marker_pos(:,1))
    if(marker_pos(i) == 0)
        blank_frame = blank_frame + 1;
    else
        blank_frame = 0;
    end
    
    if(blank_frame >= n)
        n_data_missing = 1;
        return
    end    
end %for

n_data_missing = 0;

end

