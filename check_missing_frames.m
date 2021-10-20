function [missing_n_data] = check_missing_frames(marker_pos,n)
% Returns whether marker_pos has a consecutive n amount of frames missing
blank_frame = 0;

if(isempty(marker_pos) || all(marker_pos(:)==0))
    missing_n_data = 1;
    return
end 

for i=1:length(marker_pos(:,1))
    if(marker_pos(i) == 0)
        blank_frame = blank_frame + 1;
    else
        blank_frame = 0;
    end %if
    
    if(blank_frame >= n)
        missing_n_data = 1;
        return
    end %if    
end %for

missing_n_data = 0;

end

