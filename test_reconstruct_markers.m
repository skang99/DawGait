function [rMarker] = reconstruct_markers(static_segments,dir,varargin)
% varargin contains positional data of the markers required to reconstruct
% rMarker: the total number of columns of the varargin cell array should be
% 3 * marker_count. static_segments contains the scalar values of the
% average length of the corresponding static marker's positional data: the
% length of this array should be 3 * (marker_count + 1)

%gets col 1 of marker 1
varargin{1}(:,1)

marker_count = length(varargin)
markers = horzcat(varargin)
rMarker_static = static_segs(length(static_segs)-2:length(static_segs))

recon_marker = [];
col = 1;

for i = 1:markers
    recon_marker(i) = markers{i}(:,col) + -(dir)*(static_segs(i) - rMarker_static(col))






end

