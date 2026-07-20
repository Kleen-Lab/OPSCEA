function [F,f] = plot_rotation_animation(f,tiles,sliceinfo)
%rotation animation (first few frames of movie) to help user orientation: start
%all slices from inferior view and rotate slowly to usual head-on view

numrotationframes=15;
if height(tiles.depth) > 0
    offset = 2 + height(tiles.surface);
    % Depths with no valid (non-bad) channels are skipped by plot_depths.m
    % (isempty(eNID)||isscalar(eNID)), so sliceinfo never gets populated
    % for them. Mirror that here instead of assuming every depth tile has
    % a populated sliceinfo entry.
    validdepths = false(height(tiles.depth),1);
    for dpth=1:height(tiles.depth)
        idx = dpth+offset;
        if idx <= numel(sliceinfo) && ~isempty(sliceinfo(idx).azel)
            validdepths(dpth) = true;
            sliceinfo(idx).azelorient=[linspace(sign(sliceinfo(idx).azel(1))*180,sliceinfo(idx).azel(1),numrotationframes);
                linspace(-90,sliceinfo(idx).azel(2),numrotationframes)];
        end
    end
    for rf=1:numrotationframes
        for dpth=1:height(tiles.depth)
            if ~validdepths(dpth)
                continue
            end
            idx = dpth+offset;
            depth = tiles.depth(dpth, :);
            tile(depth);
            view(sliceinfo(idx).azelorient(1,rf),sliceinfo(idx).azelorient(2,rf));
            if rf==1
                litebrain('i',.5);
            end
        end
        pause(.25);
        F(f)=capture_frame(gcf); %capture_frame prevents a writeVideo crash from inconsistent frame sizes
        f=f+1;
    end
end