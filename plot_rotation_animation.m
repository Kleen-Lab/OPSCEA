function F = plot_rotation_animation(f,tiles,sliceinfo)
%rotation animation (first few frames of movie) to help user orientation: start
%all slices from inferior view and rotate slowly to usual head-on view

numrotationframes=15;
if height(tiles.depth) > 0
    offset = 2 + height(tiles.surface);
    for dpth=1:height(tiles.depth)
        sliceinfo(dpth+offset).azelorient=[linspace(sign(sliceinfo(dpth+offset).azel(1))*180,sliceinfo(dpth+offset).azel(1),numrotationframes);
            linspace(-90,sliceinfo(dpth+offset).azel(2),numrotationframes)];
    end
    for rf=1:numrotationframes
        for dpth=1:height(tiles.depth)
            depth = tiles.depth(dpth, :);
            tile(depth);
            view(sliceinfo(dpth+offset).azelorient(1,rf),sliceinfo(dpth+offset).azelorient(2,rf));
            if rf==1
                litebrain('i',.5);
            end
        end
        pause(.25);
        F(f)=getframe(gcf);
        f=f+1;
    end
end