function plot_depths(tiles, nns, isfirstframe, em, w8s, pt, datapath, showlabels, axislim, S)
global I;
global sliceinfo;
depths = tiles.depth;
for d=1:height(depths)
    depth=depths(d, :);
    tile(depth);
    eN=depth.depths{1};
    [eNID,~,~]=intersect(find(nns),eN); %Get the specific channels for this depth, ignoring bad channels
    if isempty(eNID)
        axis off;
        if isfirstframe
            tiles.depth(d, :) = [];
        end
    elseif ~isempty(eNID)
        I.em=em;
        I.w8s=w8s;
        I.nns=nns;
        offset = 2 + height(tiles.surface);
        sliceinfo(d+offset).depthlabels=depth.labels;
        OPSCEAsurfslice(pt,S.sliceplane,em(eNID,:),w8s(eNID),datapath,[],S.cax,S.cm,S.gsp,d+offset,offset,isfirstframe)
        plot3(em(eNID,1),em(eNID,2)+((S.sliceplane=='c')),em(eNID,3),'k-'); % depth probe (line between electrodes)
        plot3(em(eNID,1),em(eNID,2)+((S.sliceplane=='c')),em(eNID,3),'k.','markersize',10); % depth electrodes (dots)
        cameratoolbar('setmode','')
        axis off; axis equal;
        hold off;

        %add light sources and set camera view
        delete(findall(gca,'type','light'));
        view(sliceinfo(d+offset).azel); %head-on angle for camlight reference
        camlight(0,30) %camlight at head-on angle and 25 degrees above azimuth (otherwise light reflects and will obscure cam view)

        if isfirstframe; hold on; % add color-coded planes parallel to slice, to highlight the plane of view during initial rotation below
            fill3(sliceinfo(d+offset).corners(1,:),sliceinfo(d+offset).corners(2,:)+.1,sliceinfo(d+offset).corners(3,:),depth.color,'edgecolor',depth.color,'facealpha',0,'edgealpha',.5,'linewidth',3); hold off;
        end

        for showp=find(cell2mat(tiles.surface.show)) % now add color-coded slices planes on any subplot(s) where you indicated showplanes=1
            if ~isempty(showp)
                surf = tiles.surface(showp, :);
                tile(surf);
                fill3(sliceinfo(d+offset).corners(1,:),sliceinfo(d+offset).corners(2,:),sliceinfo(d+offset).corners(3,:),depth.color,'edgecolor',depth.color,'facealpha',.1,'edgealpha',.5,'linewidth',3);
                tile(depth); %switch back to the slice subplot at hand
            end
        end
        axis(axislim);
        zoom(depth.zoom); % apply the specified zoom for that this view (usually similar for all depths but can depend on angle of slice, position of electrodes, etc)

        if showlabels 
            ttl=depth.labels; 
        else 
            ttl=['Depth ' num2str(d)];
        end
        ttloffset=-1; % vertical offset of title from bottom of axis, in millimeters
        text(mean(sliceinfo(d+offset).corners(1,:)),mean(sliceinfo(d+offset).corners(2,:)),axislim(5)+ttloffset,ttl,'HorizontalAlignment','center','VerticalAlignment','top','color',depth.color,'fontweight','bold','fontsize',14);
        colormap(gca,S.cm);
    end
end
end