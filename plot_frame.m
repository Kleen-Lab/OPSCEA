function plot_frame(i, isfirstframe, LL, d, sfx, nch, nns, scl, ts, ytl, chanorder, showlabels, pt, em, depthch, axislim, datapath, maxbased)
global tiles;
global S;
global loaf;

dataloaded=~isempty(LL) || ~isempty(d); 

% Make sure the current figure has a layout matching the configured grid.
% OPSCEA creates it before the frame loop, but the app's preview callbacks
% open a bare figure, and nexttile would then auto-build a 'flow' layout whose
% grid size ignores tiles.layout, throwing off every tile index and span.
if isfirstframe
    fig = gcf;
    t = findobj(fig.Children, 'flat', '-isa', 'matlab.graphics.layout.TiledChartLayout');
    if isempty(t) || ~isequal(t(1).GridSize, [tiles.layout.rows tiles.layout.cols])
        clf(fig);
        tiledlayout(fig, tiles.layout.rows, tiles.layout.cols);
    end
end

%subplot(1,1,1); %clears all axes, to start fresh each frame

if dataloaded
    w8s=LL(:,i);
    w8s(~nns)=0; %make weights for electrodes, and set NaNs (bad channels) to zero
else
    w8s=zeros(nch,1);
end

% plot the different sections
plot_ecog(tiles.ecog, d, sfx, nch, nns, i, scl, ts, ytl, chanorder, showlabels, S, isfirstframe);
plot_cb(tiles.cb, S, isfirstframe);
plot_surfaces(tiles.surface, pt, em, w8s, nns, depthch, nch, axislim, loaf, S, maxbased, isfirstframe);
plot_depths(tiles, nns, isfirstframe, em, w8s, pt, datapath, showlabels, axislim, S, maxbased);
end