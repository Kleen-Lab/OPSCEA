function plot_frame(i, isfirstframe, LL, d, sfx, nch, nns, scl, ts, ytl, chanorder, showlabels, pt, em, depthch, axislim, datapath)
global tiles;
global S;
global loaf;


%subplot(1,1,1); %clears all axes, to start fresh each frame
tiledlayout(tiles.layout.rows, tiles.layout.cols);
w8s=LL(:,i);
w8s(~nns)=0; %make weights for electrodes, and set NaNs (bad channels) to zero

% plot the different sections
plot_ecog(tiles.ecog, d, sfx, nch, nns, i, scl, ts, ytl, chanorder, showlabels, S);
plot_cb(tiles.cb, S);
plot_surfaces(tiles.surface, pt, em, w8s, nns, depthch, nch, axislim, loaf, S);
plot_depths(tiles, nns, isfirstframe, em, w8s, pt, datapath, showlabels, axislim, S);
end