function [ytl, nch, chanorder] = prepare_plotting(eleclabels, nns, d, showlabels)
global tiles;
global sliceinfo;
global loaf;

nplots = 2 + height(tiles.depth) + height(tiles.surface);
sliceinfo=[]; 
loaf.vrf=[]; 
loaf.apasrf=[]; 
loaf.normloaf=[];
sliceinfo.viewangle=zeros(nplots,3); 
sliceinfo.azel=[]; 
sliceinfo.corners=[]; 
ytl=eleclabels(nns,1);
nch=length(find(nns));

chanorder=1:size(d(nns,:),1);

% if desired, blinds user by randomizing channel order
if ~showlabels
    chanorder=randperm(size(d(nns,:),1));
end 
end