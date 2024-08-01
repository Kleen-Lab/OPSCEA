function plot_ecog(ecog, d, sfx, nch, nns, t, scl, ts, ytl, chanorder, showlabels, S)
tile(ecog);
dtoplot=d(nns,(t-S.marg+1):(t-S.marg+1)+sfx*S.iceegwin);
tstoplot=ts((t-S.marg+1):(t-S.marg+1)+sfx*S.iceegwin);
shift = repmat(-1*(1:nch)',1,size(dtoplot,2));
plot(tstoplot,dtoplot*scl+shift,'k');
ylim([-nch-1 0])
axis tight; xlabel('time (sec)')
hold on;
fill(ones(1,4)*ts(t)+[0 S.llw S.llw 0],[.5 .5 -nch-1.25 -nch-1.25],[.4 .4 .4],'facealpha',.25,'edgealpha',1); hold off; % overlay transform window
xlabel('Time (seconds)');
textfactor=min([ceil((length(find(nns))-80)/20) 4]); %scale text size
if showlabels
    set(gca,'ytick',-length(ytl):-1,'yticklabel',flipud(ytl(chanorder)),'fontsize',8-textfactor)
else
    set(gca,'ytick',[]);
    ylabel('Channels (randomized order)');
end
set(gca,'ylim',[-(nch)-1.25 .25])
text(t/sfx-.01+S.llw*.36,2,'v');
text(repmat(t/sfx+.01705+S.llw*.4,1,4),2.5:1:5.5,{'|','|','|','|'}); %draws an arrow pointing to the transform window
ttl1=title('ICEEG'); set(ttl1,'fontsize',10)
end