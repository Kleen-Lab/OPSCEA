function plot_ecog(ecogtile, d, sfx, nch, nns, t, scl, ts, ytl, chanorder, showlabels, S, isfirstframe)
% AL 7/2026: reuse the line and text handles instead of 
% recreating them, and gate the static labels/limits to first frame only
global ecoginfo

tile(ecogtile); %JK NOTE 8/2024: we should replace this in all opscea with the single line of code in the function, or at least a less generic filename (e.g. tile_opscea.m)

if ~isempty(d)
    dtoplot=d(nns,(t-round(S.marg)+1):(t-round(S.marg)+1)+sfx*S.iceegwin);
    tstoplot=ts((t-round(S.marg)+1):(t-round(S.marg)+1)+sfx*S.iceegwin);
    shift = repmat(-1*(1:nch)',1,size(dtoplot,2));
    ydata = dtoplot*scl+shift;

    if isfirstframe
        hold off;
        ecoginfo.lines_h = plot(tstoplot,ydata,'k');
        ylim([-nch-1 0])
        axis tight;
        hold on;
        ecoginfo.fill_h = fill(ones(1,4)*ts(t)+[0 S.llw S.llw 0],[.5 .5 -nch-1.25 -nch-1.25],[.4 .4 .4],'facealpha',.25,'edgealpha',1); hold off; % overlay transform window
        xlabel('Time (seconds)');
        textfactor=min([ceil((length(find(nns))-80)/20) 4]); %scale text size
        if showlabels
            set(gca,'ytick',-length(ytl):-1,'yticklabel',flipud(ytl(chanorder)),'fontsize',8-textfactor)
        else
            set(gca,'ytick',[]);
            ylabel('Channels (randomized order)');
        end
        set(gca,'ylim',[-(nch)-1.25 .25])
        ecoginfo.cursor_v_h = text(t/sfx-.01+S.llw*.36,2,'v');
        ecoginfo.cursor_bars_h = text(repmat(t/sfx+.01705+S.llw*.4,1,4),2.5:1:5.5,{'|','|','|','|'}); %draws an arrow pointing to the transform window
    else
        for c=1:numel(ecoginfo.lines_h)
            set(ecoginfo.lines_h(c),'XData',tstoplot,'YData',ydata(c,:));
        end
        xlim([tstoplot(1) tstoplot(end)]); % window keeps sliding forward in time
        set(ecoginfo.fill_h,'XData',ones(1,4)*ts(t)+[0 S.llw S.llw 0]);
        set(ecoginfo.cursor_v_h,'Position',[t/sfx-.01+S.llw*.36 2 0]);
        newx = repmat(t/sfx+.01705+S.llw*.4,1,4);
        for k=1:4
            set(ecoginfo.cursor_bars_h(k),'Position',[newx(k) 2.5+(k-1) 0]);
        end
    end
end

if isfirstframe
    ttl1=title('ICEEG'); set(ttl1,'fontsize',10)
end

end