function plot_cb(cbar, S)
tile(cbar);
hold off;
plot(1,1);
title('');
axis off;
cb=colorbar;
cb.Ticks=[0 1];
cb.Limits=[0 1];
cb.TickLabels={'0',num2str(S.cax(2))};
cb.FontSize=11;
cb.Location='west';
ylabel(cb,'z-scores','fontsize',12);
colormap(gca,S.cm(floor(size(S.cm,1)/2):end,:)); %Z-scores above baseline
end