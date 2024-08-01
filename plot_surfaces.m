function plot_surfaces(surfaces, pt, em, w8s, nns, depthch, nch, axislim, meshes, S)
for p=1:height(surfaces)
    surface = surfaces(p, :);
    tile(surface);
    hold off;
    srf=regexp(surface.surfaces,',','split'); % list the specific surfaces wanted for this subplot
    % srfalpha=regexp(surface.opacity,',','split'); % list their corresponding opacities (values from 0 to 1; 0=invisible, 1=opaque)
    if length(srf)~=length(surface.opacity)
        msgbox('Number of surface to plot does not match number of alpha designations, check excel sheet');
        return;
    end
    acceptedterms={'rcortex','lcortex','rhipp','lhipp','ramyg','lamyg','wholebrain', 'rcin', 'lcin', 'rins', 'lins'};
    for s=1:length(srf)
        srf{s}=lower(srf{s}); %convert to lower case for easier string matching
        if ~isempty(intersect(srf{s},acceptedterms)) %make sure user specified accepted terminologies for the meshes
            switch char(srf{s}) %see below for case "wholebrain"s
                case 'rcortex'; srfplot=meshes.Rcrtx;
                case 'lcortex'; srfplot=meshes.Lcrtx;
                case 'rhipp';   srfplot=meshes.Rhipp;
                case 'lhipp';   srfplot=meshes.Lhipp;
                case 'ramyg';   srfplot=meshes.Ramyg;
                case 'lamyg';   srfplot=meshes.Lamyg;
            end
        else
            disp(['ATTN: Row ' num2str(p + 2) ' defined as surface but does not contain an accepted mesh term']);
            disp(acceptedterms);
            error('');
        end
        % plot the individual heatmapped surface
        if exist('srfplot','var')
            hh=ctmr_gauss_plot_edited(srfplot,em(nns,:),w8s(nns),S.cax,0,S.cm,S.gsp);
            alpha(hh,surface.opacity{1}(s)); % Adjust opacity specified for that row
        else
            disp(['ALERT: One of the entries in row ' num2str(p + 2) ' is not a valid entry, accepts:']); 
            disp(acceptedterms);
        end
    end
    if isempty(intersect(srf{s},{'rcortex','lcortex'}))||strcmpi(srf,'wholebrain') %for glass brain (hipp and/or amyg only) and wholebrain plots
        glass1=ctmr_gauss_plot_edited(meshes.Rcrtx,em(nns,:),w8s(nns),S.cax,0,S.cm,S.gsp); alpha(glass1,.1);
        glass2=ctmr_gauss_plot_edited(meshes.Lcrtx,em(nns,:),w8s(nns),S.cax,0,S.cm,S.gsp); alpha(glass2,.1);
        plot3(em(depthch,1),em(depthch,2),em(depthch,3),'k.','markersize',10-5*(1/nch*10))
        if ~surface.show{1}
            plot3(em(nns,1),em(nns,2),em(nns,3),'k.','markersize',10-5*(1/nch*10));
        end
    else
        plot3(em(nns,1),em(nns,2),em(nns,3),'k.','markersize',10-5*(1/nch*10)) %plot electrodes
    end
    cameratoolbar('setmode','');
    litebrain(char(surface.view),.9);
    wb=strcmpi(srf,'wholebrain'); 
    if any(wb)
        alpha(glass1,surface.opacity{1}(wb));
        alpha(glass2,surface.opacity{1}(wb)); 
    end
    if strcmpi(surface.view,'i')
        view(90+meshes.isL*180,270); 
    end
    if surface.show{1}||(strcmp(surface.view,'i')&&~isempty(intersect(srf,'wholebrain'))) % messes up EC72
        view(180,270); 
    end %orients the "show planes" slice to a classic axial perspective
    axis(axislim); 
    if strcmpi(surface.view,'i')||strcmpi(surface.view,'s')||strcmpi(surface.view,'a')||strcmpi(surface.view,'p') 
        if ~strcmpi(pt,'NO181') 
            axis([axislim(1)*meshes.isL+10*meshes.isR axislim(2)*meshes.isR+10*meshes.isL  axislim(3:6)]); 
        end
    end
    zoom(surface.zoom);
    hold on; colormap(gca,S.cm); set(gca,'Clipping','off')
    clear srfplot
end
end