function OPSCEAsurfslice(subject,orientation,elecs,weights,datapath,fs_dir,cax,CM,gsp,j,offset,isfirstframe)
    
    %     (This is a subfunction created as a part of) Omni-planar and surface
    %     casting of epileptiform activity (OPSCEA) (UC Case Number SF2020-281)
    %     jointly created by Dr. Jon Kleen, Ben Speidel, Dr. Robert Knowlton,
    %     and Dr. Edward Chang is licensed for non-commercial research use at
    %     no cost by the Regents of the University of California under CC
    %     BY-NC-SA 4.0 (https://creativecommons.org/licenses/by-nc-sa/4.0/).
    %     Please contact innovation@ucsf.edu if you are interested in using
    %     OPSCEA for commercial purposes.
    
    %     The following copyright notice and citation is to be included in any
    %     publication, material or media wherein all or a part of Licensed
    %     Material is contained, “Certain materials incorporated herein are
    %     Copyright © 2016 The Regents of the University of California
    %     (REGENTS). All Rights Reserved.
    
    %     Please cite the following paper in your publications if you have used
    %     our software in your research, as well as any relevant toolboxes used
    %     herein as appropriate (img_pipe, FreeSurfer): Kleen JK, Speidel B,
    %     Baud MO, Rao VR, Ammanuel SG, Hamilton LS, Chang EF, Knowlton RC.
    %     Accuracy of omni-planar and surface casting of epileptiform activity
    %     for intracranial seizure localization. In press at Epilepsia.”
    
    % Global variables, added to improve computation time (most need to be calculated the first frame only
    global loaf; %contains data & info about the volume and surfaces
    global sliceinfo; %contains data & info about each individual slice (for each depth)
    global S; %contains data & info about the ECoG data and add'l parameters from excel sheet
    global I; %contains elecmatrix (coordinates), weights (for heatmap), and index of channels to skip
    global salphamask;
    if nargin<3; error('surfslice requires at least 3 input arguments'); end
    
    %% this section sets up paths and filenames
    fsbin = strcat(fs_dir,'/bin/');
    subjimaging = fullfile(datapath, 'imaging', subject);
    brainpathmgz = fullfile(subjimaging,'mri', 'brain.mgz');
    aparcpathmgz = fullfile(subjimaging,'mri', 'aparc+aseg.mgz');
    
    bashpath=getenv('PATH');
    newpath=strcat([bashpath, ':',fsbin]); 
    if length(newpath)>1000; newpath='/Applications/freesurfer/bin/'; end
    
    %%  this section loads in and sets up volume and surface data
    setenv('BRAINMGZ',brainpathmgz);
    if isempty(loaf.vrf); 
        v=load_mgh(brainpathmgz); 
        apas = load_mgh(aparcpathmgz); 
        vr=imrotate3(v,90,[1 0 0],'cubic','crop');
        vr=imrotate3(vr,90,[0 1 0],'cubic','crop');
        vrf=flip(vr,2);
        vrf=flip(vrf,3);
        loaf.vrf=vrf;
    else vrf=loaf.vrf;
    end
    
    if isempty(loaf.apasrf) 
        apasr=imrotate3(apas,90,[1 0 0],'nearest','crop');
        apasr=imrotate3(apasr,90,[0 1 0],'nearest','crop');
        apasrf = flip(apasr,2);
        apasrf = flip(apasrf,3);
        loaf.apasrf=apasrf; clear apasrf
    end
    
    alphamask = ones(size(loaf.apasrf)); % *note: consider making alpha mask also clip points outside of axislim boundaries
    alphamask(loaf.apasrf == 15|loaf.apasrf == 46|loaf.apasrf ==7|loaf.apasrf ==16|loaf.apasrf == 8|loaf.apasrf == 47|loaf.apasrf == 0)=0;
    
    hold on;
    if isfirstframe
        ZZ = meshgrid(-128:128)'; %coordinates for spatial calculations and plotting
    
        %slope: difference of 1st and last electrodes' y coordinates divided by
        %          difference of their x coordinates (O/A), to get oblique angle
        if size(elecs,1)>6; e1=2; e2=size(elecs,1)-1; else e1=1; e2=size(elecs,1); end % get the 2nd-the-last from each side (unless 6 or less contacts), helps in case offset or bent shaft
        m=(elecs(e1,2)-elecs(e2,2))/(elecs(e1,1)-elecs(e2,1));              
            intrcpt=128+elecs(e1,2)-m*elecs(e1,1); %algebra: y-mx=b
    
        b=intrcpt - 128;
    
        centered_elecs = elecs - [0 b 0];
        [thetas, rhos, zs] = cart2pol(centered_elecs(:,1), centered_elecs(:,2), centered_elecs(:,3));
    
        theta = circ_mean(thetas(e1:e2));
        m = tan(theta);
    
        zslice = meshgrid(1:256)';
    
        xslice = cos(theta).*meshgrid(-127.5:127.5) + elecs(e1,1) + 128.5;
        yslice = sin(theta).*meshgrid(-127.5:127.5) + elecs(e1,2) + 128.5;
    
        %XX, YY, ZZ are in coordinate space.
        XX = meshgrid(-128:128).*cos(theta) + elecs(e1,1);
        YY =sin(theta).*meshgrid(-128:128) + elecs(e1,2); %get coordinates in space
        sliceinfo(j).viewangle(1:3)=[-YY(1,end) 128 0];
        theta=theta+[pi*loaf.isLdepth(j-offset)]; 
        if (m)>10; theta=theta+[pi*loaf.isRdepth(j-offset)]; end
        azel = circ_rad2ang(theta-pi) + 20*(-1*loaf.isRdepth(j-offset) + 1*loaf.isLdepth(j-offset));
        sliceinfo(j).XX=XX; sliceinfo(j).YY=YY; sliceinfo(j).ZZ=ZZ; 
        sliceinfo(j).xslice=xslice; sliceinfo(j).yslice=yslice; sliceinfo(j).zslice=zslice;
    
        sliceinfo(j).sl=m; 
        sliceinfo(j).slicenum=intrcpt; 
    
        %create surface meshes that are split along the slice plane and plot
        %only one segment so that the sliceplane is still visible
        sliceinfo(j).lsplit = splitbrain(loaf.lpial,orientation,b,m);
        sliceinfo(j).rsplit = splitbrain(loaf.rpial,orientation,b,m);
        if ~orientation_good(sliceinfo(j).lsplit.vert, m, b, orientation) || ~orientation_good(sliceinfo(j).rsplit.vert, m, b, orientation)
            azel = azel + 180;
        end
        sliceinfo(j).azel=[azel,0]; % head-on angle minues 20 degrees for each slice to add perspective
    end
    
    hold on;
    lbrn=ctmr_gauss_plot_edited(sliceinfo(j).lsplit,I.em(I.nns,:),I.w8s(I.nns),S.cax,0,S.cm,S.gsp); 
    rbrn=ctmr_gauss_plot_edited(sliceinfo(j).rsplit,I.em(I.nns,:),I.w8s(I.nns),S.cax,0,S.cm,S.gsp); 
    
    if isfirstframe
        s=slice(vrf,sliceinfo(j).xslice,sliceinfo(j).yslice,sliceinfo(j).zslice); 
        s.Visible = 'off';
        sliceinfo(j).CData = s.CData;
    end
    bread=double(sliceinfo(j).CData); % the slice
    
    %This method tries to keep all slices the same brightness and contrast levels
    if isempty(loaf.normloaf) 
        loaf.normloaf=prctile(reshape(vrf,1,prod(size(vrf))),99.9); %99.9th percentile, since outliers can give poor max
    end
    
    %improve contrast/brightness for plotting, %otherwise can be too dark and will obscure heatmap for some patients
    bread=round(bread/loaf.normloaf*254); 
    
    cim = toaster(sliceinfo(j).XX,sliceinfo(j).YY,sliceinfo(j).ZZ,bread,elecs,weights,[0 cax(2)],CM,gsp);
    INPUT = uint8(cim); 
    sliceimage = padarray(INPUT,[1 1],0,'post'); %padding to appropriate size
    
    %Use the same slicing parameters to create an alphamask to make slice
    %transparent outside of brainvolume
    a = slice(alphamask,sliceinfo(j).xslice,sliceinfo(j).yslice,sliceinfo(j).zslice, 'nearest');
    a.Visible = 'off';
    AA = padarray(a.CData, [1 1],0, 'post');
    AAnonan=AA; AAnonan(isnan(AA))=0; 
    SE = strel('disk',2);
    alphamap = bwareaopen(imopen(AAnonan,SE),50); 
    
    [xedge, yedge, zedge] = getEdges(alphamap, sliceinfo(j).XX, sliceinfo(j).YY, sliceinfo(j).ZZ);
    sliceinfo(j).corners=[xedge fliplr(xedge);  yedge fliplr(yedge);  zedge([1 1 2 2])];   %for oblique slice planes
    
    %create surface in coordinate space that slices brain in the
    %appropriate plane and apply color and transparency data
    surface(sliceinfo(j).XX,sliceinfo(j).YY,sliceinfo(j).ZZ,'CData',sliceimage,'EdgeColor','none','FaceColor','texturemap','FaceAlpha','texturemap','EdgeAlpha',0,'AlphaData',alphamap,'specularexponent',5);
    shading flat
    
    caxis(S.cax) % colormap(S.cm); 
    alim([0.1 1])
    set(gca,'Clipping','off')
    axis vis3d
    salphamask = alphamask;
end

function status = orientation_good(verts, m, b, orientation)
    % check that the chunk we're choosing is on the right (correct) side
    % of the line.
    % returns 1 if the orientation is good, 0 otherwise
    
    status = 0;
    centroid = mean(verts);
    
    if strcmp(orientation, 'c') && centroid(2) < m*centroid(1) + b
        % for coronal cut want to choose the part behind the plane (we're
        % looking back) so we want centroid below the line
        status = 1;
    elseif strcmp(orientation, 'a') && centroid(3) > m*centroid(1) + b
        % for axial cut we want to choose the part above the plane
        % (we're looking up) so we want centroid above the line
        status = 1;
    elseif strcmp(orientation, 's') && centroid(1) > (centroid(2) - b)/m
        % for sagittal cut we want the part to the right of the plane
        % (we're looking to the right)
        status = 1;
    elseif strcmp(orientation, 'oc') && centroid(2) < (centroid(3) - b)/m
        % for oblique coronal cut we want the part to the right of the
        % plane when viewed from the side (i.e. the posterior part)
        status = 1;
    elseif strcmp(orientation, 'c') && centroid(2) < (centroid(3) - b)/m
        % for oblique coronal cut we want the part to the right of the
        % plane when viewed from the side (i.e. the posterior part)
        status = 1;
    end
end
