function OPSCEA(pt,sz,showlabels,jumpto, test)
% EXAMPLE USAGE: OPSCEA('UCSF1','01',1,0)

% pt is a string such as 'UCSF4' or 'JaneDoe', acts as a prefix for files below
% sz is a string for '01' or other number of seizure for your patient, acts
% as a secondary prefix (example above becomes UCSF4_01)
% showlabels is:
% 1 if you want to show the channel labels (default)
% 0 to hide them AND randomize channels (blinding the reader to the
% electrode locations of the trace-based ICEEG as in Kleen et al. 2021)
% jumpto allows video to start at a later point, __ seconds ahead

%     Omni-planar and surface casting of epileptiform activity (OPSCEA) (UC
%     Case Number SF2020-281) jointly created by Dr. Jon Kleen, Ben
%     Speidel, Dr. Robert Knowlton, and Dr. Edward Chang is licensed for
%     non-commercial research use at no cost by the Regents of the
%     University of California under CC BY-NC-SA 4.0
%     (https://creativecommons.org/licenses/by-nc-sa/4.0/). Please contact
%     innovation@ucsf.edu if you are interested in using OPSCEA for
%     commercial purposes.

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
%     for intracranial seizure localization. Epilepsia. 2021;62(4):947-959.”
arguments
    pt;
    sz;
    showlabels logical = true;
    jumpto double = 0;
    test logical = false;
end

datapath = getenv('KLEEN_DATA');
opsceapath = fullfile(datapath, 'opscea');%path for parameters sheet
imagingpath = fullfile(datapath, 'imaging');

ptsz=[pt '_' sz]; % prefix for filenames of specific seizure
ptpath=fullfile(opsceapath, pt); % patient's folder
ptparams = fullfile(ptpath, 'patient_params.mat');
szpath=fullfile(ptpath, ptsz); % specific seizure's folder
clipparams = fullfile(szpath, 'clip_params.mat');
disp(['Running ' pt ', seizure ' sz '...']);

%% Initiate global variables
global S; % holds general parameters
% for speed, these are filled during first frame of surfslice then re-used
global loaf;
global sliceinfo;
global I;
global tiles;


%% Import parameters

load_patient_params(ptparams);
load_clip_params(clipparams, test);

%% load ICEEG data, and the bad channels verified for that specific data
[d, sfx, nns, ntp, eleclabels, LL, jumpto, scl, ts, em, depthch, axislim] = load_prepare_data(pt, szpath, ptsz, imagingpath, jumpto);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTTING TIME!
nplots = 2 + height(tiles.depth) + height(tiles.surface);
sliceinfo=[]; loaf.vrf=[]; loaf.apasrf=[]; loaf.normloaf=[];
sliceinfo.viewangle=zeros(nplots,3); sliceinfo.azel=[]; sliceinfo.corners=[]; 
%loaf.isR=isR; loaf.isL=isL;
%% Implement isL/isR fix proposed by @aarongeller
%loaf.isRdepth=isRdepth; loaf.isLdepth=isLdepth;
clear F;
ytl=eleclabels(nns,1);
nch=length(find(nns));

chanorder=1:size(d(nns,:),1); if ~showlabels; chanorder=randperm(size(d(nns,:),1)); end % if desired, blinds user by randomizing channel order
figure('color','w','Position',[1 5 1280 700]);
frametimpoints=jumpto:S.fram:ntp-sfx*S.iceegwin; % timepoint index of each frame to be rendered
f = 1;
for i=frametimpoints
    if i==jumpto+2*S.fram
        timerem_sec=toc*length(frametimpoints); 
        disp(['Length of data is ' num2str(ntp/sfx) 'sec']); disp(datetime);
        disp([' -- VIDEO ETA: ' num2str(floor(timerem_sec/3600)) 'h ' num2str(ceil(mod(timerem_sec,3600)/60)) 'm -- ']);
        fprintf('Will be done at approx: '); disp(datetime( clock, 'InputFormat', 'HH:mm:ss' ) + seconds(timerem_sec))
    end; tic
    isfirstframe = i==jumpto;
    %subplot(1,1,1); %clears all axes, to start fresh each frame
    tiledlayout(tiles.layout.rows, tiles.layout.cols);
    w8s=LL(:,i); 
    w8s(~nns)=0; %make weights for electrodes, and set NaNs (bad channels) to zero

    % plot the different sections
    plot_ecog(tiles.ecog, d, sfx, nch, nns, i, scl, ts, ytl, chanorder, showlabels, S);
    plot_cb(tiles.cb, S);
    plot_surfaces(tiles.surface, pt, em, w8s, nns, depthch, nch, axislim, loaf, S);
    plot_depths(tiles, nns, isfirstframe, em, w8s, pt, datapath, showlabels, axislim, S);

    %rotation (first few frames of movie) to help user orientation: start
    %all slices from inferior view and rotate slowly to usual head-on view
    if isfirstframe 
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
    end
    % Save frame into an ongoing sequential structure
    F(f)=getframe(gcf); 
    f=f+1; 
    fprintf('Saved frame - '); 
    toc
end

%cd(szpath) % Save video in the same data folder for that seizure
%if showlabels; vidfilename=[ptsz '_video']; else vidfilename=[num2str(str2num(pt(3:end))*11) '_' sz]; end

if test
    viddir = fullfile(datapath, 'ictal_cinema_library', pt, 'test');
else
    viddir = fullfile(datapath, 'ictal_cinema_library',pt);
end
vidfn = [pt '_' sz];
vidfilename = fullfile(viddir, vidfn);
mkdir(viddir);
v=VideoWriter(vidfilename,'MPEG-4');
v.FrameRate = 15;
open(v);
writeVideo(v,F);
close(v);
end

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

function tile(s)
nexttile(s.tile, [s.rows s.cols]);
end

function load_patient_params(ptparams)
global tiles;

load(ptparams, 'layout', 'ecog', 'cb', 'surface', 'depth');
tiles.layout = layout; clear layout;
tiles.ecog = ecog; clear ecog;
tiles.cb = cb; clear cb;
tiles.surface = surface; clear surface;
tiles.depth = depth; clear depth;
end

function load_clip_params(clipparams, test)
global S;

load(clipparams, 'vidstart', 'vidstop', 'llw', 'iceeg_scale', 'fps', 'cax', 'gsp', 'cm', 'iceegwin', 'marg', 'slicebright');
S.VIDperiod=[vidstart vidstop];
load(clipparams, 'blstart', 'blstop')
S.BLperiod=[blstart blstop];

%transform, scaling, and display options
S.llw=llw;
S.iceeg_scale=iceeg_scale; %percentile (number >50 and <100), used here similar to gain ICEEG waveform display, usually 95
S.fps=fps;
if test
    S.fps = 0.25;
end
S.cax=cax; %color axis for heatmap
S.gsp=gsp; %gaussian spreading parameter (default 10)
params={'iceeg_scale','fps','cax','gsp'};
paramsnans=isnan([(isnan(S.iceeg_scale) | S.iceeg_scale<=50 | S.iceeg_scale>=100)   S.fps   any(isnan(S.cax)) S.gsp]);
if any(paramsnans)
    error(['ATTENTION OPSCEA USER: The "' params{paramsnans} '" term(s) is/are in an incorrect format (perhaps number instead of string), check excel seizure parameter sheet']);
end
switch cm
    case 'cmOPSCEAcool'
        cm=cmOPSCEAcool;
    case 'cmOPSCEAjet'
        cm=cmOPSCEAjet;
end
S.cm=cm; %colormap to use for heatmap
S.iceegwin=iceegwin; %how much trace-based ICEEG to view at a time in the ICEEG window
S.marg=marg; %offset of real-time LL txform from beginning of viewing window (in sec; converts to samples below)
S.slicebright=slicebright;
if isnan(S.slicebright); S.slicebright=0; end %brighten up slices (usually 0 to 50)

% additional adjustment for display window
S.VIDperiod=[S.VIDperiod(1)-S.marg   S.VIDperiod(2)+S.iceegwin-S.marg];

S.sliceplane='c'; % calculate omni-planar slice angles with respect to coronal (c) plane

S=orderfields(S); %alphabetize the structure fields for ease of use/search
clear llw iceeg_scale fps cax gsp cm iceegwin marg slicebright;
end

function [d, sfx, nns, ntp, eleclabels, LL, jumpto, scl, ts, em, depthch, axislim] = load_prepare_data(pt, szpath, ptsz, imagingpath, jumpto)
global S;
global tiles;
global loaf;

%% load ICEEG data, and the bad channels verified for that specific data
if ~exist(szpath, 'dir')
    load_sz_data(pt);
end

load(fullfile(szpath, ptsz), 'd', 'sfx');
load(fullfile(szpath, [ptsz '_badch']), 'badch');

if size(d,1)>size(d,2); d=d'; end % orient to channels by samples
[nch,ntp]=size(d);
disp(['Total length of ICEEG data: ' num2str(round(ntp/sfx)) ' sec'])
disp(['Length of data to play for video: ' num2str(diff(S.VIDperiod)) ' sec'])

% error checks for selected time periods
if any([S.VIDperiod(1) S.BLperiod(1)]<0)
    error('VIDperiod is out of bounds of file (time < 0). Check both VIDstart and BLstart times and make sure the "marg" value (subtracted from VIDstart and BLstart), cannot be < 0');
elseif any([S.VIDperiod(2) S.BLperiod(2)]>ntp)
    error('VIDperiod is beyond the length of the file. Check VIDstop and BLstop times vs. length of actual ICEEG data');
end

%% locate and load electrode file for labels and XYZ coordinates
load(fullfile(imagingpath, pt, 'elecs', 'clinical_elecs_all.mat'), 'anatomy', 'elecmatrix', 'eleclabels');
if ~exist('anatomy','var')
    anatomy=cell(size(elecmatrix,1),4);
end
if size(anatomy,1)>size(elecmatrix,1)
    anatomy(size(elecmatrix,1)+1:end)=[];
end
anat=anatomy; clear anatomy;
if size(anat,2)>size(anat,1)
    anat=anat';
end
if size(anat,2)==1
    anat(:,2)=anat(:,1);
end
if ~exist('eleclabels','var')
    eleclabels=anat(:,1);
end
em=elecmatrix; clear elecmatrix;
emnan=isnan(mean(em,2));
badch(emnan)=1;
em(emnan,:)=0;
EKGorREF=strcmpi('EKG1',anat(:,1))|strcmpi('EKG2',anat(:,1))|strcmpi('EKG',anat(:,2))|strcmpi('EKGL',anat(:,2))|strcmpi('REF',anat(:,1));
anat(EKGorREF,:)=[];
em(EKGorREF,:)=[];
eleclabels(EKGorREF,:)=[];


loaf.isR=nansum(em(:,1))>0; 
loaf.isL=loaf.isR~=1; %handy binary indicators for laterality

%% Implement isL/isR fix suggested by @aarongeller, allows the specifying of th side for all depths of bilateral implants
isRdepth = [];
isLdepth = [];

for i=1:height(tiles.depth)
    if ~isnan(tiles.depth.depths{i})
        xval_highcontact = em(tiles.depth.depths{i}(end),1);
        isRdepth(end+1) = xval_highcontact>=0;
        isLdepth(end+1) = xval_highcontact<0;
    else
        isRdepth(end+1) = nan;
        isLdepth(end+1) = nan;
    end
end
loaf.isRdepth = isRdepth;
loaf.isLdepth = isLdepth;

%% load meshes you want to plot
meshpath=fullfile(imagingpath, pt, 'Meshes');
Rcortex=load(fullfile(meshpath, [pt '_rh_pial.mat'])); 
loaf.rpial=Rcortex; 
Rcrtx=Rcortex.cortex; 
loaf.Rcrtx = Rcrtx;
clear Rcortex
Lcortex=load(fullfile(meshpath, [pt '_lh_pial.mat'])); 
loaf.lpial=Lcortex; 
Lcrtx=Lcortex.cortex; 
loaf.Lcrtx = Lcrtx;
clear Lcortex

for i=1:height(tiles.surface)
    hippentry(i)=~isempty(strfind(tiles.surface.surfaces,'hipp'));
    amygentry(i)=~isempty(strfind(tiles.surface.surfaces,'amyg')); 
end
errmsg='ATTN: MISSING A MESH, need to add this mesh file to directory (or remove/omit from frame): ';
if any(hippentry)
    Rhipp=fullfile(meshpath, 'subcortical', 'rHipp_subcort.mat'); 
    Lhipp=fullfile(meshpath, 'subcortical', 'lHipp_subcort.mat');
    if exist(Rhipp,'file')
        Rhipp=load(Rhipp); 
        Rhipp=Rhipp.cortex;
        loaf.Rhipp = Rhipp;
        Lhipp=load(Lhipp); 
        Lhipp=Lhipp.cortex; 
        loaf.Lhipp = Lhipp;
    else
        error([errmsg 'hipp']); 
    end
end
if any(amygentry)
    Ramyg=fullfile(meshpath, 'subcortical', 'rAmgd_subcort.mat'); 
    Lamyg=fullfile(meshpath, 'subcortical', 'lAmgd_subcort.mat');

    if exist(Ramyg,'file')
        Ramyg=load(Ramyg); 
        Ramyg=Ramyg.cortex;
        loaf.Ramyg = Ramyg;
        Lamyg=load(Lamyg); 
        Lamyg=Lamyg.cortex;
        loaf.Lamyg = Lamyg;
    else
        error([errmsg 'amyg']); 
    end
end

depthch=[]; 
for i=1:height(tiles.depth)
    depthch=[depthch tiles.depth.depths{i}]; 
end; clear i %identify all depth electrode channels

%% get xyz limits for plotting purposes
perim=1; % how many millimeters away from brain/electrodes boundaries to set the colorcoded plane perimeter, recommend >0 to avoid skimming brain surface (default 1mm)
axl(:,:,1)=[min([Rcrtx.vert; Lcrtx.vert]); max([Rcrtx.vert; Lcrtx.vert])]; %min and max of MESH VERTICES' x y z coordinates (2x3)
axl(:,:,2)=[min(em); max(em)]; %%min and max of ELECTRODES' x y z coordinates (2x3) and stack them (2x3x2)
axl=[min(axl(1,:,:),[],3); max(axl(2,:,:),[],3)]; % Get the minima and maxima of both
axislim=reshape(axl,1,6)+[-1 1 -1 1 -1 1]*perim; clear axl %Use the, to define the axis boundaries, and add additional perimeter (perim)

%% formatting checks, and consolidation of bad channels
ns=unique( [find(badch);   find(isnan(mean(em,2)));   find(isnan(mean(d,2)))]  ); % bad channels: those that are pre-marked, or if NaNs in their coordinates or ICEEG data traces
nns=true(nch,1); nns(ns)=0; %nns=find(nns); %consolidate bad channels and those with NaNs
%remove data channels previously clipped at end. Only include that which has electrode coordinates (intracranial)
if size(em,1)<size(d,1)
    nch=size(em,1); 
    d(nch+1:end,:)=[]; 
    nns(nch+1:end)=[]; 
    ns(ns>size(em,1))=[];
    fprintf(2, 'ALERT: Clipping off extra bad channel entries (make sure you have the right ICEEG and bad channel files loaded)\n');
end

%% ICEEG data processing and transform

% Filter out < 1 Hz (and up to nyquist) out to help decrease movement
% artifact and improve stability during ICEEG trace plotting
[b,a]=butter(2,[1 round(sfx/2-1)]/(sfx/2),'bandpass'); %% is this needed since done in preprocessing?
d(nns,:) = filtfilt(b,a,d(nns,:)')';

% Line-length transform
L=round(S.llw*sfx)-1; % number of samples to calculate line length
LL=nan(size(d));
for i=1:size(d,2)-L; LL(:,i)=sum(abs(diff(d(:,i:i+L),1,2)),2); end

%Normalize LL (to baseline period) as z-scores, "zLL"
BLstartsample=max([sfx*S.BLperiod(1) 1]); 
BLendsample=sfx*S.BLperiod(2);
for i=1:nch % z-score using channel-specific baseline
    LL(i,:)=(LL(i,:)-nanmean(LL(i,BLstartsample:BLendsample)))/nanstd(LL(i,BLstartsample:BLendsample));
end

% make a timtestamp vector
ts=0:1/sfx:size(d,2)*(1/sfx)-1/sfx;

% Extract the period of data to be used for the video (remove flanking data)
vidperiodidx=round(S.VIDperiod(1)*sfx+1):S.VIDperiod(2)*sfx;
d=d(:,vidperiodidx); ntp=length(vidperiodidx);
LL=LL(:,vidperiodidx);
ts=ts(vidperiodidx);

% Scaling ICEEG and setting windows for simultaneous trace-based display
S.iceeg_scale=S.iceeg_scale+(100-S.iceeg_scale)/2; S.iceeg_scale=[100-S.iceeg_scale S.iceeg_scale]; %conversion to two-tailed percentile
scl=2/diff(prctile(reshape(d,1,numel(d)),S.iceeg_scale));
S.marg=round(S.marg*sfx); %offset of real-time LL txform from beginning of viewing window
jumpto=S.marg+round(jumpto*sfx); %Jump ahead (sec), so video starts this much farther into the file if desired. Input argument.
S.fram=round(sfx/S.fps);
end