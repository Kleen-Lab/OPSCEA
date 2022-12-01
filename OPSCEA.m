function OPSCEA(pt,sz,showlabels,jumpto)
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

if ~exist('showlabels','var')||isempty(showlabels); showlabels=true; end %default displays ICEEG and depth labels
if ~exist('jumpto','var')||isempty(jumpto); jumpto=0; end 

opsceapath=['/Users/rchristin/Kleen-Lab/OPSCEA/'];   %path for parameters sheet
opsceadatapath=[opsceapath 'OPSCEADATA/'];   %path for OPSCEA ICEEG and imaging data
    if ~exist(opsceadatapath,'dir'); error('Directory for your data needs to be corrected'); end
cd(opsceapath);

ptsz=[pt '_' sz]; % prefix for filenames of specific seizure
ptpath=[opsceadatapath pt '/']; % patient's folder
szpath= [ptpath ptsz '/']; % specific seizure's folder
disp(['Running ' pt ', seizure ' sz '...']);

%% Initiate global variables
  global S; % holds general parameters
 % for speed, these are filled during first frame of surfslice then re-used
  global loaf; 
  global sliceinfo; 
  global I;


%% Import parameters
% for specific seizure 
[~,prm_allPtSz]=xlsread([opsceapath 'OPSCEAparams'],'params'); 
    fields_SZ=prm_allPtSz(1,:); % header for columns of seizure parameters
    prm=prm_allPtSz(strcmp(pt,prm_allPtSz(:,1))&strcmp(sz,prm_allPtSz(:,2)),:);
    if isempty(prm); error(['ATTENTION: No entry exists for ' pt ' seizure ' sz ' in the params master sheet']); end
% Import parameters for patient's specific plot (layout of video frame)
[~,plt]=xlsread([opsceapath 'OPSCEAparams'],pt); 
    fields_PLOT=plt(1,:); plt(1,:)=[]; % header for columns of plotting parameters
    plottype=plt(:,strcmpi(fields_PLOT,'plottype')); %type of plot for each subplot (accepts: iceeg, surface, depth, or colorbar)

cd 
%% prepare subplot specifications
    subplotrow=str2double(plt(:,strcmpi(fields_PLOT,'subplotrow')));
    subplotcolumn=str2double(plt(:,strcmpi(fields_PLOT,'subplotcolumn')));
    subplotstart=plt(:,strcmpi(fields_PLOT,'subplotstart')); subplotstop=plt(:,strcmpi(fields_PLOT,'subplotstop')); 
    for j=1:length(plottype); subplotnum{j,1}=str2double(subplotstart{j}):str2double(subplotstop{j});
    end
    surfaces=plt(:,strcmpi(fields_PLOT,'surfaces'));
    surfacesopacity=plt(:,strcmpi(fields_PLOT,'surfacesopacity'));
    viewangle=lower(plt(:,strcmpi(fields_PLOT,'viewangle')));

%% parcel all individual depth labels, contact #s, and colors. If no depths, make it  =[];
  depthlabels=plt(:,strcmpi(fields_PLOT,'depthlabels'));
  isdepth=strcmpi(plottype,'depth'); depths=cell(size(isdepth));
  if any(isdepth)
    depthEfirst=plt(:,strcmpi(fields_PLOT,'depthEfirst')); depthElast=plt(:,strcmpi(fields_PLOT,'depthElast')); 
    for j=1:length(depths); depths{j}=str2double(depthEfirst{j}):str2double(depthElast{j});
    end
    depthcolor=plt(:,strcmpi(fields_PLOT,'depthcolor')); 
    for j=1:length(depthcolor); splt=regexp(depthcolor{j},',','split'); depthcolor{j}=str2double(splt); end 
    pltzoom=str2double(plt(:,strcmpi(fields_PLOT,'pltzoom')));
    pltshowplanes=str2double(plt(:,strcmpi(fields_PLOT,'showplanes')))==1; %logical index of plots in which to show slice planes
  end
        
%% Get time segments within the ICEEG file to use
    VIDstart=prm(:,strcmpi(fields_SZ,'VIDstart')); VIDstop=prm(:,strcmpi(fields_SZ,'VIDstop')); %chunk of data (seconds into ICEEG data file) to use from the whole ICEEG data clip for the video
    S.VIDperiod=[str2double(VIDstart{1}) str2double(VIDstop{1})];
    BLstart=prm(:,strcmpi(fields_SZ,'BLstart')); BLstop=prm(:,strcmpi(fields_SZ,'BLstop')); %chunk of data (seconds into ICEEG data file) to use for baseline (for z-score step)
    S.BLperiod=[str2double(BLstart{1}) str2double(BLstop{1})];

%transform, scaling, and display options
S.llw=str2double(prm{strcmp('llw',fields_SZ)}); %default linelength window (in seconds)
S.iceeg_scale=prm{strcmp('iceeg_scale',fields_SZ)}; %percentile (number >50 and <100), used here similar to gain ICEEG waveform display, usually 95
    if ischar(S.iceeg_scale); S.iceeg_scale=str2double(S.iceeg_scale); end 
S.fps=str2double(prm{strcmp('fps',fields_SZ)});             %frames per sec of ICEEG (default 15)
S.cax=str2double(regexp(prm{strcmp('cax',fields_SZ)},',','split'));         %color axis for heatmap
S.gsp=str2double(prm{strcmp('gsp',fields_SZ)}); %gaussian spreading parameter (default 10)
    params={'iceeg_scale','fps','cax','gsp'}; 
    paramsnans=isnan([(isnan(S.iceeg_scale) | S.iceeg_scale<=50 | S.iceeg_scale>=100)   S.fps   any(isnan(S.cax)) S.gsp]); 
    if any(paramsnans); error(['ATTENTION OPSCEA USER: The "' params{paramsnans} '" term(s) is/are in an incorrect format (perhaps number instead of string), check excel seizure parameter sheet']); 
    end 
  cm=prm{strcmp('cm',fields_SZ)};
  switch cm; case 'cmOPSCEAcool'; cm=cmOPSCEAcool; 
             case 'cmOPSCEAjet'; cm=cmOPSCEAjet; 
  end
S.cm=cm; %colormap to use for heatmap
S.iceegwin=str2double(prm{strcmp('iceegwin',fields_SZ)}); %how much trace-based ICEEG to view at a time in the ICEEG window
S.marg=str2double(prm{strcmp('marg',fields_SZ)}); %offset of real-time LL txform from beginning of viewing window (in sec; converts to samples below)
S.slicebright=str2double(prm{strcmp('slicebright',fields_SZ)}); if isnan(S.slicebright); S.slicebright=0; end %brighten up slices (usually 0 to 50)


% additional adjustment for display window
S.VIDperiod=[S.VIDperiod(1)-S.marg   S.VIDperiod(2)+S.iceegwin-S.marg]; 
S.fields=fields_SZ; clear fields

S.prm=prm; clear prm
S.prm_allPtSz=prm_allPtSz; clear prm_allPtSz
S=orderfields(S); %alphabetize the structure fields for ease of use/search

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S.sliceplane='c'; % calculate omni-planar slice angles with respect to coronal (c) plane

%% load ICEEG data, and the bad channels verified for that specific data
load([szpath ptsz])
load([szpath ptsz '_badch']); 
if size(d,1)>size(d,2); d=d'; end % orient to channels by samples
[nch,ntp]=size(d); f=1; 
disp(['Length of data to play for video is ' num2str(round(ntp/sfx)) 'sec'])

% error checks for selected time periods
if any([S.VIDperiod(1) S.BLperiod(1)]<0)
    error('VIDperiod is out of bounds of file (time < 0). Check both VIDstart and BLstart times and make sure the "marg" value (subtracted from VIDstart and BLstart), cannot be < 0'); 
elseif any([S.VIDperiod(2) S.BLperiod(2)]>ntp)
    error('VIDperiod is beyond the length of the file. Check VIDstop and BLstop times vs. length of actual ICEEG data'); 
end 

%% locate and load electrode file for labels and XYZ coordinates
    load([ptpath 'Imaging/Elecs/Electrodefile.mat']); 
    if ~exist('anatomy','var'); anatomy=cell(size(elecmatrix,1),4); end
    if size(anatomy,1)>size(elecmatrix,1); anatomy(size(elecmatrix,1)+1:end)=[]; end
    anat=anatomy; clear anatomy; if size(anat,2)>size(anat,1); anat=anat'; end
    if size(anat,2)==1; anat(:,2)=anat(:,1); end; 
    if ~exist('eleclabels','var'); eleclabels=anat(:,1); end
   em=elecmatrix; clear elecmatrix; emnan=isnan(mean(em,2)); badch(emnan)=1; em(emnan,:)=0; EKGorREF=strcmpi('EKG1',anat(:,1))|strcmpi('EKG2',anat(:,1))|strcmpi('EKG',anat(:,2))|strcmpi('EKGL',anat(:,2))|strcmpi('REF',anat(:,1)); anat(EKGorREF,:)=[]; em(EKGorREF,:)=[]; eleclabels(EKGorREF,:)=[]; 
   

isR=nansum(em(:,1))>0; isL=isR~=1; %handy binary indicators for laterality

%% Implement isL/isR fix suggested by @aarongeller, allows the specifying of th side for all depths of bilateral implants
  isRdepth = [];
  isLdepth = [];
  
  for i=1:length(depths)
    if ~isnan(depths{i})
        xval_highcontact = em(depths{i}(end),1);
        isRdepth(end+1) = xval_highcontact>=0;
        isLdepth(end+1) = xval_highcontact<0;
    else
        isRdepth(end+1) = nan;
        isLdepth(end+1) = nan;
    end
  end

%% load meshes you want to plot
meshpath='Imaging/Meshes/';
Rcortex=load([ptpath meshpath pt '_rh_pial.mat']); loaf.rpial=Rcortex; Rcrtx=Rcortex.cortex; clear Rcortex
Lcortex=load([ptpath meshpath pt '_lh_pial.mat']); loaf.lpial=Lcortex; Lcrtx=Lcortex.cortex; clear Lcortex

for i=1:length(surfaces); hippentry(i)=~isempty(strfind(surfaces{i},'hipp')); amygentry(i)=~isempty(strfind(surfaces{i},'amyg')); end; 
errmsg='ATTN: MISSING A MESH, need to add this mesh file to directory (or remove/omit from frame): ';
  if any(hippentry); Rhipp=[ptpath meshpath 'subcortical/rHipp_subcort.mat']; Lhipp=Rhipp; Lhipp(end-16)='l'; if exist(Rhipp,'file'); Rhipp=load(Rhipp); Rhipp=Rhipp.cortex;  Lhipp=load(Lhipp); Lhipp=Lhipp.cortex;     else; error([errmsg 'hipp']); end; end
  if any(amygentry); Ramyg=[ptpath meshpath 'subcortical/rAmgd_subcort.mat']; Lamyg=Ramyg; Lamyg(end-16)='l'; if exist(Ramyg,'file'); Ramyg=load(Ramyg); Ramyg=Ramyg.cortex;  Lamyg=load(Lamyg); Lamyg=Lamyg.cortex;     else; error([errmsg 'amyg']); end; end
drows=find(strcmp(plottype,'depth'))'; ndepths=length(drows);

depthch=[]; for i=1:length(drows); depthch=[depthch depths{drows(i)}]; end; clear i %identify all depth electrode channels

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
if size(em,1)>size(d,1); nch=size(em,1); d(nch+1:end,:)=[]; LL(nch+1:end,:)=[]; nns(nch+1:end)=[]; ns(ns>size(em,1))=[]; 
   fprintf(2, 'ALERT: Clipping off extra bad channel entries (make sure you have the right ICEEG and bad channel files loaded)\n');
end


%% ICEEG data processing and transform

% Filter out < 1 Hz (and up to nyquist) out to help decrease movement
% artifact and improve stability during ICEEG trace plotting
[b,a]=butter(2,[1 round(sfx/2-1)]/(sfx/2),'bandpass'); 
d(nns,:) = filtfilt(b,a,d(nns,:)')';

% Line-length transform
L=round(S.llw*sfx)-1; % number of samples to calculate line length
LL=nan(size(d)); 
for i=1:size(d,2)-L; LL(:,i)=sum(abs(diff(d(:,i:i+L),1,2)),2); end

%Normalize LL (to baseline period) as z-scores, "zLL"
BLstartsample=max([sfx*S.BLperiod(1) 1]); BLendsample=sfx*S.BLperiod(2); 
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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTTING TIME!
sliceinfo=[]; loaf.vrf=[]; loaf.apasrf=[]; loaf.normloaf=[]; sliceinfo.viewangle=zeros(size(plt,1),3); sliceinfo.azel=[]; sliceinfo.corners=[]; loaf.isR=isR; loaf.isL=isL; 
%% Implement isL/isR fix proposed by @aarongeller
loaf.isRdepth=isRdepth; loaf.isLdepth=isLdepth;
clear F; 
ytl=eleclabels(nns,1); 
nch=length(find(nns)); 

chanorder=1:size(d(nns,:),1); if ~showlabels; chanorder=randperm(size(d(nns,:),1)); end % if desired, blinds user by randomizing channel order
figure('color','w','Position',[1 5 1280 700]); 
frametimpoints=jumpto:S.fram:ntp-sfx*S.iceegwin; % timepoint index of each frame to be rendered
for i=frametimpoints; 
    if i==jumpto+2*S.fram; timerem_sec=toc*length(frametimpoints); disp(['Length of data is ' num2str(ntp/sfx) 'sec']); disp(datetime); 
        disp([' -- VIDEO ETA: ' num2str(floor(timerem_sec/3600)) 'h ' num2str(ceil(mod(timerem_sec,3600)/60)) 'm -- ']); 
        fprintf('Will be done at approx: '); disp(datetime( clock, 'InputFormat', 'HH:mm:ss' ) + seconds(timerem_sec))
    end; tic 
    isfirstframe = i==jumpto; 
    subplot(1,1,1); %clears all axes, to start fresh each frame
    w8s=LL(:,i); w8s(~nns)=0; %make weights for electrodes, and set NaNs (bad channels) to zero
    for j=1:size(plt,1); 
        subplot(subplotrow(j),subplotcolumn(j),subplotnum{j}); 
    switch upper(plottype{j,1})
      case 'ECOG' %plot the raw data
            dtoplot=d(nns,(i-S.marg+1):(i-S.marg+1)+sfx*S.iceegwin);
            tstoplot=ts((i-S.marg+1):(i-S.marg+1)+sfx*S.iceegwin);
            shift = repmat(-1*(1:nch)',1,size(dtoplot,2)); 
            plot(tstoplot,dtoplot*scl+shift,'k');
              ylim([-nch-1 0])
              axis tight; xlabel('time (sec)')
              hold on;
            fill(ones(1,4)*ts(i)+[0 S.llw S.llw 0],[.5 .5 -nch-1.25 -nch-1.25],[.4 .4 .4],'facealpha',.25,'edgealpha',1); hold off; % overlay transform window
              xlabel('Time (seconds)'); 
              textfactor=min([ceil((length(find(nns))-80)/20) 4]); %scale text size
              if showlabels; set(gca,'ytick',-length(ytl):-1,'yticklabel',flipud(ytl(chanorder)),'fontsize',8-textfactor); else; set(gca,'ytick',[]); ylabel('Channels (randomized order)'); end
              set(gca,'ylim',[-(nch)-1.25 .25])
              text(i/sfx-.01+S.llw*.36,2,'v'); text(repmat(i/sfx+.01705+S.llw*.4,1,4),2.5:1:5.5,{'|','|','|','|'}) %draws an arrow pointing to the transform window
              ttl1=title('ICEEG'); set(ttl1,'fontsize',10)
      case 'COLORBAR' 
                hold off; plot(1,1); title(''); axis off; cb=colorbar; cb.Ticks=[0 1]; cb.Limits=[0 1]; cb.TickLabels={'0',num2str(S.cax(2))}; cb.FontSize=11; cb.Location='west'; 
                ylabel(cb,'z-scores','fontsize',12); colormap(gca,S.cm(floor(size(S.cm,1)/2):end,:)); %Z-scores above baseline
      case 'SURFACE' % plotting surfaces only
          hold off; srf=regexp(surfaces{j},',','split'); % list the specific surfaces wanted for this subplot
          srfalpha=regexp(surfacesopacity{j},',','split'); % list their corresponding opacities (values from 0 to 1; 0=invisible, 1=opaque)
          if length(srf)~=length(srfalpha); msgbox('Number of surface to plot does not match number of alpha designations, check excel sheet'); return; end
          acceptedterms={'rcortex','lcortex','rhipp','lhipp','ramyg','lamyg','wholebrain'};
            for s=1:length(srf)
                srf{s}=lower(srf{s}); %convert to lower case for easier string matching
              if ~isempty(intersect(srf{s},acceptedterms)) %make sure user specified accepted terminologies for the meshes
                switch srf{s}; %see below for case "wholebrain"s
                    case 'rcortex'; srfplot=Rcrtx; 
                    case 'lcortex'; srfplot=Lcrtx; 
                    case 'rhipp';   srfplot=Rhipp; 
                    case 'lhipp';   srfplot=Lhipp; 
                    case 'ramyg';   srfplot=Ramyg; 
                    case 'lamyg';   srfplot=Lamyg; 
                end
              else
                  disp(['ATTN: Row ' num2str(j) ' defined as surface but does not contain an accepted mesh term']); disp(acceptedterms); error('');
              end
              % plot the individual heatmapped surface
              if exist('srfplot','var'); hh=ctmr_gauss_plot_edited(srfplot,em(nns,:),w8s(nns),S.cax,0,S.cm,S.gsp); 
                                         alpha(hh,str2double(srfalpha{s})); % Adjust opacity specified for that row
              else; disp(['ALERT: One of the entries in row ' num2str(j) ' is not a valid entry, accepts:']); disp(acceptedterms); 
              end
            end
                if isempty(intersect(srf{s},{'rcortex','lcortex'}))||strcmpi(srf,'wholebrain') %for glass brain (hipp and/or amyg only) and wholebrain plots
                    glass1=ctmr_gauss_plot_edited(Rcrtx,em(nns,:),w8s(nns),S.cax,0,S.cm,S.gsp); alpha(glass1,.1); 
                    glass2=ctmr_gauss_plot_edited(Lcrtx,em(nns,:),w8s(nns),S.cax,0,S.cm,S.gsp); alpha(glass2,.1); 
                    plot3(em(depthch,1),em(depthch,2),em(depthch,3),'k.','markersize',10-5*(1/nch*10))
                    if ~pltshowplanes(j); plot3(em(nns,1),em(nns,2),em(nns,3),'k.','markersize',10-5*(1/nch*10)); end 
                else plot3(em(nns,1),em(nns,2),em(nns,3),'k.','markersize',10-5*(1/nch*10)) %plot electrodes
                end
                cameratoolbar('setmode',''); 
                litebrain(viewangle{j},.9); 
                wb=strcmpi(srf,'wholebrain'); if any(wb); alpha(glass1,srfalpha{wb}); alpha(glass2,srfalpha{wb}); end
                if strcmp(viewangle{j},'i'); view(90+isL*180,270); end
                if pltshowplanes(j)||(strcmp(viewangle{j},'i')&&~isempty(intersect(srf,'wholebrain'))); view(180,270); end %orients the "show planes" slice to a classic axial perspective
                axis(axislim); if strcmpi(viewangle{j},'i')||strcmpi(viewangle{j},'s')||strcmpi(viewangle{j},'a')||strcmpi(viewangle{j},'p');  if ~strcmpi(pt,'NO181'); axis([axislim(1)*isL+10*isR axislim(2)*isR+10*isL  axislim(3:6)]); end; end
                zoom(pltzoom(j)); 
                hold on; colormap(gca,S.cm); set(gca,'Clipping','off')
                clear srfplot
      case 'DEPTH' %plot depth electrode with parallel slice (plus surface behind it)
        eN=depths{j}; [eNID,~,~]=intersect(find(nns),eN); %Get the specific channels for this depth, ignoring bad channels
        if      isempty(eNID); axis off; if isfirstframe; drows(drows==j)=[]; end
        elseif ~isempty(eNID); I.em=em; I.w8s=w8s; I.nns=nns; sliceinfo(j).depthlabels=depthlabels{j}; 
          OPSCEAsurfslice(pt,S.sliceplane,em(eNID,:),w8s(eNID),opsceadatapath,[],S.cax,S.cm,S.gsp,j,isfirstframe)
          plot3(em(eNID,1),em(eNID,2)+((S.sliceplane=='c')),em(eNID,3),'k-'); % depth probe (line between electrodes)
          plot3(em(eNID,1),em(eNID,2)+((S.sliceplane=='c')),em(eNID,3),'k.','markersize',10); % depth electrodes (dots)
          cameratoolbar('setmode','')
          axis off; axis equal; 
          hold off; 
          
          %add light sources and set camera view
            delete(findall(gca,'type','light')); 
            view(sliceinfo(j).azel); %head-on angle for camlight reference
            camlight(0,30) %camlight at head-on angle and 25 degrees above azimuth (otherwise light reflects and will obscure cam view)
            
            if isfirstframe; hold on; % add color-coded planes parallel to slice, to highlight the plane of view during initial rotation below
                fill3(sliceinfo(j).corners(1,:),sliceinfo(j).corners(2,:)+.1,sliceinfo(j).corners(3,:),depthcolor{j},'edgecolor',depthcolor{j},'facealpha',0,'edgealpha',.5,'linewidth',3); hold off; 
            end
            
          for showp=find(pltshowplanes) % now add color-coded slices planes on any subplot(s) where you indicated showplanes=1
              if ~isempty(showp); subplot(subplotrow(showp),subplotcolumn(showp),subplotnum{showp}); 
              fill3(sliceinfo(j).corners(1,:),sliceinfo(j).corners(2,:),sliceinfo(j).corners(3,:),depthcolor{j},'edgecolor',depthcolor{j},'facealpha',.1,'edgealpha',.5,'linewidth',3)
                    subplot(subplotrow(j),subplotcolumn(j),subplotnum{j}) %switch back to the slice subplot at hand
              end
          end
          axis(axislim); zoom(pltzoom(j)); % apply the specified zoom for that this view (usually similar for all depths but can depend on angle of slice, position of electrodes, etc)
          
          if showlabels; ttl=depthlabels{j}; else ttl=['Depth ' num2str(j-(find(isdepth==1,1)-1))]; end
          ttloffset=-1; % vertical offset of title from bottom of axis, in millimeters
          text(mean(sliceinfo(j).corners(1,:)),mean(sliceinfo(j).corners(2,:)),axislim(5)+ttloffset,ttl,'HorizontalAlignment','center','VerticalAlignment','top','color',depthcolor{j},'fontweight','bold','fontsize',14) 
          
          colormap(gca,S.cm)

        end
    end; pause(.5) 
    end
    
    %rotation (first few frames of movie) to help user orientation: start
    %all slices from inferior view and rotate slowly to usual head-on view
    if isfirstframe; numrotationframes=15; 
      if ~isempty(drows)
        for dpth=drows
            sliceinfo(dpth).azelorient=[linspace(sign(sliceinfo(dpth).azel(1))*180,sliceinfo(dpth).azel(1),numrotationframes); 
                linspace(-90,sliceinfo(dpth).azel(2),numrotationframes)]; 
        end
        for rf=1:numrotationframes; for dpth=drows; subplot(subplotrow(dpth),subplotcolumn(dpth),subplotnum{dpth}); view(sliceinfo(dpth).azelorient(1,rf),sliceinfo(dpth).azelorient(2,rf)); if rf==1; litebrain('i',.5); end; end; pause(.25)
           F(f)=getframe(gcf); f=f+1; 
        end
      end
    end
    
    % Save frame into an ongoing sequential structure
    F(f)=getframe(gcf); f=f+1; fprintf('Saved frame - '); toc
end

cd(szpath) % Save video in the same data folder for that seizure
if showlabels; vidfilename=[ptsz '_video']; else vidfilename=[num2str(str2num(pt(3:end))*11) '_' sz]; end
v=VideoWriter(vidfilename,'MPEG-4'); 
v.FrameRate = 15; 
open(v); 
writeVideo(v,F); 
close(v); 
    

