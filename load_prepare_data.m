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