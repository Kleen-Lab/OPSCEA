function [d, sfx, LL, ntp, jumpto, scl, ts, nns] = load_prepare_clip(pt, szpath, ptsz, jumpto, em)
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

scl=2/diff(prctile(reshape(d,1,numel(d)),S.iceeg_scale)); % scl=1e4/diff(prctile(make1d(d),S.ecog_scale)); % Raw ECoG scaling for display
% jumpto=S.marg+round(jumpto*sfx); %Jump ahead (sec), so video starts this much farther into the file if desired. Input argument.
end