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
[ytl, nch, chanorder] = prepare_plotting(eleclabels, nns, d, showlabels);

figure('color','w','Position',[1 5 1280 700]);
frametimpoints=jumpto:S.fram:ntp-sfx*S.iceegwin; % timepoint index of each frame to be rendered

clear F;
f = 1;
for i=frametimpoints
    if i==jumpto+2*S.fram
        timerem_sec=toc*length(frametimpoints);
        disp(['Length of data is ' num2str(ntp/sfx) 'sec']); disp(datetime);
        disp([' -- VIDEO ETA: ' num2str(floor(timerem_sec/3600)) 'h ' num2str(ceil(mod(timerem_sec,3600)/60)) 'm -- ']);
        fprintf('Will be done at approx: '); disp(datetime( clock, 'InputFormat', 'HH:mm:ss' ) + seconds(timerem_sec))
    end; tic
    isfirstframe = i==jumpto;

    % Save each frame into the ongoing sequential structure F
    plot_frame(i, isfirstframe, LL, d, sfx, nch, nns, scl, ts, ytl, chanorder, showlabels, pt, em, depthch, axislim, datapath);
    
    if isfirstframe
        F = plot_rotation_animation(f, tiles, sliceinfo);
    end

    F(f) = getframe(gcf); 
    f = f+1; 
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
vidfn = [pt '_' sz '_shu_2'];
vidfilename = fullfile(viddir, vidfn);
mkdir(viddir);
v=VideoWriter(vidfilename,'Motion JPEG AVI');
v.FrameRate = 15;
open(v);
writeVideo(v,F);
close(v);
end