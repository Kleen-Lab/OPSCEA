function load_clip_params(clipparams, test, app)
global S;

if exist('app', 'var')
    [vidstart, vidstop, blstart, blstop, llw, iceeg_scale, fps, cax, gsp, cm, iceegwin, marg, slicebright, etype] = process_clip_param_inputs(app);
    % Need to also load sfx
    
    pt = app.PatientIDEditField.Value;
    sz = app.ClipIDEditField.Value;
    ptsz = [pt '_' sz];
    datapath = getenv('KLEEN_DATA');
    opsceapath = fullfile(datapath, 'opscea'); %path for parameters sheet
    ptpath = fullfile(opsceapath, pt); % patient's folder
    szpath = fullfile(ptpath, ptsz); % specific seizure's folder
else
    load(clipparams, 'vidstart', 'vidstop', 'llw', 'iceeg_scale', 'fps', 'cax', 'gsp', 'cm', 'iceegwin', 'marg', 'slicebright');
    load(clipparams, 'blstart', 'blstop') 

    szpath = replace(clipparams, "/clip_params.mat", "");
    splitpath = split(szpath, '/');
    ptsz = splitpath{end};

end

load(fullfile(szpath, ptsz), 'sfx');

S.VIDperiod=[vidstart vidstop];
S.BLperiod=[blstart blstop];

% transform, scaling, and display options
% need defaults for if clip_params.mat is not fully configured yet (e.g.
% outputs from EEGHILITE only include vidstart, vidstop, blstart, blstop)
S.llw = exist('llw', 'var') && llw || 0.25;
S.iceeg_scale= exist('iceeg_scale', 'var') && iceeg_scale || 97; %percentile (number >50 and <100), used here similar to gain ICEEG waveform display, usually 95
S.fps = exist('fps', 'var') && fps || 5;
if test
    S.fps = 0.25;
end

%color axis for heatmap
if exist('cax', 'var')
    S.cax = cax;
else
    S.cax = [-20 20];
end

S.gsp = exist('gsp', 'var') && gsp || 25; %gaussian spreading parameter (default 10)
params={'iceeg_scale','fps','cax','gsp'};
paramsnans=isnan([(isnan(S.iceeg_scale) | S.iceeg_scale<=50 | S.iceeg_scale>=100)   S.fps   any(isnan(S.cax)) S.gsp]);
if any(paramsnans)
    error(['ATTENTION OPSCEA USER: The "' params{paramsnans} '" term(s) is/are in an incorrect format (perhaps number instead of string), check excel seizure parameter sheet']);
end
if ~exist('cm', 'var')
    cm = 'cmOPSCEAjet';
end
switch cm
    case 'cmOPSCEAcool'
        cm=cmOPSCEAcool;
    case 'cmOPSCEAjet'
        cm=cmOPSCEAjet;
end
S.cm = cm; %colormap to use for heatmap
S.iceegwin = exist('iceegwin', 'var') && iceegwin || 5; %how much trace-based ICEEG to view at a time in the ICEEG window
S.marg = exist('marg', 'var') && marg || 1; %offset of real-time LL txform from beginning of viewing window (in sec; converts to samples below)
S.slicebright = exist('slicebright', 'var') && slicebright || 25;
if isnan(S.slicebright); S.slicebright=0; end %brighten up slices (usually 0 to 50)

% additional adjustment for display window
S.VIDperiod=[S.VIDperiod(1)-S.marg   S.VIDperiod(2)+S.iceegwin-S.marg];

S.sliceplane='c'; % calculate omni-planar slice angles with respect to coronal (c) plane

S=orderfields(S); %alphabetize the structure fields for ease of use/search

% Scaling ICEEG and setting windows for simultaneous trace-based display
S.iceeg_scale=S.iceeg_scale+(100-S.iceeg_scale)/2; 
S.iceeg_scale=[100-S.iceeg_scale S.iceeg_scale]; %conversion to two-tailed percentile

S.marg=round(S.marg*sfx); %offset of real-time LL txform from beginning of viewing window
S.fram=round(sfx/S.fps);

clear llw iceeg_scale fps cax gsp cm iceegwin marg slicebright;
end