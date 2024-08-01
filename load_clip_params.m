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