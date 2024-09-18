function [vidstart, vidstop, blstart, blstop, llw, iceeg_scale, fps, cax, gsp, cm, iceegwin, marg, slicebright, etype] = process_clip_param_inputs(app)
vidstart = double(app.ClipIDConfigTable.Data{"vidstart", 2}{1});
vidstop = double(app.ClipIDConfigTable.Data{"vidstop", 2}{1});
blstart = double(app.ClipIDConfigTable.Data{"blstart", 2}{1});
blstop = double(app.ClipIDConfigTable.Data{"blstop", 2}{1});
llw = double(app.ClipIDConfigTable.Data{"llw", 2}{1});
iceeg_scale = double(app.ClipIDConfigTable.Data{"iceeg_scale", 2}{1});
fps = double(app.ClipIDConfigTable.Data{"fps", 2}{1});
cax = double(split(app.ClipIDConfigTable.Data{"cax", 2}{1}, ","));
gsp = double(app.ClipIDConfigTable.Data{"gsp", 2}{1});
cm = app.ClipIDConfigTable.Data{"cm", 2}{1};
iceegwin = double(app.ClipIDConfigTable.Data{"iceegwin", 2}{1});
marg = double(app.ClipIDConfigTable.Data{"marg", 2}{1});
slicebright = double(app.ClipIDConfigTable.Data{"slicebright", 2}{1});
etype = app.ClipIDConfigTable.Data{"etype", 2}{1};