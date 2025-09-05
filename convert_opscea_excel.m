function convert_opscea_excel()
dataroot = getenv('KLEEN_DATA');
opscea_path = fullfile(dataroot, 'opscea');
file = fullfile(opscea_path, 'OPSCEAparams.xlsx');
opts = detectImportOptions(file);
params_opts = opts;
params_opts.Sheet = 'params';
clip_params = readtable(file, params_opts);

nclips = height(clip_params);
for row = 1:nclips
    % extract clip params
    params = clip_params(row, :);
    patient = char(params.patient);
    if strcmpi(patient, 'EC266')
        continue
    end
    clipid = char(params.sz);
    vidstart = str2double(params.VIDstart);
    vidstop = str2double(params.VIDstop);
    blstart = str2double(params.Blstart);
    if isnan(blstart)
        blstart = params.Blstart;
    end
    blstop = str2double(params.Blstop);
    if isnan(blstop)
        blstop = params.Blstop;
    end
    llw = str2double(params.llw);
    iceeg_scale = str2double(params.iceeg_scale);
    fps = str2double(params.fps);
    cax = str2double(split(cell2mat(params.cax), ','));
    gsp = str2double(params.gsp);
    cm = char(params.cm);
    iceegwin = str2double(params.iceegwin);
    marg = str2double(params.marg);
    slicebright = str2double(params.slicebright);
    etype = char(params.etype);
    chunk = params.chunk;
    blstretch = str2double(split(cell2mat(params.BLstretch), ','));
    zeach = str2double(params.zeach);
    clip_params_file = fullfile(opscea_path, patient, [patient '_' clipid], 'clip_params.mat');
    save(clip_params_file, 'patient', 'clipid', 'vidstart', 'vidstop', 'blstart', 'blstop',...
        'llw', 'iceeg_scale', 'fps', 'cax', 'gsp', 'cm', 'iceegwin', 'marg',...
        'slicebright', 'etype', 'chunk', 'blstretch', 'zeach');

    % extract params
    if strcmpi(patient, 'EC287')
        continue
    end
    patient_params_file = fullfile(opscea_path, patient, 'patient_params.mat');
    patient_opts = opts;
    patient_opts.Sheet = patient;
    patient_params = readtable(file, patient_opts);
    ecog = patient_params(strcmp(patient_params.plottype, 'ecog'), 2:5);
    ecog.Properties.VariableNames = ["row", "col", "start", "stop"];
    ecog.row = str2double(ecog.row);
    ecog.col = str2double(ecog.col);
    ecog.start = str2double(ecog.start);
    if isstring(ecog.stop)
        ecog.stop = str2double(ecog.stop);
    end
    cb = patient_params(strcmp(patient_params.plottype, 'colorbar'), 2:5);
    cb.Properties.VariableNames = ["row", "col", "start", "stop"];
    cb.row = str2double(cb.row);
    cb.col = str2double(cb.col);
    cb.start = str2double(cb.start);
    if isstring(cb.stop)
        cb.stop = str2double(cb.stop);
    end
    surface = patient_params(strcmp(patient_params.plottype, 'surface'), 2:10);
    surface.Properties.VariableNames = ["row", "col", "start", "stop", "zoom", ...
                                        "surfaces", "opacity", "view", "show"];
    surface.row = str2double(surface.row);
    surface.col = str2double(surface.col);
    surface.start = str2double(surface.start);
    if isstring(surface.stop)
        surface.stop = str2double(surface.stop);
    end
    if isstring(surface.zoom)
        surface.zoom = str2double(surface.zoom);
    end
    for i=1:height(surface)
        surf = surface(i, :);
        surface.opacity{i} = double(split(string(surf.opacity), ','));
        surface.show{i} = ~isempty(surf.show{1});
    end
    depth = patient_params(strcmp(patient_params.plottype, 'depth'), [2:6, 11:14]);
    depth.Properties.VariableNames = ["row", "col", "start", "stop", "zoom", ...
                                      "labels", "first", "last", "color"];
    depth.row = str2double(depth.row);
    depth.col = str2double(depth.col);
    depth.start = str2double(depth.start);
    if isstring(depth.stop)
        depth.stop = str2double(depth.stop);
    end
    if isstring(depth.zoom)
        depth.zoom = str2double(depth.zoom);
    end
    depth.first = str2double(depth.first);
    depth.last = str2double(depth.last);
    for i=1:height(depth)
        dpth = depth(i, :);
        depth.depths{i} = dpth.first:dpth.last;
    end
    depth.first = [];
    depth.last = [];
    color = split(string(depth.color), ',');
    if isequal(size(color), [3 1])
        color = color';
    end
    depth.color = double(color);
    sprows = [ecog.row; cb.row; surface.row; depth.row];
    spcols = [ecog.col; cb.col; surface.col; depth.col];
    spstarts = [ecog.start; cb.start; surface.start; depth.start];
    spstops = [ecog.stop; cb.stop; surface.stop; depth.stop];
    widths = spstops - spstarts + 1;
    [~, d] = rat(widths./spcols);
    layout.rows = max(lcm(sprows, sprows(end)));
    try
        layout.cols = max(lcm(d, d(end)));
    catch ME
        rethrow(ME);
    end
    ecog = convert_subplot_to_nexttile(ecog, layout);
    cb = convert_subplot_to_nexttile(cb, layout);
    surface = convert_subplot_to_nexttile(surface, layout);
    depth = convert_subplot_to_nexttile(depth, layout);
    save(patient_params_file, 'layout', 'ecog', 'cb', 'surface', 'depth');
end
end

function table = convert_subplot_to_nexttile(table, layout)
for i=1:height(table)
    plot = table(i, :);
    row_ratio = layout.rows/plot.row;
    col_ratio = layout.cols/plot.col;
    width_orig = plot.stop - plot.start + 1;
    [n, d] = rat(width_orig/plot.col);
    n_preceding_orig_rows = ceil(plot.start/plot.col)-1;

    table.rows(i) = row_ratio;
    table.cols(i) = layout.cols/d*n;
    table.tile(i) = n_preceding_orig_rows*row_ratio*layout.cols+mod(plot.start-1, plot.col)*col_ratio+1;
end
table.row = [];
table.col = [];
table.start = [];
table.stop = [];
end