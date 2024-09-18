function load_patient_params(ptparams, app)
global tiles;

if exist('app', 'var')
    [layout, ecog, cb, surface, depth] = process_patient_param_inputs(app);
else
    load(ptparams, 'layout', 'ecog', 'cb', 'surface', 'depth');
end

tiles.layout = layout; clear layout;
tiles.ecog = ecog; clear ecog;
tiles.cb = cb; clear cb;
tiles.surface = surface; clear surface;
tiles.depth = depth; clear depth;
end