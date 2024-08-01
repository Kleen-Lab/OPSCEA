function load_patient_params(ptparams)
global tiles;

load(ptparams, 'layout', 'ecog', 'cb', 'surface', 'depth');
tiles.layout = layout; clear layout;
tiles.ecog = ecog; clear ecog;
tiles.cb = cb; clear cb;
tiles.surface = surface; clear surface;
tiles.depth = depth; clear depth;
end