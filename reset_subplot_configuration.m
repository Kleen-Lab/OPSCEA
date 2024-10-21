function reset_subplot_configuration(app)

% Reset patient ID
app.PatientIDEditField.Value = "";

% Reset plot layout
app.ConfigurePlotLayoutEditField.Value = "6,100";

% Reset ECoG & colorbar configuration
elementNames = {"ECoG"; "Colorbar"};
rows = [6; 1];
cols = [20; 5];
tile = [81; 75];
app.ConfigurePlotElementsTable.Data = table(elementNames, tile, rows, cols);
app.ConfigurePlotElementsTable.Data.Properties.RowNames = ["ecog" "cb"];

% Reset surface and depth electrodes
rows = [];
cols = [];
tile = [];
zoom = [];
surfaces = [];
opacity = [];
view = [];
show = false(0, 0);
remove = false(0, 0);
app.SurfaceElectrodesTable.Data = table(tile, rows, cols, zoom, surfaces, opacity, view, show, remove);

rows = [];
cols = [];
tile = [];
zoom = [];
labels = [];
color = [];
depths = [];
remove = false(0, 0);
app.DepthElectrodesTable.Data = table(tile, rows, cols, zoom, labels, color, depths, remove);