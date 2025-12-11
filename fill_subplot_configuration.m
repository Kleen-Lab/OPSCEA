function fill_subplot_configuration(patientParams, app)

%%%% Layout configuration %%%%
app.ConfigurePlotLayoutEditField.Value = sprintf(...
    "%d,%d",...
    patientParams.layout.rows,...
    patientParams.layout.cols...
);
    
%%%% ECoG and colorbar configuration %%%% 
elementNames = {"ECoG"; "Colorbar"};
rows = [int32(patientParams.ecog.rows); int32(patientParams.cb.rows)];
cols = [int32(patientParams.ecog.cols); int32(patientParams.cb.cols)];
tile = [int32(patientParams.ecog.tile); int32(patientParams.cb.tile)];
app.ConfigurePlotElementsTable.Data = table(elementNames, tile, rows, cols);
app.ConfigurePlotElementsTable.Data.Properties.RowNames = ["ecog" "cb"];

%%%% Surface electrodes %%%%
% Initialize column for removing electrodes
[numElectrodes, ~] = size(patientParams.surface);
removeElectrodesColumn = false(numElectrodes, 1);

% Fix display for opacity column
opacityColumn = strings(numElectrodes, 1);
for i=1:numElectrodes
        opacityColumn(i) = strjoin(arrayfun(@num2str, patientParams.surface.opacity{i}, 'UniformOutput', false), ',');
end

surfaceElectrodesTable = table( ...
    patientParams.surface.tile, ...   
    patientParams.surface.rows, ...
    patientParams.surface.cols, ...
    patientParams.surface.zoom, ...
    patientParams.surface.surfaces, ...
    opacityColumn, ...
    patientParams.surface.view, ...
    patientParams.surface.show, ...
    removeElectrodesColumn...
  );
app.SurfaceElectrodesTable.Data = surfaceElectrodesTable;
app.SurfaceElectrodesTable.Data.Properties.VariableNames = [
    "tile", "rows", "cols", "zoom", "surfaces", "opacity", "view", "show", "remove"
];

%%%% Depth electrodes %%%%
% Initialize column for removing electrodes
[numElectrodes, ~] = size(patientParams.depth);
removeElectrodesColumn = false(numElectrodes, 1);

% Fix display for color column
colorColumn = strings(numElectrodes, 1);
for i=1:numElectrodes
    colorColumn(i) = strjoin(arrayfun(@num2str, transpose(patientParams.depth.color(i, :)), 'UniformOutput', false), ',');
end

% Fix display for depths column
depthColumn = strings(numElectrodes, 1);
for i=1:numElectrodes
    % depthColumn(i) = strjoin(string(patientParams.depth.depths{i}), ',');
    depth_i = patientParams.depth.depths{i};
    depthColumn(i) = string(depth_i(1)) + ':' + string(depth_i(end));
end

if size(patientParams.depth, 1) > 0
    depthElectrodesTable = table( ...
        patientParams.depth.tile, ...
        patientParams.depth.rows, ...
        patientParams.depth.cols, ...
        patientParams.depth.zoom, ...
        patientParams.depth.labels, ...
        colorColumn, ...
        depthColumn, ...
        removeElectrodesColumn...
      );
    app.DepthElectrodesTable.Data = depthElectrodesTable;
else
    % Set depth electrode table to be blank
    rows = [];
    cols = [];
    tile = [];
    zoom = [];
    labels = [];
    color = [];
    depths = [];
    remove = false(0, 0);
    app.DepthElectrodesTable.Data = table(tile, rows, cols, zoom, labels, color, depths, remove);
end
app.DepthElectrodesTable.Data.Properties.VariableNames = [
    "tile", "rows", "cols", "zoom", "labels", "color", "depths", "remove"
];