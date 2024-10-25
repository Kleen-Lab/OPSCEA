function [layout, ecog, cb, surface, depth] = process_patient_param_inputs(app)
% Patient params: layout
splitStr = split(app.ConfigurePlotLayoutEditField.Value, ",");
layout.rows = str2double(splitStr(1));
layout.cols = str2double(splitStr(2));
            
% Patient params: ecog
ecogFields = app.ConfigurePlotElementsTable.Data{"ecog", 2:4};
rows = double(ecogFields(2));
cols = double(ecogFields(3));
tile = double(ecogFields(1));
ecog = table(rows, cols, tile);

% Patient params: cb 
cbFields = app.ConfigurePlotElementsTable.Data{"cb", 2:4};
rows = double(cbFields(2));
cols = double(cbFields(3));
tile = double(cbFields(1));
cb = table(rows, cols, tile);

% Patient params: surface electrodes
surfaceElectrodesTable = app.SurfaceElectrodesTable.Data;
rows = double(surfaceElectrodesTable.rows);
cols = double(surfaceElectrodesTable.cols);
tile = double(surfaceElectrodesTable.tile);
zoom = double(surfaceElectrodesTable.zoom);
surfaces = surfaceElectrodesTable.surfaces;
view = surfaceElectrodesTable.view;
show = surfaceElectrodesTable.show;

% Need to convert comma-delimited string to array
[numRows, ~] = size(app.SurfaceElectrodesTable.Data);
opacity = cell(numRows, 1);
for i=1:numRows
    opacity{i, 1} = transpose(double(split(surfaceElectrodesTable.opacity(i), ",")));
end

surface = table(zoom, surfaces, opacity, view, show, rows, cols, tile);

% Patient params: depth electrodes
depthElectrodesTable = app.DepthElectrodesTable.Data;
rows = double(depthElectrodesTable.rows);
cols = double(depthElectrodesTable.cols);
tile = double(depthElectrodesTable.tile);
zoom = double(depthElectrodesTable.zoom);
labels = depthElectrodesTable.labels;

% Need to convert comma-delimited string to array
[numRows, ~] = size(app.DepthElectrodesTable.Data);
depths = cell(numRows, 1);
color = cell(numRows, 3);
for i=1:numRows
    depths{i, 1} = transpose(double(split(depthElectrodesTable.depths(i), ",")));
    color(i, :) = num2cell(double(split(depthElectrodesTable.color(i), ",")));
end
color = double(string(color));

depth = table(zoom, labels, color, depths, rows, cols, tile);