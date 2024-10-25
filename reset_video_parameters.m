function reset_video_parameters(app)

% Reset clip ID
app.ClipIDEditField.Value = "";

% Reset table for video clip configuration
paramNames = {
    "Video Start"; 
    "Video Stop";
    "Baseline Start"; 
    "Baseline Stop";
    "Line Length Window (sec)"; 
    "ICEEG Scale";
    "FPS"; 
    "Color Axis";
    "GSP"; 
    "Color Map";
    "ICEEG Window"; 
    "Margin (sec)";
    "Slice Brightness"; 
    "Electrode Type"
};
defaultValues = {
    0;
    0;
    0; 
    0;
    0.25; 
    97; 
    5; 
    "-20,20";
    25; 
    "cmOPSCEAjet";
    5; 
    1;
    25; 
    "clinical";
};
app.ClipIDConfigTable.Data = table(paramNames, defaultValues);
app.ClipIDConfigTable.Data.Properties.RowNames = [
    "vidstart"
    "vidstop"
    "blstart"
    "blstop"
    "llw"
    "iceeg_scale"
    "fps"
    "cax"
    "gsp"
    "cm"
    "iceegwin"
    "marg"
    "slicebright"
    "etype"
];
