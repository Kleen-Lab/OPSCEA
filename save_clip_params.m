function finalClipID = save_clip_params(app, savedir)
    patientID = app.PatientIDEditField.Value;
    
    % Check that patientID is not empty
    if strrep(patientID, " ", "") == ""
        msgbox(sprintf("Please specify a patient ID"));
        finalClipID = '';
        return;
    end

    % Check that patient ID is valid (e.g. a folder has been
    % created for this patient). TODO(steph): is this needed? 
    % Should we always guarantee that valid patients will have
    % an existing folder (that may or may not be empty?)
    if exist('savedir', 'var')
        baseFolder = savedir;
    else
        baseFolder = fullfile(getenv("KLEEN_DATA"), 'opscea');
    end
    patientFolder = fullfile(baseFolder, patientID);
    if exist(patientFolder, "dir") ~= 7
        %msgbox(sprintf("Patient %s not found. Please specify another patient or create a new subplot configuration for this patient.", patientID));
        mkdir(patientFolder)
        finalClipID = '';
        return;
    end

    % Get list of folders for existing clips
    fileList = dir(patientFolder);
    folders = {fileList([fileList.isdir]).name};
    foldersToCheck = folders(startsWith(folders, patientID));
    idsOnly = strings(length(foldersToCheck));
    for i=1:length(foldersToCheck)
        idsOnly(i) = strrep(foldersToCheck(i), sprintf("%s_", patientID), "");
    end

    clipID = strrep(app.ClipIDEditField.Value, " ", "");
    if clipID == ""
        % If empty clip ID, set it to be the next highest number
        numericalIdsOnly = double(idsOnly);
        numericalIdsOnly = numericalIdsOnly(~isnan(numericalIdsOnly));
        if isempty(numericalIdsOnly)
            maxID = 0;
        else
            maxID = double(max(numericalIdsOnly));
        end

        finalClipID = maxID + 1;

        % Format with leading 0s
        if finalClipID < 10
            finalClipID = sprintf("0%d", finalClipID);
        else
            finalClipID = string(finalClipID);
        end

    elseif ismember(clipID, idsOnly)
        % Check that user intends to overwrite/update existing data
        message = sprintf("Video params with clip ID %s exists for patient %s. Do you want to overwrite the existing data?", clipID, patientID);
        btn1 = "Overwrite";
        btn2 = "Cancel";
        choice = questdlg(message, "Warning", btn1, btn2, btn2);
        if choice == "Cancel"
            msgbox("Configuration not saved.");
            finalClipID = '';
            return;
        else
            finalClipID = clipID;
        end
    
    else
        finalClipID = clipID;
    end
    
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
    patient = patientID;
    clipid = finalClipID;
    
    patientClipFolder = sprintf("%s_%s", patientID, finalClipID);
    folderPath = sprintf("%s/%s", patientFolder, patientClipFolder);
    if exist(folderPath, "dir") ~= 7
        mkdir(folderPath)
    end
    
    filename = sprintf("%s/%s_params.mat", folderPath, patientClipFolder);
    save(filename, "vidstart", "vidstop", "blstart", "blstop", "llw", "iceeg_scale",...
        "fps", "cax", "gsp", "cm", "iceegwin", "marg", "slicebright", "etype",...
        "patient", "clipid");

end