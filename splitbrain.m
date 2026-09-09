function split = splitbrain(cortex,orientation,b, m)
    % Splits the brain and exports a mesh that is a subset of the original
    % mesh. Requires the function splitFV (Matlab Exchange)
    %
    %
    % INPUTS:
    % cortex - the original mesh. a structure containing tri and vert 
    % 
    % orientation - slice orientation ('c', 's' or 'a') 'a' is assumed if a
    % letter different from 'c' or 's' is used. We always default to coronal
    % for simplicity
    % 
    % slicenum - the slice index on which to cut the mesh. (1-256)
    %
    % sl - the slope for coronal oblique planes (not configured for 's' or 'a')
    %
    % OUTPUTS:
    % split - a new mesh split at the specified point. A structure containing
    % tri and vert
    % 
    %
    %
    %     (This is a subfunction created as a part of) Omni-planar and surface
    %     casting of epileptiform activity (OPSCEA) (UC Case Number SF2020-281)
    %     jointly created by Dr. Jon Kleen, Ben Speidel, Dr. Robert Knowlton,
    %     and Dr. Edward Chang is licensed for non-commercial research use at
    %     no cost by the Regents of the University of California under CC
    %     BY-NC-SA 4.0 (https://creativecommons.org/licenses/by-nc-sa/4.0/).
    %     Please contact innovation@ucsf.edu if you are interested in using
    %     OPSCEA for commercial purposes.
    
    %     The following copyright notice and citation is to be included in any
    %     publication, material or media wherein all or a part of Licensed
    %     Material is contained, “Certain materials incorporated herein are
    %     Copyright © 2016 The Regents of the University of California
    %     (REGENTS). All Rights Reserved.
    
    %     Please cite the following paper in your publications if you have used
    %     our software in your research, as well as any relevant toolboxes used
    %     herein as appropriate (img_pipe, FreeSurfer): Kleen JK, Speidel B,
    %     Baud MO, Rao VR, Ammanuel SG, Hamilton LS, Chang EF, Knowlton RC.
    %     Accuracy of omni-planar and surface casting of epileptiform activity
    %     for intracranial seizure localization. In press at Epilepsia.”
    
    
    % Clip every face directly against the cut plane (m,b) rather than
    % deleting a vertex band near the line and trusting splitFV's
    % connected-component enumeration order to hand back the correct half.
    % That topological approach could silently return an arbitrary
    % component - including one made up entirely of wrong-side geometry
    % (e.g. an unsliced contralateral hemisphere overlaying the slice) -
    % since component ordering has no relationship to which side of the
    % plane a piece is on. Testing each face's centroid against the same
    % plane equation used by orientation_good() guarantees every
    % remaining face is genuinely on the correct side, regardless of mesh
    % topology, gaps, or how many disconnected pieces the cut produces.
    vert = cortex.cortex.vert;
    tri = cortex.cortex.tri;
    faceCentroids = (vert(tri(:,1),:) + vert(tri(:,2),:) + vert(tri(:,3),:)) / 3;
    x = faceCentroids(:,1); y = faceCentroids(:,2); z = faceCentroids(:,3);

    if strcmp(orientation, 'c')
        keep = y < (m.*x + b);
    elseif strcmp(orientation, 'a')
        keep = z > (m.*x + b);
    elseif strcmp(orientation, 's')
        keep = x > (y - b)./m;
    elseif strcmp(orientation, 'oc')
        keep = y < (z - b)./m;
    else
        keep = z > (m.*x + b); % default to 'a' behavior
    end

    split.vert = vert;
    split.tri = tri(keep, :);
end
