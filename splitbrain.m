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
    
    
    %mind the gap (in mm): prevents triangles that span the gap from inducing unwanted mesh
    thegap = max([abs(5.*m) 10]);
    
    % NOTE; occasionally produces unwanted side effects for certain
    % slices where a physical gap is inserted between the slice and cortex.
    % Adjusting thegap can help for individual cases but will need a more 
    % unified solution in a future iteration of OPSCEA software.
    % -------------------------------------------------------------------
    % Attempting a fix for this issue that was proposed by @aarongeller
    %   Adds the orientation_good function 
    
    if strcmp(orientation, 'c')
        if (min(cortex.cortex.vert(:,2)) < b) && (b < max(cortex.cortex.vert(:,2)))
            [idx, ~] = sort(find(abs((cortex.cortex.vert(:,2) - (m.*cortex.cortex.vert(:,1) + b + thegap))) <= thegap)); 
            % get indices of verts with y between (mx + b) and (mx + b + 2*thegap)
        else 
            % This means intercept lies outside A-P extent of the mesh
            disp('Intercept lies outside A-P extent of the mesh');
            [idx, ~] = sort(find(abs((cortex.cortex.vert(:,2) - (m.*cortex.cortex.vert(:,1) + b + thegap))+2*thegap)<=thegap));
        end
        
    elseif strcmp(orientation, 'a')
        if (min(cortex.cortex.vert(:,3)) < b) && (b < max(cortex.cortex.vert(:,3)))
            [idx, ~] = sort(find(abs((cortex.cortex.vert(:,3) - (m.*cortex.cortex.vert(:,1) + b - thegap))) <= thegap)); 
            % get indices of verts with z between (mx + b - 2*thegap) and (mx + b)
        else 
            % This means intercept lies outside S-I extent of the mesh
            disp('Intercept lies outside S-I extent of the mesh');
            [idx, ~] = sort(find(abs((cortex.cortex.vert(:,3)-(m.*cortex.cortex.vert(:,1) + b + thegap))+2*thegap)<=thegap));
        end
    
    elseif strcmp(orientation, 's')
        [idx, ~] = sort(find(abs((cortex.cortex.vert(:,2) - (m.*cortex.cortex.vert(:,1) + b - thegap))) <= thegap)); 
        % get indices of verts with y between (mx + b - 2*thegap) and (mx + b)
    
    elseif strcmp(orientation, 'oc')
        [idx, ~] = sort(find(abs((cortex.cortex.vert(:,3) - (m.*cortex.cortex.vert(:,2) + b + thegap))) <= thegap)); 
        % get indices of verts with z between (my + b + 2*thegap) and (my + b)
    end
    
    if isempty(idx)
        split.vert = cortex.cortex.vert;
        split.tri = cortex.cortex.tri;
    else
        mesh.tri = delete_verts(cortex.cortex.tri, idx);        
        mesh.vert = cortex.cortex.vert;
    
        disp('Generating partial mesh for slice view...')
        FVout = splitFV(mesh.tri,mesh.vert);
        fv.vert = FVout(1).vertices;
        fv.tri = FVout(1).faces;
        split = fv;
    end
end

function tri = delete_verts(tri, idx)
    while(~isempty(intersect(tri(:,2), idx)))
        [~,ia2,~] = intersect(tri(:,2), idx); %looks for the tris that match the indices of the verts on the line
        tri(ia2,:) = []; %remove the tris from the mesh dataset that have already been counted
    end
end
