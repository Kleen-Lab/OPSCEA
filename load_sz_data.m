function load_sz_data(pt_id, clip_id)
    arguments
        pt_id (1, :) char;
        clip_id (1, :) char;
    end

    setenv("OPSCEA_PT_ID", pt_id);
    setenv("OPSCEA_CLIP_ID", clip_id)
    
    eeghilite = EEGHILITE;
    uiwait(eeghilite);

    
end