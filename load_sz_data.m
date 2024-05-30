function [d, sfx, badch] = load_sz_data(pt_id)
    dt_prompts = {
        'Start date of the ICEEG recording: ', ...
        'Start datetime of the seizure: ', ...
        'End datetime of the seizure: '
    };

    dt_title = 'Seizure time information';

    dt_definput = {'yyyy-MM-dd', 'yyyy-MM-dd_HH:mm:ss', 'yyyy-MM-dd_HH:mm:ss'};

    dt_fieldsize = [1 50; 1 50; 1 50];

    dt = inputdlg(dt_prompts, dt_title, dt_fieldsize, dt_definput);

    iceeg_date = dt{1};
    start_dt = dt{2};
    end_dt = dt{3};

    [d, sfx] = fetch_opscea_clip(pt_id, start_dt, end_dt, iceeg_date);
    badch = [];

    [d, sfx] = preprocess_clip(d, sfx);
end