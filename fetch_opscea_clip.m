function [d, sfx] = fetch_opscea_clip(pt_id,start_dt,end_dt,iceeg_dt,python_path,py_module_path)
    arguments
        pt_id (1,:) char;
        start_dt (1,19) char;
        end_dt (1,19) char;
        iceeg_dt char = '';
        python_path (1,:) char = '/scratch/fetch_clip/.venv/bin/python';
        py_module_path (1, :) char = '/scratch/fetch_clip';
    end

    %setup python
    pyenv(Version=python_path, ExecutionMode="OutOfProcess");
    insert(py.sys.path, int64(0), py_module_path);
    mod = py.importlib.import_module('opscea_fetch');
    
    %fetch the clip
    if isempty(iceeg_dt)
        pystruct = mod.fetch_clip(pt_id, start_dt, end_dt);
    else
        pystruct = mod.fetch_clip(pt_id, start_dt, end_dt, iceeg_dt);
    end
    %unpack and convert python data to matlab objects
    py_d = pystruct{1};
    py_sfx = pystruct{2};
    d = double(py_d);
    sfx = double(py_sfx);
end