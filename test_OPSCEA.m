function test_OPSCEA()
    patient_list = readtable("OPSCEAparams.xlsx", "Sheet", "params", 'DataRange', 'A2:A240');
    seizure_list = readtable("OPSCEAparams.xlsx", "Sheet", "params", "DataRange", "B2:B240");
    
    for i=35:height(patient_list)
        pt = char(patient_list{i,1});
        sz = char(seizure_list{i,1});
        if strcmp(pt, 'TD143_singlebrain');continue;end;
        OPSCEA(pt, sz, 1, 0)
    end
end