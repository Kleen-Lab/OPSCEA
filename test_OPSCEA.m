function test_OPSCEA()
    patient_list = readtable("OPSCEAparams.xlsx", "Sheet", "params", 'DataRange', 'A2:A240');
    seizure_list = readtable("OPSCEAparams.xlsx", "Sheet", "params", "DataRange", "B2:B240");
    
    for i=1:height(patient_list)
        pt = string(patient_list{i,1});
        sz = string(seizure_list{i,1});
        OPSCEA(pt, sz, 1, 0)
    end
end