function reset_loaded_seizure_data(app)

app.d = [];
app.LL = [];
app.sfx = 0;
app.ntp = 0;
app.scl = 1;
app.ts = [];
app.jumpto = 1;
app.nns = [];
app.currentLoadedClip = '';

app.DataLoadedInfoLabel.Text = 'Seizure data currently loaded for: None';