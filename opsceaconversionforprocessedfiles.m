


curdir=cd; 
fld=uigetdir(curdir,'Select seizure folder (DeSpiker folder of processed files should be within this)'); 
if fld==0; return; else; cd(fld); end
fld=[fld '/'];

pathfolders=regexp(fld,'/','split');
pt_sz=pathfolders{end-1};

edffile = dir('*.edf');
edffile=edffile.name;
if size(edffile,1)>1; error('more than one EDF file present'); end

dsfld=[fld 'DeSpiker.' edffile(1:end-4) '/'];

movefile([dsfld 'ppEEG.mat'],[fld pt_sz '.mat'])
movefile([dsfld 'bad_chs.mat'],[fld pt_sz '_badch.mat'])


load([fld pt_sz '.mat']);
d=ppEEG; sfx=fs; clear ppEEG fs study
save([fld pt_sz '.mat'],'d','sfx');

load([fld pt_sz '_badch.mat']);
badch=logical(bad_chs); clear bad_chs
save([fld pt_sz '_badch.mat'],'badch');
