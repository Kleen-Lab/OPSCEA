function [d, sfx, nns, ntp, eleclabels, LL, jumpto, scl, ts, em, depthch, axislim] = load_prepare_empty_data(pt, imagingpath)
global S;
global tiles;
global loaf;

% Instantiate variables that depend on seizure data to be empty
d = [];
LL = [];
sfx = 0;
ntp = 0;
scl = 1;
ts = [];
jumpto = 0;

%% locate and load electrode file for labels and XYZ coordinates
load(fullfile(imagingpath, pt, 'elecs', 'clinical_elecs_all.mat'), 'anatomy', 'elecmatrix', 'eleclabels');
if ~exist('anatomy','var')
    anatomy=cell(size(elecmatrix,1),4);
end
if size(anatomy,1)>size(elecmatrix,1)
    anatomy(size(elecmatrix,1)+1:end)=[];
end
anat=anatomy; clear anatomy;
if size(anat,2)>size(anat,1)
    anat=anat';
end
if size(anat,2)==1
    anat(:,2)=anat(:,1);
end
if ~exist('eleclabels','var')
    eleclabels=anat(:,1);
end
em=elecmatrix; clear elecmatrix;
EKGorREF=strcmpi('EKG1',anat(:,1))|strcmpi('EKG2',anat(:,1))|strcmpi('EKG',anat(:,2))|strcmpi('EKGL',anat(:,2))|strcmpi('REF',anat(:,1));
anat(EKGorREF,:)=[];
em(EKGorREF,:)=[];
eleclabels(EKGorREF,:)=[];

% When there is seizure data present, nns will depend on that data instead
nch = length(eleclabels);
nns = true(nch,1);

loaf.isR=nansum(em(:,1))>0; 
loaf.isL=loaf.isR~=1; %handy binary indicators for laterality

%% Implement isL/isR fix suggested by @aarongeller, allows the specifying of th side for all depths of bilateral implants
isRdepth = [];
isLdepth = [];

for i=1:height(tiles.depth)
    if ~isnan(tiles.depth.depths{i})
        xval_highcontact = em(tiles.depth.depths{i}(end),1);
        isRdepth(end+1) = xval_highcontact>=0;
        isLdepth(end+1) = xval_highcontact<0;
    else
        isRdepth(end+1) = nan;
        isLdepth(end+1) = nan;
    end
end
loaf.isRdepth = isRdepth;
loaf.isLdepth = isLdepth;

%% load meshes you want to plot
meshpath=fullfile(imagingpath, pt, 'Meshes');
Rcortex=load(fullfile(meshpath, [pt '_rh_pial.mat'])); 
loaf.rpial=Rcortex; 
Rcrtx=Rcortex.cortex; 
loaf.Rcrtx = Rcrtx;
clear Rcortex
Lcortex=load(fullfile(meshpath, [pt '_lh_pial.mat'])); 
loaf.lpial=Lcortex; 
Lcrtx=Lcortex.cortex; 
loaf.Lcrtx = Lcrtx;
clear Lcortex

for i=1:height(tiles.surface)
    hippentry(i)=~isempty(strfind(tiles.surface.surfaces,'hipp'));
    amygentry(i)=~isempty(strfind(tiles.surface.surfaces,'amyg')); 
end
errmsg='ATTN: MISSING A MESH, need to add this mesh file to directory (or remove/omit from frame): ';
if any(hippentry)
    Rhipp=fullfile(meshpath, 'subcortical', 'rHipp_subcort.mat'); 
    Lhipp=fullfile(meshpath, 'subcortical', 'lHipp_subcort.mat');
    if exist(Rhipp,'file')
        Rhipp=load(Rhipp); 
        Rhipp=Rhipp.cortex;
        loaf.Rhipp = Rhipp;
        Lhipp=load(Lhipp); 
        Lhipp=Lhipp.cortex; 
        loaf.Lhipp = Lhipp;
    else
        error([errmsg 'hipp']); 
    end
end
if any(amygentry)
    Ramyg=fullfile(meshpath, 'subcortical', 'rAmgd_subcort.mat'); 
    Lamyg=fullfile(meshpath, 'subcortical', 'lAmgd_subcort.mat');

    if exist(Ramyg,'file')
        Ramyg=load(Ramyg); 
        Ramyg=Ramyg.cortex;
        loaf.Ramyg = Ramyg;
        Lamyg=load(Lamyg); 
        Lamyg=Lamyg.cortex;
        loaf.Lamyg = Lamyg;
    else
        error([errmsg 'amyg']); 
    end
end

depthch=[]; 
for i=1:height(tiles.depth)
    depthch=[depthch tiles.depth.depths{i}]; 
end; clear i %identify all depth electrode channels

%% get xyz limits for plotting purposes
perim=1; % how many millimeters away from brain/electrodes boundaries to set the colorcoded plane perimeter, recommend >0 to avoid skimming brain surface (default 1mm)
axl(:,:,1)=[min([Rcrtx.vert; Lcrtx.vert]); max([Rcrtx.vert; Lcrtx.vert])]; %min and max of MESH VERTICES' x y z coordinates (2x3)
axl(:,:,2)=[min(em); max(em)]; %%min and max of ELECTRODES' x y z coordinates (2x3) and stack them (2x3x2)
axl=[min(axl(1,:,:),[],3); max(axl(2,:,:),[],3)]; % Get the minima and maxima of both
axislim=reshape(axl,1,6)+[-1 1 -1 1 -1 1]*perim; clear axl %Use the, to define the axis boundaries, and add additional perimeter (perim)
end