function litebrain(direc,bri) %Jon Kleen 2018
% direc is direction from where the light should originate. multiple ok.
% bri is brightness from 0 (none) to 1 (max)
%
% Omni-planar and surface casting of epileptiform activity (OPSCEA)
% 
% Dr. Jon Kleen, 2017

if nargin<1; direc = 'd';end
if nargin<2; bri = 0.8; end
bri=max([bri 0]); bri=min([bri 1]); %cannot be less than 0 or more than 1
for j=1:length(direc)
switch lower(direc(j))
    case 'l'; v=[270   0]; psn=[-1 0 0];  %left
    case 'r'; v=[ 90   0]; psn=[1  0 0];  %right
    case 'a'; v=[180   0]; psn=[0  1 0];  %anterior
    case 'p'; v=[  0   0]; psn=[0 -1 0];  %posterior
    case 's'; v=[  0  90]; psn=[0  0 1];  %superior
    case 'i'; v=[180 -90]; psn=[0 0 -1];  %inferior
    case 'd'; v=[  0 -90]; psn=[0 0 -1];  %default diagnonal (antero-superior)
end
V(j,:)=v;
P(j,:)=psn;
end
%merges to a single angle for both, in case multiple entered
V=mean(V,1);      P=mean(P,1); 
% apply
if strcmpi(direc,'i'); view(V(1),V(2)); else; view(P); end %make sure anterior at top for cardinal inferior view 
if bri; set(light,'position',P,'color',bri*[1 1 1]); end


%% Extra tips: 
% view inferior but rotated so that in the figure it is positioned:
%             - anterior to the right: view(90,270)
%             - anterior to the left: view(270,270)
%             - anterior to the bottom: view(0,270)
