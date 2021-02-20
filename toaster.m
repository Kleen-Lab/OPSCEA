function [cim]=toaster(XX,YY,ZZ,im,em,w8,cax,cm,gsp)

% recolors a 2D grayscale slice with a heatmap according to electrode
% locations in 3D space (similar proximity principles as ctmr_gauss_plot)
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

%interpolate grid points for centered coordinates
x=(XX(1:end-1,:)+XX(2:end,:))/2; x=(x(:,1:end-1)+x(:,2:end))/2;
y=(YY(1:end-1,:)+YY(2:end,:))/2; y=(y(:,1:end-1)+y(:,2:end))/2;
z=(ZZ(1:end-1,:)+ZZ(2:end,:))/2; z=(z(:,1:end-1)+z(:,2:end))/2;

%vectorize
im=reshape(im,1,256^2); x=reshape(x,1,256^2); y=reshape(y,1,256^2); z=reshape(z,1,256^2); 

c=zeros(length(im),1);
for i=1:length(em(:,1))
    bx=abs(x-em(i,1));    by=abs(y-em(i,2));    bz=abs(z-em(i,3)); 
    c=c+(w8(i)*exp((-(bx.^2+by.^2+bz.^2))/gsp))'; %gaussian fall off
end

hcm=size(cm,1)/2; %use the same colormap used for surfaces
c=c-cax(1); c=c/diff(cax); %scales weights to the range of cax input
c(c<0)=0; c(c>1)=1; % hedge any values out of the cax range
c=round(c*hcm+hcm); % using 1/2 of the colormap (eg. if looking at absolute weights)
cidx=cm(c,:); %pixels' weighted values now act as colormap indices
cidx=reshape(cidx,256,256,3);
im=repmat(im',1,3)/255; %normalize equivalent to colormap (0 to 1)
im=reshape(im,256,256,3); %return to a matrix, and set up for uint8 format
cim=round(cidx.*double(im)*255);%combine them
end