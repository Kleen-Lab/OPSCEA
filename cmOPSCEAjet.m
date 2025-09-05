function [cm]=cmOPSCEAjet
% colormap file using adapted cool.m function from Matlab, 
% adjusted for a more even spread across colors
% Omni-planar and surface casting of epileptiform activity (OPSCEA)
% 
% Jon Kleen, 2017

cm=[ones(64,3); jet(64)]; 
cm(65:76,:)=[]; 
cm(65:80,:)=[linspace(1,cm(80,1),16); linspace(1,cm(80,2),16); linspace(1,cm(80,3),16)]'; 
scm=size(cm,1); cm([1:2:round(scm*.7) round(scm*.85):2:end],:)=[];
cm=[ones(8,3); cm]; 