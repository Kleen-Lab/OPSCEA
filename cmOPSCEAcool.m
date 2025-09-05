function [cm]=cmOPSCEAcool
% colormap file using adapted cool.m function from Matlab
%
% Omni-planar and surface casting of epileptiform activity (OPSCEA)
% 
% Dr. Jon Kleen, 2017

cm=-cool(64)+1;
a=[0:1/32:1]'; 
a(end)=[]; 
cm(1:32,:)=cm(1:32,:).*[a a a]; 
cm=-[zeros(size(cm,1)+1,3); cm]+1;
end
