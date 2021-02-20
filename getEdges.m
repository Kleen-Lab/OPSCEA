function [xedge, yedge, zedge] = getEdges(alphamap, XX, YY, ZZ)
%find the first row that has a nonzero
%
% Omni-planar and surface casting of epileptiform activity (OPSCEA)
% 
% Dr. Jon Kleen, 2017

i=1;
allblack = 1;
while allblack==1

    if size(find(alphamap(i,:)==1))>0
        firstrow = i;
        allblack = 0;
    end
    
    i = i+1;
end

%find the first column that has a nonzero
j=1;
allblack=1;
while allblack==1

    if size(find(alphamap(:,j)==1))>0
        firstcol = j;
        allblack = 0;
    end
    
    j = j+1;
end

%same thing coming from the bottom
i = size(alphamap,1);
allblack = 1;
while allblack==1

    if size(find(alphamap(i,:)==1))>0
        lastrow = i;
        allblack = 0;
    end
    
    i = i-1;
end

j=size(alphamap,2);
allblack=1;
while allblack==1

    if size(find(alphamap(:,j)==1))>0
        lastcol = j;
        allblack = 0;
    end
    
    j = j-1;
end

% corner1 = [firstrow, firstcol];
% corner2 = [firstrow, lastcol];
% corner3 = [lastrow, firstcol];
% corner4 = [lastrow, lastcol];

xedge = ([XX(firstrow, firstcol) XX(firstrow, lastcol)]);
yedge = ([YY(firstrow, firstcol) YY(firstrow, lastcol)]);
zedge = ([ZZ(lastrow, firstcol) ZZ(firstrow, firstcol)]);

end
