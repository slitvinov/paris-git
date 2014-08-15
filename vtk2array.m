function A = vtk2array(folderloc, time)
%vtk2array
%   Returns a three dimensional array of the color function at a specific
%   point in time.  To use, simply input the VTK folder and desired
%   timestep.  This should work for variables other than the color
%   function, but only if the value is described at the center of the grid
%   cell.  Please feel free to contact me with any comments/concerns.
%   
%   Sample Call:
%   folderloc = '/home/jakedynes/Projects/paris-devel/Tests/VOF/out/VTK/';
%   time = '00000';
%   A = vtk2array(folderloc, time)
%
%   Jake Dynes, August 2014 (dynesj1@montclair.edu)
    
for p=0:7
    
    % import data from the current processor
    intProcessor=strcat('0000',int2str(p));
    strFileloc=strcat(folderloc, 'VOF', time, '-', intProcessor, '.vtk');
    matData=importdata(strFileloc);
    
    % some setup from first processor
    if p==0

        % get grid size (specified in 'input')
        intXgriddim=(str2double(matData.textdata{5,2})-1)*2;
        intYgriddim=(str2double(matData.textdata{5,3})-1)*2;
        intZgriddim=(str2double(matData.textdata{5,4})-1)*2;

        % create VOF array (dimensions)
        matVOF=zeros(intZgriddim,intXgriddim,intYgriddim);
    end
    
    % create raw data matrices (points and values)      
    intTotalCellswithGhost=str2double(matData.textdata{5,2})*str2double(matData.textdata{5,3})*str2double(matData.textdata{5,4});
    matPoints=zeros(intTotalCellswithGhost,3);

    for i=1:intTotalCellswithGhost
        for j=1:3    
            matPoints(i,j)=str2double(matData.textdata{i+6,j});
        end
    end

    matValues=zeros(intTotalCellswithGhost,1);
    intStartCell=6+intTotalCellswithGhost+3;
    
    for i=1:intTotalCellswithGhost
        matValues(i,1)=str2double(matData.textdata{intStartCell+i,1});
    end
    
    % more setup if on first processor
    if p==0
    
        % find domain dimensions
        xDomain=intXgriddim*matPoints(1,1);
        yDomain=intYgriddim*matPoints(1,2);
        zDomain=intZgriddim*matPoints(1,3);
    
    end
    
    % sort out ghost values
    intNumInVOF=(intXgriddim/2)*(intYgriddim/2)*(intZgriddim/2);
    matCorrectVOF=zeros(intNumInVOF,1);
    j=1;
    for i=1:intTotalCellswithGhost
        xOK=0;
        yOK=0;
        zOK=0;
        
        % tedious but necessary formatting for current processor
        if p==0 % front left bottom processor
            if and(matPoints(i,1)>0,matPoints(i,1)<xDomain)
                xOK=1;
            end
            if and(matPoints(i,2)>0,matPoints(i,2)<yDomain)
                yOK=1;
            end
            if and(matPoints(i,3)>0,matPoints(i,3)<zDomain)
                zOK=1;
            end
        elseif p==1 % front left top processor   
            if and(matPoints(i,1)>0,matPoints(i,1)<xDomain)
                xOK=1;
            end
            if and(matPoints(i,2)>0,matPoints(i,2)<yDomain)
                yOK=1;
            end
            if and(matPoints(i,3)>zDomain,matPoints(i,3)<zDomain*2)
                zOK=1;
            end
        elseif p==2 % back left bottom processor
            if and(matPoints(i,1)>0,matPoints(i,1)<xDomain)
                xOK=1;
            end
            if and(matPoints(i,2)>yDomain,matPoints(i,2)<yDomain*2)
                yOK=1;
            end
            if and(matPoints(i,3)>0,matPoints(i,3)<zDomain)
                zOK=1;
            end
        elseif p==3 % back left top processor
            if and(matPoints(i,1)>0,matPoints(i,1)<xDomain)
                xOK=1;
            end
            if and(matPoints(i,2)>yDomain,matPoints(i,2)<yDomain*2)
                yOK=1;
            end
            if and(matPoints(i,3)>zDomain,matPoints(i,3)<zDomain*2)
                zOK=1;
            end
        elseif p==4 % front right bottom processor
            if and(matPoints(i,1)>xDomain,matPoints(i,1)<xDomain*2)
                xOK=1;
            end
            if and(matPoints(i,2)>0,matPoints(i,2)<yDomain)
                yOK=1;
            end
            if and(matPoints(i,3)>0,matPoints(i,3)<zDomain)
                zOK=1;
            end
        elseif p==5 % front right top processor
            if and(matPoints(i,1)>xDomain,matPoints(i,1)<xDomain*2)
                xOK=1;
            end
            if and(matPoints(i,2)>0,matPoints(i,2)<yDomain)
                yOK=1;
            end
            if and(matPoints(i,3)>zDomain,matPoints(i,3)<zDomain*2)
                zOK=1;
            end
        elseif p==6 % back right bottom processor
            if and(matPoints(i,1)>xDomain,matPoints(i,1)<xDomain*2)
                xOK=1;
            end
            if and(matPoints(i,2)>yDomain,matPoints(i,2)<yDomain*2)
                yOK=1;
            end
            if and(matPoints(i,3)>0,matPoints(i,3)<zDomain)
                zOK=1;
            end
        elseif p==7 % back right top processor
            if and(matPoints(i,1)>xDomain,matPoints(i,1)<xDomain*2)
                xOK=1;
            end
            if and(matPoints(i,2)>yDomain,matPoints(i,2)<yDomain*2)
                yOK=1;
            end
            if and(matPoints(i,3)>zDomain,matPoints(i,3)<zDomain*2)
                zOK=1;
            end
        end
        
        A=[xOK,yOK,zOK];
        if all(A)
            matCorrectVOF(j)=matValues(i,1);
            j=j+1;            
        end        
        
    end
    
    % now siphon the values from matCorrectVOF into matVOF
    % again, tedious but necessary formatting
    l=1;
    if p==0
        for k=intZgriddim:-1:(intZgriddim/2)+1
            for j=1:(intYgriddim/2)
                for i=1:(intXgriddim/2)
                    matVOF(k,i,j)=matCorrectVOF(l);
                    l=l+1;
                end
            end
        end
    elseif p==1
        for k=(intZgriddim/2):-1:1
            for j=1:(intYgriddim/2)
                for i=1:(intXgriddim/2)
                    matVOF(k,i,j)=matCorrectVOF(l);
                    l=l+1;
                end
            end
        end
    elseif p==2
        for k=intZgriddim:-1:(intZgriddim/2)+1
            for j=(intYgriddim/2)+1:intYgriddim
                for i=1:(intXgriddim/2)
                    matVOF(k,i,j)=matCorrectVOF(l);
                    l=l+1;
                end
            end
        end
    elseif p==3
        for k=(intZgriddim/2):-1:1
            for j=(intYgriddim/2)+1:intYgriddim
                for i=1:(intXgriddim/2)
                    matVOF(k,i,j)=matCorrectVOF(l);
                    l=l+1;
                end
            end
        end
    elseif p==4
        for k=intZgriddim:-1:(intZgriddim/2)+1
            for j=1:(intYgriddim/2)
                for i=(intXgriddim/2)+1:intXgriddim
                    matVOF(k,i,j)=matCorrectVOF(l);
                    l=l+1;
                end
            end
        end
    elseif p==5
        for k=(intZgriddim/2):-1:1
            for j=1:(intYgriddim/2)
                for i=(intXgriddim/2)+1:intXgriddim
                    matVOF(k,i,j)=matCorrectVOF(l);
                    l=l+1;
                end
            end
        end
    elseif p==6
        for k=intZgriddim:-1:(intZgriddim/2)+1
            for j=(intYgriddim/2)+1:intYgriddim
                for i=(intXgriddim/2)+1:intXgriddim
                    matVOF(k,i,j)=matCorrectVOF(l);
                    l=l+1;
                end
            end
        end
    elseif p==7
        for k=(intZgriddim/2):-1:1
            for j=(intYgriddim/2)+1:intYgriddim
                for i=(intXgriddim/2)+1:intXgriddim
                    matVOF(k,i,j)=matCorrectVOF(l);
                    l=l+1;
                end
            end
        end
    end
    
end

A=matVOF;

end

