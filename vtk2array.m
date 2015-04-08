function A = vtk2array(folderloc, time, numprocesses)
%vtk2array
%   Returns a three dimensional array of the color function at a specific
%   point in time.  To use, simply input the VTK folder, desired
%   timestep, and number of processes.  This should work for variables other than the color
%   function, but only if the value is described at the center of the grid
%   cell.  Please feel free to contact me with any comments/concerns.
%   
%   Note: numprocesses MUST be a perfect cube!
%
%   Sample Call:
%   folderloc = '/home/jakedynes/paris-stable-8proc-ex/Tests/VOF/out/VTK/';
%   time = '00000';
%   numprocesses = 8;
%   A = vtk2array(folderloc, time, numprocesses)
%   
%   Jake Dynes, March 2015 (dynesj1@montclair.edu)
%
%   GPL Licence
%
%     This file is part of PARIS.
%
%     PARIS is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     PARIS is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with PARIS.  If not, see <http://www.gnu.org/licenses/>.
%
    
for p=0:numprocesses-1
    % import data from the current processor
    if p<10
        intProcessor=strcat('0000',int2str(p));
    elseif p<100 && p>9
        intProcessor=strcat('000',int2str(p));
    elseif p<1000 && p>99
        intProcessor=strcat('00',int2str(p));
    elseif p<10000 && p>999
        intProcessor=strcat('0',int2str(p));
    elseif p<100000 && p>9999
        intProcessor=strcat('',int2str(p));
    end
    strFileloc=strcat(folderloc, 'VOF', time, '-', intProcessor, '.vtk');
    matData=importdata(strFileloc);
    proccube=nthroot(numprocesses,3);

    % some setup from first processor
    if p==0

        % get grid size (specified in 'input')
        intXgriddim=(str2double(matData.textdata{5,2})-1)*proccube;
        intYgriddim=(str2double(matData.textdata{5,3})-1)*proccube;
        intZgriddim=(str2double(matData.textdata{5,4})-1)*proccube;
        procXsize=intXgriddim/proccube;
        procYsize=intYgriddim/proccube;
        procZsize=intZgriddim/proccube;

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
        xDomain=intXgriddim*((matPoints(1,1)*2)/proccube);
        yDomain=intYgriddim*((matPoints(1,2)*2)/proccube);
        zDomain=intZgriddim*((matPoints(1,3)*2)/proccube);
    
    end
    
    % sort out ghost values
    intNumInVOF=(intXgriddim/proccube)*(intYgriddim/proccube)*(intZgriddim/proccube);
    matCorrectVOF=zeros(intNumInVOF,1);
    xzero=(p-mod(p,proccube^2))/proccube^2;
    yzero=floor(mod(p,proccube^2)/proccube);
    zzero=mod(p,proccube);
    j=1;
    
    for i=1:intTotalCellswithGhost
        xOK=0;
        yOK=0;
        zOK=0;
        
        if and(matPoints(i,1)>xzero*xDomain,matPoints(i,1)<(xzero+1)*xDomain)
            xOK=1;
        end
        if and(matPoints(i,2)>yzero*yDomain,matPoints(i,2)<(yzero+1)*yDomain)
            yOK=1;
        end
        %if and(matPoints(i,3)>0,matPoints(i,3)<zDomain)
        if and(matPoints(i,3)>zzero*zDomain,matPoints(i,3)<(zzero+1)*zDomain)
            zOK=1;
        end
               
            A=[xOK,yOK,zOK];
            if all(A)
                matCorrectVOF(j)=matValues(i,1);
                j=j+1;            
            end  
    end
    
    %create vof array for current processor
    matprocVOF=zeros(intZgriddim/proccube, intXgriddim/proccube, intYgriddim/proccube);
    l=1;

    for c=intZgriddim/proccube:-1:1
        for b=1:intYgriddim/proccube
            for a=1:intXgriddim/proccube
                matprocVOF(c,a,b)=matCorrectVOF(l);
                l=l+1;
            end
        end
    end
    
    % stick the processor VOF array in the right spot in the overall VOF
    % array
    zzero=proccube-mod(p,proccube);
    yzero=floor(mod(p,proccube^2)/proccube);
    xzero=floor(p/proccube^2);
    
    matVOF(procZsize*(zzero-1)+1:procZsize*zzero,...
        (xzero*procXsize)+1:(xzero+1)*procXsize,...
        (yzero*procYsize)+1:(yzero+1)*procYsize)=matprocVOF;
end

A=matVOF;

end
