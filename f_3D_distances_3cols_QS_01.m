function [D, ID] = f_3D_distances_3cols_QS_01(pixsize, zsize, coord1, coord2, coord3)
%   Calculates 3D distances between sets of 3D coordinates

%   INPUT
%   pixsize: pixel size
%   zsize: z-step size
%   coord1: [x, y, z] coordinates of FISH channel 1
%   coord2: [x, y, z] coordinates of FISH channel 2
%   coord3: [x, y, z] coordinates of FISH channel 3 (optional)

%   OUTPUT
%   D: Distance between mutual nearest neighbors (D(1-2) if 2 channels
%   or [D(1-2), D(1-3), D(2-3)] if 3 channels)
%   IDs: FISH IDs ([1, 2] if 2 channels or [1, 2, 3] if 3 channels)


% Calculates all pairwise distances between channels 1 & 2
alldistances12 = [];
for icoord1 = 1:size(coord1,1)
	for icoord2 = 1:size(coord2,1)
         xdist = (coord2(icoord2,1) - coord1(icoord1,1)) * pixsize;
         ydist = (coord2(icoord2,2) - coord1(icoord1,2)) * pixsize;
         zdist = (coord2(icoord2,3) - coord1(icoord1,3)) * zsize;     
         distance = sqrt(xdist^2 + ydist^2 + zdist^2);
         distances = [distance, icoord1, icoord2];
         alldistances12 = [alldistances12; distances];
	end
end      
    
% If 3 channels, calculates all pairwise distances between channels 1 & 3 and 2 & 3
if nargin == 5
    alldistances13 = [];
    for icoord1 = 1:size(coord1,1)
        for icoord3 = 1:size(coord3,1) 
              xdist = (coord3(icoord3,1) - coord1(icoord1,1)) * pixsize;
              ydist = (coord3(icoord3,2) - coord1(icoord1,2)) * pixsize;
              zdist = (coord3(icoord3,3) - coord1(icoord1,3)) * zsize; 
              distance = sqrt(xdist^2 + ydist^2 + zdist^2);
              distances = [distance, icoord1, icoord3];
              alldistances13 = [alldistances13; distances];

         end
    end
    alldistances23 = [];
    for icoord2 = 1:size(coord2,1)
         for icoord3 = 1:size(coord3,1)
              xdist = (coord3(icoord3,1) - coord2(icoord2,1)) * pixsize;
              ydist = (coord3(icoord3,2) - coord2(icoord2,2)) * pixsize;
              zdist = (coord3(icoord3,3) - coord2(icoord2,3)) * zsize;       
              distance = sqrt(xdist^2 + ydist^2 + zdist^2);
              distances = [distance, icoord2, icoord3];
              alldistances23 = [alldistances23; distances];    
         end
    end
end

% Creates structure with all distances for each channel pair
sdistances = struct;
sdistances(1).dist = alldistances12;
if nargin == 5
    sdistances(2).dist = alldistances13;
    sdistances(3).dist = alldistances23;
end

% For each channel pair, calculates mutual nearest neighbor
for j = 1:size(sdistances,2)     
     listdistances = sdistances(j).dist;
     sdistances(j).mutdist = [];
     ncoord1 = max(listdistances(:,2));
     for ID1 = 1:ncoord1
          distancesID1 = listdistances(listdistances(:,2) == ID1,:);   
          mindisID1 = distancesID1(distancesID1(:,1) == min(distancesID1(:,1)),:);
          ID2 = mindisID1(1,3);
          distancesID2 = listdistances(listdistances(:,3) == ID2,:);
          mindisID2 = distancesID2(distancesID2(:,1) == min(distancesID2(:,1)),:);
          ID1b = mindisID2(1,2);
          if ID1 == ID1b
              mutnear = [mindisID1(1,1), ID1, ID2];
              sdistances(j).mutdist = double([sdistances(j).mutdist; mutnear]);
          end
     end
end

% Output if 2 channels
if nargin == 4
    D = sdistances(j).mutdist(:,1);
    ID = sdistances(j).mutdist(:,2:3);
else

% Triple mutual nearest neighbors if 3 channels
idx12common1 = ismember(sdistances(1).mutdist(:,2), sdistances(2).mutdist(:,2)); % Common 1 in 12 & 13
precleaned12 = sdistances(1).mutdist(idx12common1,:); % Precleaned 12
idx23common2 = ismember(sdistances(3).mutdist(:,2), precleaned12(:,3)); % Common 2 in 12 and 23
precleaned23 = sdistances(3).mutdist(idx23common2,:); % Precleaned 23
idx13common3 = ismember(sdistances(2).mutdist(:,3), precleaned23(:,3)); % Common 3 in 13 and 23
cleaned13 =  sdistances(2).mutdist(idx13common3,:); % Fully cleaned 13
idx12 = ismember(precleaned12(:,2), cleaned13(:,2));
cleaned12 = precleaned12(idx12, :);
idx23 = ismember(precleaned23(:,2), cleaned12(:,3));
cleaned23 = precleaned23(idx23,:);
    
% Output if 3 channels 
D = [cleaned12(:,1), cleaned13(:,1), cleaned23(:,1)]; 
ID = [cleaned13(:,2), cleaned12(:,3), cleaned13(:,3)];

end
end


