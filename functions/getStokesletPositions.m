% Function to initialise the geometry of the system.
% Outputs an Nx2 array with N (x,y) positisons for each of the stokeslets.

function stks = getStokesletPositions(rho1,geometry,U01,UC1)

% Channel parameters
Lt    = geometry.channel_parameters(1);
Lm    = geometry.channel_parameters(2);
Lb    = geometry.channel_parameters(3);
theta = geometry.channel_parameters(4);
Ptx   = geometry.channel_parameters(6);
Pty   = geometry.channel_parameters(7);

% Appendage parameters
dsep = geometry.appendage_parameters(1);
psi  = geometry.appendage_parameters(2);
PRAx = geometry.appendage_parameters(3);
PRAy = geometry.appendage_parameters(4);

%% Set the positions of the stokeslets

stks_channel = geometry_poisuelle(rho1,Lt,Lm,Lb,theta,Ptx,Pty); % Set the channel geometry.
stks_appendages1 = geometry_cylinderPair(rho1,dsep,psi,PRAx,PRAy,0); % Set the left appendage pair geometry.
stks = [stks_channel;stks_appendages1]; % Combine all structures.

%% Set the corresponding boundary velocities

[nStok,~] = size(stks); % Get number of stokeslets required at given linear density.
BdryVelo = zeros(nStok,2); % Preallocate an array to store the boundary velocities (for BVP).

%% No-slip boundaries -- stks(:,3) == 1, this code currently does nothing, so is commented out.
%ind = find(stks(:,5)==1);
%BdryVelo(ind,:) = 0; % Set zero-flow on the channel boundaries

%% Poisuelle boundary flow -- stks(:,3) == 2
% Flow going as 1-(r/a)^2, r = distance from channel center, a = half channel-width.
ind = find(stks(:,3)==2); % Get relevant Stokeslets for which this boundary condition holds
nTemp = length(ind); % Get the number of Stokeslets contained in this set.
BdryVelo(ind,:) = poisuelleFlow(nTemp,U01); % Prescibe the boundary velocity at these points.

%% No-slip boundaries -- stks(:,3) == 3
% Can add other boundary conditions for the bottom of a channel here.

%% No-slip boundaries -- stks(:,3) == 4,5,6,7
% Boundary conditions following those found in Nawroth 2017, see SM.

ind = find(stks(:,3)==4); nTemp = length(ind);
%BdryVelo(ind,:) = surfaceFlow(nTemp,2*pi*(358)/360,UC1);
BdryVelo(ind,:) = surfaceFlow(nTemp,2*pi*(358)/360,UC1);

ind = find(stks(:,3)==5); nTemp = length(ind);
%BdryVelo(ind,:) = surfaceFlow(nTemp,2*pi*(88)/360,UC1);
BdryVelo(ind,:) = surfaceFlow(nTemp,2*pi*(88)/360,UC1);

ind = find(stks(:,3)==6); nTemp = length(ind);
BdryVelo(ind,:) = surfaceFlow(nTemp,2*pi*(182)/360,UC1);

ind = find(stks(:,3)==7); nTemp = length(ind);
BdryVelo(ind,:) = surfaceFlow(nTemp,2*pi*(92)/360,UC1);

%% Append the BdryVelo array to stks
stks = [stks,BdryVelo];

temp = stks(stks(:,3)==4,:);
temp = [temp;stks(stks(:,3)==5,:)];

temp(:,3) = temp(:,3) + 2;

temp(:,1) = -temp(:,1);
temp(:,4) = -temp(:,4);

stks = [stks;temp];

%% Safety check, remove repeating points

for i = 1:length(stks(:,1))
    for j = 1:length(stks(:,1))
        d = norm(stks(i,1:2)-stks(j,1:2));
        if (d == 0 && i~=j)
            disp('Removal active -- check geometry')
            stks(j,:) = [];
        end
    end
end

end