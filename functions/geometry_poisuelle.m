% Function to initialise the geometry of the channel with a contraction.
% Outputs an Nx3 array with N (x,y) positisons for each of the stokeslets, with the third argument a one-hot encoded boundary condition.
% boundary codes: (:,3) = 1, no-slip/pen. (:,3) = 2, Poisuelle. (:,4) = 3, free.

function [channel_stks] = geometry_poisuelle(rho2,Lt2,Lm2,Lb2,theta2,Ptx2,Pty2)
    
    removeNum = 1;

    % Construct the right boundary
    WRTl = linspace(0,Lt2,floor(Lt2*rho2)); % Parameterise the right wall top section.
    WRMl = linspace(0,Lm2,floor(Lm2*rho2)); % Parameterise the right wall middle section.
    WRBl = linspace(0,Lb2,floor(Lb2*rho2)); % Parameterise the right wall bottom section.

    % Translate them into the right places in space
    WRT = [Ptx2 - WRTl*sin(0); Pty2 - WRTl*cos(0);ones(floor(Lt2*rho2),1)']'; % Get the top segment coordinates.
    WRM = [Ptx2 - WRMl*sin(theta2); Pty2 - WRTl(end)*cos(0) - WRMl*cos(theta2);ones(floor(Lm2*rho2),1)']'; % Get the middle segment coordinates.
    WRB = [Ptx2 - WRMl(end)*sin(theta2) - WRBl*sin(0);Pty2 - WRTl(end)*cos(0) - WRMl(end)*cos(theta2) - WRBl*cos(0);ones(floor(Lb2*rho2),1)']'; % Get the bottom segment coordinates.
    
    % Combine them into one array
    wallR = [WRT;WRM(1+removeNum:end-removeNum,:);WRB]; % Combine all the segments to give the right wall.

    % Construct the left boundary.
    wallL = wallR; % Copy the left boundary coordinates.
    wallL(:,1) = -wallL(:,1); % Reflect x-coordinate.

    % Construct the top boundary.
    wallT = [linspace(-Ptx2,Ptx2,floor(2*Ptx2*rho2)); Pty2*ones(floor(2*Ptx2*rho2),1)'; 2*ones(floor(2*Ptx2*rho2),1)']';
    wallT = wallT(2:end-1,:);

    % Combine to give the full position array.
    temp_channel_stks = [wallL;wallR;wallT];

    % Set the loction of the stokeslets to be the same as where their boundary conditions are.
    channel_stks = [temp_channel_stks(:,1:2),temp_channel_stks(:,3)];

end

