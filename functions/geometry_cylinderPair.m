% Function to initialise the geometry of the appendages.

function [stks_cyl] = geometry_cylinderPair(rho2,dsep2,psi2,PRAx2,PRAy2,mod1)
    Npts = floor(rho2*2*pi)+1; % Get the number of stokeslets on each cylinder
    theta = linspace(0,2*pi,Npts); theta = theta(1:end-1); % Get the angles to parameterise the surfaces.
    circle = [0.95*cos(theta);0.95*sin(theta);zeros(floor(rho2*2*pi),1)']'; % This is the co-ords of one surface of an appendage.
    stks_cyl1 = circle + [PRAx2+(1+dsep2/2)*cos(psi2),PRAy2+(1+dsep2/2)*sin(psi2),4+2*mod1]; % Shift the first circle to + side of the separation.
    stks_cyl2 = circle + [PRAx2-(1+dsep2/2)*cos(psi2),PRAy2-(1+dsep2/2)*sin(psi2),5+2*mod1]; % Shift the second circle to - side of the separation.
    stks_cyl  = [stks_cyl1;stks_cyl2]; % Combine both coordinate sets into one array for output.
end

