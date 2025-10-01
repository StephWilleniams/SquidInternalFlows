% Function to solve the inverse problem to get the stokeslet forces.
% Here, the location of the boundary conditions is NOT localized to the position of the stokeslet

function [iS] = getForces(stks1,eps_reg1,mu1)

    %% Create the stokeslet system

    ntot = length(stks1(:,1)); % Get the number of Stokeslets
    S = zeros(2*ntot); % Preallocate the full linear stokeslet.

    % Loop pairwise
    for ii = 1:ntot
        for jj = 1:ntot

            dx = stks1(ii,1)-stks1(jj,1); % x-distance
            dy = stks1(ii,2)-stks1(jj,2); % y-distance
            R = sqrt(dx^2 + dy^2 + eps_reg1^2) + eps_reg1; % Regularized distance
            rho = (R+eps_reg1)/(R*(R-eps_reg1)); % Used to simplify Stokeslet calculation

            Sub = zeros(2,2); % Preallocate Stokeslet from force jj on position ii

            %Calculate the Stokeslet
            Sub(1,1) = -log(R) + eps_reg1*rho ...
                     + dx^2  * rho/R; % Effect of x-component of ii on flow x-component as position jj
            
            Sub(1,2) = dx*dy * rho/R; % Effect of y-component of ii on flow x-component as position jj
            
            Sub(2,1) =  dx*dy * rho/R; % Effect of y-component of ii on flow y-component as position jj
            
            Sub(2,2) = -log(R) + eps_reg1*rho ...
                     + dy^2  * rho/R; % Effect of y-component of ii on flow y-component as position jj
            
            S(2*ii-1:2*ii,2*jj-1:2*jj) = Sub; % Store the position in the full system Stokeslet

        end % End jj loop
    end % End ii loop

    %% Explicit viscosity term

    S = S/(4*pi*mu1);

    %% Add in the force constraints columns.

    IN2 = zeros(2*ntot,2);
    IN2(1:2:end,1) = -1; IN2(2:2:end,2) = -1;
    S = [S,IN2];
    I2N = [-IN2',zeros(2,2)]; S = [S;I2N];

    %% Add in the torque constraints column.

    Rcol = zeros(2*ntot,1);
    Rcol(1:2:end) = -stks1(:,2);
    Rcol(2:2:end) =  stks1(:,1);
    S    = [S,[Rcol;0;0]];
    Rcol = [-Rcol',0,0,0];
    S    = [S;Rcol];

    %% Invert the linear system.

    iS = inv(S);

end

