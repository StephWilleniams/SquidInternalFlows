% Function to calculate the surface flow on the appendages.

function [Uflow] = surfaceFlow(nTemp1,rot1,Uc1)

scaler = Uc1; % Get the maximum strength of the appendage flow

theta = linspace(2*pi/nTemp1,2*pi,nTemp1); % Parameterise the surface of the appendage on which the flow is applied
tangent = [-sin(theta);cos(theta)]'; % Calculate the tangent vectors along the surface
magnitude = zeros(nTemp1,1); % Preallocate space for the magnitude scaler
rot1 = mod(rot1-3*pi/4, 2*pi); % Get the value of the rotation within a range [0,2*pi]

if rot1 <= pi/2 && rot1 ~= 0 % Case 1 -- rotation by < 90deg

    a = find(theta > rot1);
    b = find(theta(a) < rot1+3*pi/2); b = a(b);
    
    angles = 4*(theta(b)-rot1)/3;
    magnitude(b) = -sin(angles);

elseif rot1 >= pi/2 % Case 2 -- rotation by > 90deg (handles array overlaps)

    a = find(theta < rot1 - pi/2);
    angles = 4*(2*pi + theta(a) - rot1)/3;
    magnitude(a) = -sin(angles);

    b = find(theta > rot1);
    angles = 4*(theta(b) - rot1)/3 ;
    magnitude(b) = -sin(angles);

elseif rot1 == pi/2 % Case 3 -- rotation by 90deg (handles array overlaps)

    a = find(theta > pi/2);

    angles = 4*(theta(a)-rot1)/3;
    magnitude(a) = -sin(angles);

else % Case 3 -- rotation by 0deg (handles array overlaps)

    a = find(theta < 3*pi/2);

    angles = 4*theta(a)/3;
    magnitude(a) = -sin(angles);

end

% Combine all calculations to give the desired output flows
Uflow = zeros(nTemp1,2); % Preallocate output array
Uflow(:,1) = scaler*magnitude.*tangent(:,1);
Uflow(:,2) = scaler*magnitude.*tangent(:,2);

end

