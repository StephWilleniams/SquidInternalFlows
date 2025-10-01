% Function to solve the prescribe a Poisuelle flow across the channel entry.

function [Uflow] = poisuelleFlow(nTemp1,U02)

if mod(nTemp1,2) == 1 % Odd sized array
    N = floor(nTemp1/2);
    Uflow = zeros(nTemp1,2); % Preallocate an output array
    Uflow(1:N,2) = U02*(1 - (nTemp1/2:-1:1).^2/((nTemp1/2)^2)); % Calculate the flow profile across half the channel
    Uflow(N+1,2) = U02;
    Uflow(2+N:end,2) = flip(Uflow(1:N,2)); % Flip the first half
    %disp('odd-array')
else % Even sized array
    Uflow = zeros(nTemp1,2); % Preallocate an output array
    Uflow(1:nTemp1/2,2) = U02*(1 - (nTemp1/2:-1:1).^2/((nTemp1/2)^2)); % Calculate the flow profile across half the channel
    Uflow(1+nTemp1/2:end,2) = flip(Uflow(1:nTemp1/2,2)); % Flip the first half
    %disp('even-array')
end

end

