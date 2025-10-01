% Title: Code to get the fluid flow from a set of Regularized Stokeslets.
% This flow is applied to a time-dependent flow for some particles to get the trajectories.
% Those trajectories spending sufficiently long near the appendages are "captured".
% Author: Stephen Williams.

function CalculateTrajectory(index)

% Get a parallel pool for workers
if isempty(gcp('nocreate'))
    parpool('Threads')
end
% Get the number of available cores in the pool
numCores = feature('NumCores'); 
% Error if there aren't enough available cores - slurm can under-allocate at this step
if numCores < 8
    error('Not enough cores to complete execution.')
end

%% Add the function files need to run
addpath('functions/')
addpath('classes/')

%% Set up the timer to report runtime
tic;

%% Set parameters
% Parameters:
NDL = 45; % Non-dimensional length scale
rho = 10;
eps_reg = 0.5/rho; % Regularisaton parameter.
mu = 8.9E-3; % Viscocity of water
U0 = -100/45; UC0 = 600/45; % Required for default set up

% Set the system geometry
system = geometry; % Preallocate a channel geometry
system.channel_parameters(1:4) = [(367+500)/NDL,366/NDL,200/NDL,200/NDL];
system.channel_parameters(5) = system.channel_parameters(1)+sin(pi/4)*system.channel_parameters(2)+system.channel_parameters(3); % Total height of system simulated.
system.channel_parameters(6) = 500/NDL; 
system.channel_parameters(7) = system.channel_parameters(1)+system.channel_parameters(2)/2; % Position y of top point of right boundary.
system.appendage_parameters(1:4) = [0.95*90/NDL,0.82-pi/2,125/NDL,11];

%% Set channel geometry
stks = getStokesletPositions(rho,system,U0,UC0);

%% Solve for the forces
[iS] = getForces(stks,eps_reg,mu); % This could be done with a pre-load, for a minor speed up.

%% Simulate the particle motion

pfnum = 24; % Total number of loops to run for each thread
taskPerCore = ceil(pfnum/numCores); % Expected tasks to perform per thread
LHS = false; % Optional: generate LHS locally for the thread
if LHS
    nvars = 4; % Number of variables in the LHS
    scaling = 2*(lhsdesign(pfnum,nvars)-0.5); % Get the parameter scaling
else % Preload the parameter set, used for Sobol analysis
    scaling = readmatrix('inputs/params_2pow12.txt'); %#ok<UNRCH> % Get the parameter scaling
    maxIndex_scaling = length(scaling(:,1));
    minRange = 1+(index-1)*pfnum;
    maxRange = min(index*pfnum,maxIndex_scaling);
    scaling = scaling(minRange:maxRange,:);
end
nparticles = 50; % Number of particle trajectories to simulate
tmin = 0; tmax = 30; ntsteps = 301; % Time paramaters
delta = 20; % Boundary force potential scaling

values_x = cell(numCores);
values_y = cell(numCores);
for core = 1:numCores
    values_x{core} = zeros(taskPerCore,ntsteps,nparticles);
    values_y{core} = zeros(taskPerCore,ntsteps,nparticles);
end

%%

parfor core = 1:numCores

    localScaling = scaling; % Lower communication to overhead with a local copy.
    lowerPF = 1 + taskPerCore*(core-1);
    upperPF = taskPerCore*(core);

    tmp_x_pos = zeros(taskPerCore,ntsteps,nparticles);
    tmp_y_pos = zeros(taskPerCore,ntsteps,nparticles);

    for pfind = lowerPF:upperPF

        local_index = 1 + pfind - lowerPF;

        % Set up the current iteration of the systen
        sc = localScaling(pfind,:); % Extract the values the parameters need to be scaled by for this loop
        U0    = sc(1)*75/NDL + 125/NDL; % Ventilation flow parameter
        omega = sc(2)*2.5 + 3; % Pouiselle flow parameter
        UC    = sc(3) + 600/NDL; % Cilia flow strength parameter
        y0    = sc(4)*2 + 16; % Initial distance upstream from the appendages: appendages at 11 N.D. units.

        % Get the solver conditions
        ic = [linspace(-10,0,nparticles);y0*ones(1,nparticles)]; % Initial conditions
        tempstks = stks; % Get a local copy of the stks array to be updated within this thread of the code

        % Modify the U_c values
        ind = find(stks(:,3)==4); nTemp = length(ind);
        tempstks(ind,4:5) = surfaceFlow(nTemp,2*pi*(358)/360,UC);
        ind = find(stks(:,3)==5); nTemp = length(ind);
        tempstks(ind,4:5) = surfaceFlow(nTemp,2*pi*(88)/360,UC);
        tempstks(stks(:,3)==6,4:5) = [-1,1].*tempstks(stks(:,3)==4,4:5);
        tempstks(stks(:,3)==7,4:5) = [-1,1].*tempstks(stks(:,3)==5,4:5);

        % Solve for the time series of the particle motion
        [~,y] = ode45(@(t,y) flowdefunct(t,y,tempstks,iS,eps_reg,omega,U0,delta,system,mu),linspace(tmin,tmax,ntsteps),ic); % Solve the trajectories under the flow

        % Extract and store solution
        x_pos = y(:,1:2:end); y_pos = y(:,2:2:end); % Extract the positions into an array for plotting

        tmp_x_pos(local_index,:,:) = x_pos;
        tmp_y_pos(local_index,:,:) = y_pos;
        
    end

    values_x{core} = tmp_x_pos;
    values_y{core} = tmp_y_pos;

end

% Unwrap the parallel storage for outputting
outputs_x = [];
outputs_y = [];
for i = 1:numCores
    outputs_x = [outputs_x;values_x{i}];
    outputs_y = [outputs_y;values_y{i}];
end

% Finalise and output
save(['outputs/data_' num2str(index)],'scaling','outputs_x','outputs_y')
runtimecount = toc; disp(['Total runtime: ' num2str(runtimecount) ' seconds.'])

end

%% Functions needed for ODE45

function dydt = flowdefunct(t1,y1,stks1,iS1,eps_reg1,omega,U0,delta1,system1,mu1)

R = 1; % Radius of appendages
r = 0.5/45; % Radius of particles

a = find(stks1(:,3) == 2); % Find the Pousielle boundary sections
Ut = -U0*(1+cos(t1*omega))/2; % Set the maximum flow rate in
stks1(a,4:5) = poisuelleFlow(length(a),Ut); % Re-set the Pousielle boundary sections

Uflow = calculateFlowVector(stks1,iS1,y1,eps_reg1,mu1); % Get the flow from fluid interactions

% Get the appendage centers
AP1 = [system1.appendage_parameters(3) + (1+system1.appendage_parameters(1)/2)*cos(system1.appendage_parameters(2)), ...
    system1.appendage_parameters(4) + (1+system1.appendage_parameters(1)/2)*sin(system1.appendage_parameters(2))];
AP2 = [system1.appendage_parameters(3) - (1+system1.appendage_parameters(1)/2)*cos(system1.appendage_parameters(2)), ...
    system1.appendage_parameters(4) - (1+system1.appendage_parameters(1)/2)*sin(system1.appendage_parameters(2))];
AP3 = [-1,1].*AP1; AP4 = [-1,1].*AP2;
app_pos = [AP1;AP2;AP3;AP4]; % Collect up all of the the positions
x = [y1(1:2:end)';y1(2:2:end)'];
[~,I] = min([vecnorm(x(1:2,:)-AP1',2) ; ...
    vecnorm(x(1:2,:)-AP2',2) ; ...
    vecnorm(x(1:2,:)-AP3',2) ; ...
    vecnorm(x(1:2,:)-AP4',2) ]);

% Get the perturbation from solid interactions
dU =  delta1 * ((R + r) > vecnorm(x - app_pos(I,:)')) .* ...
    ((R + r)  - vecnorm(x - app_pos(I,:)')) .* ...
    (x - app_pos(I,:)')  ./ vecnorm(x - app_pos(I,:)');

% Get the motion resulting from these
dydt = zeros(length(y1),1);
dydt(1:2:end) = Uflow(:,1) + dU(1,:)';
dydt(2:2:end) = Uflow(:,2) + dU(2,:)';

end