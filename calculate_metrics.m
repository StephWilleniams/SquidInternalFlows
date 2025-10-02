
% Load in the data

all_data_x = [];
all_data_y = [];

for i = 1:11

    load(['outputs/data_' num2str(i) '.mat']);

    all_data_x = [all_data_x;outputs_x];
    all_data_y = [all_data_y;outputs_y];

end

%%

parameters

% Get the appendage centers
AP1 = [system.appendage_parameters(3) + (1+system.appendage_parameters(1)/2)*cos(system.appendage_parameters(2)), ...
       system.appendage_parameters(4) + (1+system.appendage_parameters(1)/2)*sin(system.appendage_parameters(2))];
AP2 = [system.appendage_parameters(3) - (1+system.appendage_parameters(1)/2)*cos(system.appendage_parameters(2)), ... 
       system.appendage_parameters(4) - (1+system.appendage_parameters(1)/2)*sin(system.appendage_parameters(2))];
AP3 = [-1,1].*AP1; AP4 = [-1,1].*AP2;

[si1,si2,si3] = size(all_data_x);

dist = zeros(si1,si3,si2);

runs = 1:si1;
N = 1:si3;
times = 1:si2;

threshold_distance = 2.5;

for run = runs
    for n = N
        for time = times

            % (X) Is this being implemented right???
            pos = [all_data_x(run,time,n),all_data_y(run,time,n)];
            d1 = norm(AP1 - pos);
            d2 = norm(AP2 - pos);
            d3 = norm(AP3 - pos);
            d4 = norm(AP4 - pos);
            dist(run,n,time) = min([d1,d2,d3,d4]);
        end
    end    
end

save('outputs/dist','dist','-V7.3')

load('outputs/dist')

params = readmatrix('inouts/params.txt');

dist = dist(1:length(params(:,1)),:,:);

[si1,si2,si3] = size(dist);

runs = 1:si1;
pIndicies = 1:si2;
times = 1:si3;

dTresh = 2.5;
close_enough = zeros(si1,si2,si3);

for run = runs
    for partIndex = pIndicies
        for t = times
            if dist(run,partIndex,t) <= dTresh
                close_enough(run,partIndex,t) = 1;
            end
        end
    end
end

%%
% 
% data = reshape(dist(1,:,:),[si2 si3]);
% 
% imagesc(data); hold on

%%
% 
% data = reshape(close_enough(end,:,:),[si2 si3]);
% 
% imagesc(data,'alphaData',0.5)

%%

output1 = zeros(length(runs),1);
output2 = zeros(length(runs),1);

for run = runs
    temp_output1 = 0;
    temp_output2 = 0;
    for partIndex = pIndicies
        indicator = reshape(close_enough(run,partIndex,:),[length(close_enough(run,partIndex,:)) 1]);
        temp_output1 = temp_output1 + sum(indicator);
        ch_len = 0;ma_len = 0;
        for jj = 1:length(indicator)
            if indicator(jj) == 1
                ch_len = ch_len + 1;
            else
                ma_len = max(ma_len,ch_len);
                ch_len = 0;
            end
        end
        temp_output2 = max(temp_output2,ma_len);
    end
    output1(run) = temp_output1;
    output2(run) = temp_output2; 
end

dt = 30/300;
output1 = output1 * dt;
output2 = output2 * dt;

if 1

    % Total time
    %writematrix(scaling,'data_outputs/scalings.txt');
    % Open the text file for writing
    fid = fopen('outputs/outputs1.txt', 'w');
    % Write the variable to the file
    fprintf(fid, '%g\n', output1);
    % Close the file
    fclose(fid);

    % Max time
    % Open the text file for writing
    fid = fopen('outputs/outputs2.txt', 'w');
    % Write the variable to the file
    fprintf(fid, '%g\n', output2);
    % Close the file
    fclose(fid);

end
