function run_danio_feed(rerun_flow,show_sim)
% Wrapper to run a single simulation of danio predator-prey interaction


if nargin < 2
    show_sim = 0;
end


%% Define paths

% Root directory for saving and loading simulation data
% (set up to deal with Matt's 2 computers)
if isdir('/Users/mmchenry/Dropbox/Projects/Holzmann/sims/danio')
    sim_path = '/Users/mmchenry/Dropbox/Projects/Holzmann/sims/danio';
else
    sim_path = '/Volumes/Flow HD/Dropbox/Projects/Holzmann/sims/danio';
end


%% Set parameter values

% Default parameter values
[sim,pred,prey] = params_danio;


%% Generate (or load) flow field

% Prompt to create new flow field, if necessary
if 1 % ~isempty(dir([sim_path filesep 'danio_flow.mat']));
    
    if nargin < 1
        answer = questdlg('Generate new flow field data?',...
                          'Flow data',...
                          'No','Yes','Cancel','No');
    elseif rerun_flow
        answer = 'Yes';
    else
        answer = 'No';
    end
                  
    if strcmp(answer,'Cancel')
        return
        
    elseif strcmp(answer,'Yes')
        % Load flow data (x2, y2, Field1)
        load([sim_path filesep 'Speed_Field.mat'])
        
        % Create % save flow field structure
        fl = flow_daniofeed(sim,pred,x2,y2,Field1);
        save([sim_path filesep 'danio_flow'],'fl')
    end
    
end

% Load 'fl' flow field structure
load([sim_path filesep 'danio_flow.mat'])  

% Animate flow field
%vis_flow(fl,sim)


%% Run simulation

% Run simulation
tic
r = siffer(sim,prey,fl);

% Report time
t_lapse = toc;
disp(' ')
disp(['Sim completed in ' num2str(t_lapse) 's'])

if show_sim
    % Animate simulation results
    vis_sim(fl,r)
end

% Visualize results
figure;

subplot(2,1,1)
plot(1000.*r.t,1000.*r.pos(1,:),'k-',1000.*r.t,1000.*r.pred_pos(1,:),'r--')
xlabel('time (ms)');ylabel('position (mm)');
legend('prey','pred','Location','West')

if ~isnan(r.cap) && r.cap
    title('Captured!');
else
    title('Escaped!');
end

subplot(2,1,2)
plot(1000.*r.t,r.D(1,:),'r-',1000.*r.t,r.PF(1,:),'b-', ...
    1000.*r.t,r.AR(1,:),'g-',1000.*r.t,r.Th(1,:),'m-')
legend('D','PF','AR','Th','Location','SouthWest')
xlabel('time (ms)');ylabel('Force in x-direction (N)');


end

