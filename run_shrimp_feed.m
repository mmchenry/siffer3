function run_shrimp_feed
% Wrapper to run a single simulation of predation on a shrimp (same
% conditions as V0.9.5 simulation.

%% Define paths

% Root directory for saving and loading simulation data
sim_path = '/Users/mmchenry/Dropbox/Projects/Holzmann/sims';

% Path to CFD data
cfd_path = 'flow_field.mat';


%% Set parameter values

% Default parameter values
[sim,pred,prey] = default_params;

% Visualize flow
vis_flow = 0;

% Vizulaize simulation results
vis_sim = 1;


%% Calculate (or load) predator's flow field

% Flow field is loaded, if previously saved
if ~isempty(dir([sim_path filesep 'feed_field.mat']))
    
    % Echo data loading
    disp(' ');disp('Loading flow data . . .');disp(' ');
    
    % Load 'fl' structure (flow data)
    load([sim_path filesep 'feed_field.mat']);
    
% Otherwise, calculate the flow field
else
    
    % Load CFD flow velocity data ('f')
    load(cfd_path)
    
    % Create feeding field
    fl = feedField(sim,pred,f,vis_flow);
    
    % Save feeding flow field data
    save([sim_path filesep 'feed_field.mat'],'fl')
    
end


%% Run simulation
tic

% Run simulation
r = siffer(sim,pred,prey,fl);

% Save feeding flow field data
save([sim_path filesep 'sim_results.mat'],'r')

% Report time
t_lapse = toc;
disp(['Sim complete in ' num2str(t_lapse) 's'])


%% Visualize results

if vis_sim
    
figure

subplot(3,1,1)

ax1 = gca;
h1 = line(fl.t,1e3*fl.pos(:,1),'Color','b','parent',ax1,'linestyle','--');
set(ax1,'XColor','b','YColor','b')
xlabel('time(s)');
ylabel('Predator pos. (mm)');
ax2 = axes('Position',get(ax1,'Position'),...
          'XAxisLocation','top',...
          'YAxisLocation','right',...
          'Color','none',...
          'YColor','black');
%set(ax2,'YColor','black')
ylabel('Speed(m/s)');

h2 = line(fl.t,fl.gape,'Color','black','parent',ax2);
h3 = line(fl.t,1e3*fl.gape_spd,'Color','b');
h = [h1(1) h2 h3];
legend(h,'x pos (mm)','gape(mm)','x spd (m/s)');
title('Kinematic characteristics');

subplot(3,1,2)
%plot(r.t,r.drag(:,1)+r.PF(:,1),'k-')
xlabel('time (s)')
ylabel('Total x force (N)')

subplot(3,1,3)
plot(1000.*fl.pos(:,1),1000.*fl.gape./2,'r-',...
    1000.*fl.pos(:,1),-1000.*fl.gape./2,'r-')
hold on
plot(1000.*r.pos(:,1),1000.*r.pos(:,2),'b.')
hold off
xlabel('x - coord (mm)')
ylabel('y - coord (mm)')
axis equal


end


end

