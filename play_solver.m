function play_solver
% Comparing numerical results with analytical solutions 

% Solver options
options = odeset('RelTol',1e-1);

% %% Simple integration with constant acceleration
% 
% % Acceleration and duration values
% accel   = -10;
% dur     = .05;
% 
% % Function that runs the numerical solver
% [t,X] = solver1(accel,dur,options);
% 
% % Analytical prediction
% predX = 0.5.*accel.*t.^2;
% 
% figure
% subplot(3,1,1)
% plot(t,X(:,1),'k',t,predX,'r--')
% title('Constant acceleration')
% 
% clear accel dur predX t X
% 
% 
% %% Drag with a velocity ramp
% 
% % Acceleration of flow, duration, mass
% aFlow = -10;
% dur   = 0.05; 
% m     = 6.2832e-5;
% 
% % 0.5*rho*S*Cd
% k     = 0.07539;
% 
% % Function that runs the numerical solver
% [t,X] = solver2(aFlow,dur,k,m,options);
% 
% 
% % Analytical prediction
% predX = -k*aFlow.* t.^3 .*abs(aFlow.*t) ./ (12.*m);
% 
% subplot(3,1,2)
% plot(t,X(:,1),'k',t,predX,'r--')
% title('Drag with velocity ramp')
% 
% clear accel dur predX t X
% 
% 
% %% ODE: low Re Drag equation
% % Doesn't solve well without scaling the parameter values
% 
% % Define scaling constants
% sL = 20e-3;
% sT = 10^-3;
% sM = 6.2832e-5;
% 
% % Acceleration of flow, duration, mass
% aFlow = -10 /sL *sT^2;
% dur   = 0.05 /sT; 
% m     = 6.2832e-5 /sM;
% 
% % 0.5*rho*S*Cd
% k     = 0.07539 /sM * sL;
% 
% % Function that runs the numerical solver
% [t,X] = solver3(aFlow,dur,k,m,options);
% 
% t = t .*sT;
% X(:,1) = X(:,1) .*sL;
% X(:,2) = X(:,2) .*sL ./sT;
% 
% % Prediction fof t = 0.05 (determined in mathematica:siff_calculations.nb)
% predX = -.0125+t.*0;
% 
% subplot(3,1,3)
% plot(t,X(:,1),'k',t,predX,'r--')
% title('ODE (low Re)')
% 
% clear accel dur predX t X

% %% ODE: high Re Drag equation
% % Doesn't solve well without scaling the parameter values
% 
% % Initial position
% X0 = 9.2e-3;
% 
% % Define scaling constants
% sL = 20e-3;
% sT = 10^-3;
% sM = (6.2832e-5)*10^6;
% 
% % Acceleration of flow, duration, mass
% r     = 1e-3 /sL;
% L     = 20e-3 /sL;
% aFlow = -33.33 /sL *sT^2;
% dur   = 0.05 /sT; 
% m     = 6.2832e-5 /sM;
% rho   = 1000 /sM * sL^3;
% 
% % Already scaled parameters
% Cd    = 1.2;
% S     = 2*pi*r*L;
% k     = 0.5*rho*S*Cd;
% 
% % Function that runs the numerical solver
% [t,X] = solver4(aFlow,dur,k,m,options);
% 
% t = t .*sT;
% X(:,1) = X0 + X(:,1) .*sL;
% X(:,2) = X(:,2) .*sL ./sT;
% 
% % Prediction fof t = 0.05 (determined in mathematica:siff_calculations.nb)
% 
% figure
% subplot(3,1,1)
% plot(t,X(:,1),'k')
% title('ODE (High Re)')
% 
% clear accel dur predX t X


%% ODE: Pressure force
% Doesn't solve well without scaling the parameter values

% Initial position
X0 = 9.2e-3;

% Define scaling constants
sL = 20e-3;
sT = 10^-3;
sM = (6.2832e-5)*10^6;

% Acceleration of flow, duration, mass
r     = 1e-3 /sL;
L     = 20e-3 /sL;
aFlow = -33.33 /sL *sT^2;
dur   = 0.05 /sT; 
m     = 6.2832e-5 /sM;
rho   = 1000 /sM * sL^3;

% Already scaled parameters
Cd    = 1.2;
S     = 2*pi*r*L;
k     = 0.5*rho*S*Cd;

% Function that runs the numerical solver
[t,X] = solver5(aFlow,dur,k,m,options);

t = t .*sT;
X(:,1) = X0 + X(:,1) .*sL;
X(:,2) = X(:,2) .*sL ./sT;

% Prediction fof t = 0.05 (determined in mathematica:siff_calculations.nb)

figure
subplot(3,1,1)
plot(t,X(:,1),'k')
title('ODE (High Re)')

clear accel dur predX t X


end


function [t,X] = solver1(accel,dur,options)

% Solve governing equation
[t,X] = ode45(@gov_eqn,[0 dur],[0; 0],options);

    function dX = gov_eqn(t,X)
        % ODE of the dynamics of the system
        
        % Prey COM position 
        posCOM = X(1);
        
        % Prey COM velocity
        velCOM = X(2);
        
        % Body acceleration
        accelCOM = accel;
   
        % Output: x speed
        dX(1,1) = velCOM;
        % Output: x accel
        dX(2,1) = accelCOM;     
    end
end

function [t,X] = solver2(aFlow,dur,k,m,options)

% Solve governing equation
[t,X] = ode45(@gov_eqn,[0 dur],[0; 0],options);

    function dX = gov_eqn(t,X)
        % ODE of the dynamics of the system
        
        % Prey COM position 
        posCOM = X(1);
        
        % Prey COM velocity
        velCOM = X(2);
        
        % Flow velocity
        velFlow = aFlow * t;
        
        % Drag
        D = -k*velFlow*abs(velFlow);
        
        % Body acceleration
        accelCOM = D/m;
   
        % Output: x speed
        dX(1,1) = velCOM;
        
        % Output: x accel
        dX(2,1) = accelCOM;     
    end
end

function [t,X] = solver3(aFlow,dur,k,m,options)

% Solve governing equation
[t,X] = ode45(@gov_eqn,[0 dur],[0; 0],options);

    function dX = gov_eqn(t,X)
        % ODE of the dynamics of the system
        
        % Prey COM position 
        posCOM = X(1);
        
        % Prey COM velocity
        velCOM = X(2);
        
        % Flow velocity
        velFlow = aFlow * t;
        
        % Relative velocity
        relFlow = velFlow - velCOM;
        
        % Drag
        %D = k*relFlow*abs(relFlow);
        D = k*relFlow;
        
        % Body acceleration
        accelCOM = D/m;
   
        % Output: x speed
        dX(1,1) = velCOM;
        
        % Output: x accel
        dX(2,1) = accelCOM;     
    end
end

function [t,X] = solver4(aFlow,dur,k,m,options)

% Solve governing equation
[t,X] = ode45(@gov_eqn,[0 dur],[0; 0],options);

    function dX = gov_eqn(t,X)
        % ODE of the dynamics of the system
        
        % Prey COM position 
        posCOM = X(1);
        
        % Prey COM velocity
        velCOM = X(2);
        
        % Flow velocity
        velFlow = aFlow * t;
        
        % Relative velocity
        relFlow = velFlow - velCOM;
        
        % Drag
        %D = k*relFlow*abs(relFlow);
        D = k*relFlow*abs(relFlow);
        
        % Body acceleration
        accelCOM = D/m;
   
        % Output: x speed
        dX(1,1) = velCOM;
        
        % Output: x accel
        dX(2,1) = accelCOM;     
    end
end

function [t,X] = solver5(aFlow,dur,k,m,options)

% Solve governing equation
[t,X] = ode45(@gov_eqn,[0 dur],[0; 0],options);

    function dX = gov_eqn(t,X)
        % ODE of the dynamics of the system
        
        % Prey COM position 
        posCOM = X(1);
        
        % Prey COM velocity
        velCOM = X(2);
        
        % Flow velocity
        velFlow = aFlow * t;
        
        % Relative velocity
        relFlow = velFlow - velCOM;
        
        % Drag
        %D = k*relFlow*abs(relFlow);
        D = k*relFlow*abs(relFlow);
        
        % Body acceleration
        accelCOM = D/m;
   
        % Output: x speed
        dX(1,1) = velCOM;
        
        % Output: x accel
        dX(2,1) = accelCOM;     
    end
end