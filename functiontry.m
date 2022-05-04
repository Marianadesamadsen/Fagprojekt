%% Script for trying functions
clear 
clc
close all

%% Converting factors 

h2min = 60;      % Convert from h   to min 
min2h = 1/h2min; % Convert from min to h 
U2mU  = 1e3;     % Convert from U   to mU 
mU2U  = 1/U2mU;  % Convert from mU  to U 


%% Inisializing parameters

p = [49,47,20.1,0.0106,0.0081,0.0022,1.33,253,47,5]; % Parameters at steady state
Gs = 108; % [mg/dL]: Steady state blood glucose concentration
xs=zeros(1,7);
d=0;
u=[0,0];
t0 =  0;       
tf = 18*h2min; 
Ts=5;
N = (tf - t0)/Ts;
tspan=Ts*(0:N);

%% Inisializing models

simModel = @MVPmodel; % The differential equations 
outputModel = @CGMsensor; % [mg/dL] Subcutaneous glucose concentration (Under the skin)
simMethod = @ExplicitEuler; % Solving the differential equations

%% Trying to solve differential equation based on explicit euler:

[T,X]=simMethod(simModel, tspan, xs, u, d, p)


%% Trying to solve differential equation based on ode15

[T,X] = ode15s(@MVPmodel,tspan,xs,[],u,d,pf);







% Computing OpenLoopsSimulation
x0 = zeros(7,1);

[T_matrix,X_vector] = OpenLoopSimulation(x0, tspan, U, D, p, simModel, simMethod, opts);




%% Computing the steady state


x=ExplicitEuler_vores1(@MVPmodel_vores1,tspan, xk, D, U, pf);

% Computing OpenLoopsSimulation
x0 = zeros(7,1);

[T_matrix,X_vector] = OpenLoopSimulation(x0, tspan, U, D, p, simModel, simMethod, opts);

%% 2

% Solving with ode15s
[T,X] = ode15s(@MVPmodel_vores1,tspan,x0,[],U,D,pf);

% steady at (66,6)

% Solving with Euler
[TE,XE] = ExplicitEuler_vores1(@MVPmodel_vores1,tspan,x0,U,D,pf);

% steady at (6,735)

%% 24 af 5 min interval 

tspan2 = [1 288];

[Z_vector2,Y_vector2,X_vector2] = OpenLoopSimulation(@MVPmodel_vores1, pf, @CGMsensor_vores, pg, @CHMsensor_vores, ph,x0, tspan2, U, D);

%s1=zeros(8,289);

%for i=1:289
 %   s1(:,1)=MVPmodel_vores1(tspan,X_vector2(i,:),U,D,pf);
%end

%steady_con=MVPmodel_vores1(tspan,X_vector2(6,249),U,D,pf);
%steady_mes=MVPmodel_vores1(tspan,X_vector2(7,246),U,D,pf);

glucose_concentration=X_vector2(6,:);
glucose_measured=X_vector2(7,:);

figure (1)
hold on
plot(glucose_concentration)
plot(glucose_measured)








