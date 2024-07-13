 % This is the main PN SIRS Model code for

% This code will not work if you do not 
% first run code for simple_pn_SIRv1_pdf.m and COMMON_PRE.m 

%% Descrtiption and Directions
% First, choose section to run by changing Run_Section to desired number.
% Then go to section and change variables listed in description to desired
% level. Hitting command+G on mac will skip between each section.
% 1 = run basic SIRS PN model with 1 set of parameters and ode set to same parameters
    % beta, gamma, delta, tau (called timidivi, time steps per unit time),
    % Initial Conditions, a total time for simulations are what are params
    % to be set
    % go to line 53 to change parameters and look over section

% 2 = run SIRS PN model over range of parameters and compares RRMSE
    % beta, gamma can be any size grid between [0,1] (using bglength and spacer,
    % delta needs to be grid of size 5, tau (called timidivi, time steps per unit time),
    % Initial Conditions, a total time for simulations are params to be set
    % Go to line __ to change parameters and look over section

% 3 = finds RRMSE for various spread of parameter values
    % Set tau_spread, where tau is spread of parameters between [0,1]
    % different from how tau is defined in the paper
    % Go to line __ to change parameters and look over section

% 4 = comparing different types of rounding methods for single param value
    % all params to be set are same as section 1
    % Go to line __ to change parameters and look over section

% 5 = run SIRS PN model over range of parameters and Rounding methods and compares RRMSE
    % all params set same as Run_Section 2
    % Go to line __ to change parameters and look over section

% 6 = Load parameters from supercomputer runs to compare mean RRMSE of
    % different time steps
    % Go to line __ to change parameters and look over section

% 7 = Computation time plot
    % Go to line __ to change parameters and look over section

clear all; clc; 
global global_info m0Susceptible m0Infected m0Recovered %beta
tic %for timing program run time

Run_Section=1;

%% Run_Section 1, run 1 SIRS PN model and compare to ODE
if Run_Section==1

timidivi=20;                    %tau in paper
delta=1; beta=1; gamma=1;       %parameter values 

try
delta_T=delta/timidivi;         %parameters scaled for tau value
beta_T=beta/timidivi; 
gamma_T=gamma/timidivi;

MAX_ITERATIONS = 100*timidivi;  %Max time of simulation  

S_0=1000;                       %initial token value of place S (initial population of susceptible)
I_0=10;                         %initial token value of place I (initial population of infected)
R_0=10;                         %initial token value of place R (initial population of recovered)
m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

%loop for PN model (necessary to find variable arc weights after firings)
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;
while le(Number_of_iterations, MAX_ITERATIONS)
    %disp(['Iteration: ', int2str(Number_of_iterations)]);
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 
    
    %Find arc weights for each iteration using standard+residual rounding
    %method
    if beta_T==0 || m0Infected==0
        global_info.pSus_tInf = zero;
        pSus_tInf_resEr=(beta_T * m0Susceptible *  m0Infected + pSus_tInf_resEr)-round(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
        global_info.tInf_pInf = zero;
        tInf_pInf_resEr=(beta_T * m0Susceptible *  m0Infected+tInf_pInf_resEr)-round(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
    elseif m0Infected<=(1/(beta_T))
        global_info.pSus_tInf = round(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr);
        pSus_tInf_resEr=(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr)-round(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
        global_info.tInf_pInf = round(beta_T * m0Susceptible *   m0Infected + tInf_pInf_resEr);
        tInf_pInf_resEr=(beta_T * m0Susceptible *   m0Infected +tInf_pInf_resEr)-round(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
    else 
        beta2=1/(m0Infected+1);
        global_info.pSus_tInf = round(beta2 * m0Susceptible * m0Infected + pSus_tInf_resEr);
        pSus_tInf_resEr=(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr)-round(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr);
        global_info.tInf_pInf = round(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);
        tInf_pInf_resEr=(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr)-round(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);
    
    end

    global_info.pInf_tRec = round(gamma_T * (m0Infected) + pInf_tRec_resEr);
    pInf_tRec_resEr=(gamma_T * (m0Infected) + pInf_tRec_resEr)-round(gamma_T * (m0Infected) + pInf_tRec_resEr);
    global_info.tRec_pRec = round(gamma_T * (m0Infected) + tRec_pRec_resEr);
    tRec_pRec_resEr=(gamma_T * (m0Infected) + tRec_pRec_resEr)-round(gamma_T * (m0Infected) + tRec_pRec_resEr);
    global_info.pRec_tSus = round(delta_T * (m0Recovered) + pRec_tSus_resEr);
    pRec_tSus_resEr=(delta_T * (m0Recovered) + pRec_tSus_resEr)-round(delta_T * (m0Recovered) + pRec_tSus_resEr);
    global_info.tSus_pSus = round(delta_T * (m0Recovered) + tSus_pSus_resEr);
    tSus_pSus_resEr=(delta_T * (m0Recovered)+ tSus_pSus_resEr)-round(delta_T * (m0Recovered) + tSus_pSus_resEr);
    
 pns = pnstruct('simple_pn_SIRv1_pdf');
    dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
catch

    disp(['Error at delta=', int2str(delta), ' , gamma=',int2str(gamma), ' , beta=',int2str(beta), ' , timesteps=',int2str(timidivi)]);
    
    if delta==0
    delta=delta; 
    else
    delta=delta-10^-7; 
    end
   
    if beta==0
    beta=beta; 
    else
    beta=beta-10^-7;
    end
    
    if gamma==0
    gamma=gamma;
    else
    gamma=gamma-10^-7;
    end
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

MAX_ITERATIONS = 100*timidivi;

% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

% Initial Population
S_0=1000;
I_0=10;
R_0=10;

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

% the loop   
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;
while le(Number_of_iterations, MAX_ITERATIONS)
    %disp(['Iteration: ', int2str(Number_of_iterations)]);
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 

    %Find arc weights for each iteration using standard+residual rounding
    %method
    if beta_T==0 || m0Infected==0
        global_info.pSus_tInf = zero;
        pSus_tInf_resEr=(beta_T * m0Susceptible *  m0Infected + pSus_tInf_resEr)-round(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
        global_info.tInf_pInf = zero;
        tInf_pInf_resEr=(beta_T * m0Susceptible *  m0Infected+tInf_pInf_resEr)-round(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
    elseif m0Infected<=(1/(beta_T))
        global_info.pSus_tInf = round(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr);
        pSus_tInf_resEr=(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr)-round(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
        global_info.tInf_pInf = round(beta_T * m0Susceptible *   m0Infected + tInf_pInf_resEr);
        tInf_pInf_resEr=(beta_T * m0Susceptible *   m0Infected +tInf_pInf_resEr)-round(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
    else 

        %beta2=(1/(m0Infected+??????));
        beta2=1/(m0Infected+1);
        global_info.pSus_tInf = round(beta2 * m0Susceptible * m0Infected + pSus_tInf_resEr);
        pSus_tInf_resEr=(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr)-round(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr);
        global_info.tInf_pInf = round(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);
        tInf_pInf_resEr=(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr)-round(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);
    
    end

    global_info.pInf_tRec = round(gamma_T * (m0Infected) + pInf_tRec_resEr);
    pInf_tRec_resEr=(gamma_T * (m0Infected) + pInf_tRec_resEr)-round(gamma_T * (m0Infected) + pInf_tRec_resEr);
    global_info.tRec_pRec = round(gamma_T * (m0Infected) + tRec_pRec_resEr);
    tRec_pRec_resEr=(gamma_T * (m0Infected) + tRec_pRec_resEr)-round(gamma_T * (m0Infected) + tRec_pRec_resEr);
    global_info.pRec_tSus = round(delta_T * (m0Recovered) + pRec_tSus_resEr);
    pRec_tSus_resEr=(delta_T * (m0Recovered) + pRec_tSus_resEr)-round(delta_T * (m0Recovered) + pRec_tSus_resEr);
    global_info.tSus_pSus = round(delta_T * (m0Recovered) + tSus_pSus_resEr);
    tSus_pSus_resEr=(delta_T * (m0Recovered)+ tSus_pSus_resEr)-round(delta_T * (m0Recovered) + tSus_pSus_resEr);
    

    pns = pnstruct('simple_pn_SIRv1_pdf');
    %dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    %figure()
    %plotp(Sim_Results, {'pSusceptible','pInfected','pRecovered'});

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
end

y0=[S_0;I_0;R_0];           %ICs for ODE
delta_ODE=delta; 
beta_ODE=beta; 
gamma_ODE=gamma;
params=[delta_ODE;beta_ODE;gamma_ODE];
[t,y]=ode15s(@SIR,[1:1:(MAX_ITERATIONS/timidivi)],y0,[],params);    %ODE model run

% Plot Results
figure()
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MSusceptible, 'b--o','LineWidth',2);
hold on
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MInfected,'c--o','LineWidth',2);
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MRecovered, 'r--o','LineWidth',2);

plot(t,y(:,1), 'y','LineWidth',2)
plot(t,y(:,2),'m','LineWidth',2)
plot(t,y(:,3), 'k','LineWidth',2);
legend('Susceptible_{PN}','Infected_{PN}','Recovered_{PN}',...
    'Susceptible_{ODE}','Infected_{ODE}','Recovered_{ODE}');
title("SIR Petri Net vs ODE", 'FontSize', 24)
subtitle("beta="+beta*timidivi+', gamma='+gamma*timidivi+' delta='+delta*timidivi+", Susc. MSE="+err_S...
   +", Infe. MSE="+err_I+", Reco. MSE="+err_R, 'FontSize', 14)
hold off

end %end of section 1



















%% Run Section 2, runs PN and ODE model over range of params
if Run_Section==2

timidivi_col_S=[];
timidivi_col_I=[];
timidivi_col_R=[];

timidivi=40;                        %tau in paper

bglength=10;                        %size of grid for beta and gamma

spacer=linspace(0,1,bglength);      %grid of params for beta and gamma

MSE_s=[zeros(bglength,bglength,5)]; %storage array for place s (susceptible) RRMSE
MSE_i=[zeros(bglength,bglength,5)]; %storage array for place s (susceptible) RRMSE
MSE_r=[zeros(bglength,bglength,5)]; %storage array for place s (susceptible) RRMSE

 deltactr=1;
  for delta=[0,logspace(-3,0,4)]    %grid of params for delta (must be of size 5 for figs)
     betactr=1;
     for beta=spacer   
         gammactr=1;
         for gamma=spacer

try
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;   %scale params based on tau

MAX_ITERATIONS = 100*timidivi;

% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

% Initial Population
S_0=1000;
I_0=10;
R_0=10;

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

% the loop   
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;

while le(Number_of_iterations, MAX_ITERATIONS)
   
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 
   
    %Find arc weights for each iteration using standard+residual rounding
    %method
    if beta_T==0 || m0Infected==0
        global_info.pSus_tInf = zero;
        pSus_tInf_resEr=(beta_T * m0Susceptible *  m0Infected + pSus_tInf_resEr)-round(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
        global_info.tInf_pInf = zero;
        tInf_pInf_resEr=(beta_T * m0Susceptible *  m0Infected+tInf_pInf_resEr)-round(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
    elseif m0Infected<=(1/(beta_T))
        global_info.pSus_tInf = round(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr);
        pSus_tInf_resEr=(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr)-round(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
        global_info.tInf_pInf = round(beta_T * m0Susceptible *   m0Infected + tInf_pInf_resEr);
        tInf_pInf_resEr=(beta_T * m0Susceptible *   m0Infected +tInf_pInf_resEr)-round(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
    else 

       
        beta2=1/(m0Infected+1);
        global_info.pSus_tInf = round(beta2 * m0Susceptible * m0Infected + pSus_tInf_resEr);
        pSus_tInf_resEr=(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr)-round(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr);
        global_info.tInf_pInf = round(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);
        tInf_pInf_resEr=(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr)-round(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);
    
    end

    global_info.pInf_tRec = round(gamma_T * (m0Infected) + pInf_tRec_resEr);
    pInf_tRec_resEr=(gamma_T * (m0Infected) + pInf_tRec_resEr)-round(gamma_T * (m0Infected) + pInf_tRec_resEr);
    global_info.tRec_pRec = round(gamma_T * (m0Infected) + tRec_pRec_resEr);
    tRec_pRec_resEr=(gamma_T * (m0Infected) + tRec_pRec_resEr)-round(gamma_T * (m0Infected) + tRec_pRec_resEr);
    global_info.pRec_tSus = round(delta_T * (m0Recovered) + pRec_tSus_resEr);
    pRec_tSus_resEr=(delta_T * (m0Recovered) + pRec_tSus_resEr)-round(delta_T * (m0Recovered) + pRec_tSus_resEr);
    global_info.tSus_pSus = round(delta_T * (m0Recovered) + tSus_pSus_resEr);
    tSus_pSus_resEr=(delta_T * (m0Recovered)+ tSus_pSus_resEr)-round(delta_T * (m0Recovered) + tSus_pSus_resEr);

    pns = pnstruct('simple_pn_SIRv1_pdf');
    dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
catch

    disp(['Error at delta=', int2str(delta), ' , gamma=',int2str(gamma), ' , beta=',int2str(beta), ' , timesteps=',int2str(timidivi)]);
    
    if delta==0
    delta=delta; 
    else
    delta=delta-10^-7; 
    end
   
    if beta==0
    beta=beta; 
    else
    beta=beta-10^-7;
    end
    
    if gamma==0
    gamma=gamma;
    else
    gamma=gamma-10^-7;
    end
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

MAX_ITERATIONS = 100*timidivi;

%Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

% Initial Population
S_0=1000;
I_0=10;
R_0=10;

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

% the loop   
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;
while le(Number_of_iterations, MAX_ITERATIONS)
    
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 

    if beta_T==0 || m0Infected==0
        global_info.pSus_tInf = zero;
        pSus_tInf_resEr=(beta_T * m0Susceptible *  m0Infected + pSus_tInf_resEr)-round(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
        global_info.tInf_pInf = zero;
        tInf_pInf_resEr=(beta_T * m0Susceptible *  m0Infected+tInf_pInf_resEr)-round(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
    elseif m0Infected<=(1/(beta_T))
        global_info.pSus_tInf = round(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr);
        pSus_tInf_resEr=(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr)-round(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
        global_info.tInf_pInf = round(beta_T * m0Susceptible *   m0Infected + tInf_pInf_resEr);
        tInf_pInf_resEr=(beta_T * m0Susceptible *   m0Infected +tInf_pInf_resEr)-round(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
    else 

        beta2=1/(m0Infected+1);
        global_info.pSus_tInf = round(beta2 * m0Susceptible * m0Infected + pSus_tInf_resEr);
        pSus_tInf_resEr=(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr)-round(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr);
        global_info.tInf_pInf = round(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);
        tInf_pInf_resEr=(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr)-round(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);
    
    end

    global_info.pInf_tRec = round(gamma_T * (m0Infected) + pInf_tRec_resEr);
    pInf_tRec_resEr=(gamma_T * (m0Infected) + pInf_tRec_resEr)-round(gamma_T * (m0Infected) + pInf_tRec_resEr);
    global_info.tRec_pRec = round(gamma_T * (m0Infected) + tRec_pRec_resEr);
    tRec_pRec_resEr=(gamma_T * (m0Infected) + tRec_pRec_resEr)-round(gamma_T * (m0Infected) + tRec_pRec_resEr);
    global_info.pRec_tSus = round(delta_T * (m0Recovered) + pRec_tSus_resEr);
    pRec_tSus_resEr=(delta_T * (m0Recovered) + pRec_tSus_resEr)-round(delta_T * (m0Recovered) + pRec_tSus_resEr);
    global_info.tSus_pSus = round(delta_T * (m0Recovered) + tSus_pSus_resEr);
    tSus_pSus_resEr=(delta_T * (m0Recovered)+ tSus_pSus_resEr)-round(delta_T * (m0Recovered) + tSus_pSus_resEr);
    

    pns = pnstruct('simple_pn_SIRv1_pdf');
    %dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
end

%ODE function run
y0=[S_0;I_0;R_0];
delta_ODE=delta; 
beta_ODE=beta; 
gamma_ODE=gamma;
params=[delta_ODE;beta_ODE;gamma_ODE];
[t,y]=ode15s(@SIR,[1:1:(MAX_ITERATIONS/timidivi)],y0,[],params);


% RRMSE Calculation
 err_S = (rmse(MSusceptible(1:timidivi:MAX_ITERATIONS) , y(:,1)')/sqrt(sumsqr(y(:,1)')))*100;
 err_I = (rmse(MInfected(1:timidivi:MAX_ITERATIONS) , y(:,2)')/sqrt(sumsqr( y(:,2)')))*100;
 err_R = (rmse(MRecovered(1:timidivi:MAX_ITERATIONS) , y(:,3)')/sqrt(sumsqr(y(:,3)')))*100;
 % % 
 % % timidivi_col_S=[timidivi_col_S err_S];
 % % timidivi_col_I=[timidivi_col_I err_I];
 % % timidivi_col_R=[timidivi_col_R err_R];
% %end

disp(['still going: ', int2str(deltactr)]);

        % MSE_s(gammactr,betactr)=err_S;
        % MSE_i(gammactr,betactr)=err_I;
        % MSE_iBeta(gammactr,betactr)=beta;
        % MSE_iGamma(gammactr,betactr)=gamma;
        % MSE_r(gammactr,betactr)=err_R;

        %RRMSE storage
        MSE_s(gammactr,betactr,deltactr)=err_S;
        MSE_i(gammactr,betactr,deltactr)=err_I;
        MSE_r(gammactr,betactr,deltactr)=err_R;

           gammactr=gammactr+1;
         end
        betactr=betactr+1;
     end
     deltactr=deltactr+1;
  end

%Plot the results of grid param value RRMSE vs ODE
% casxismin=min([min(MSE_s(:)),min(MSE_i(:)),min(MSE_r(:))]);
% trial=1;
% if trial==1
% casxismax=max([max(MSE_s(:)),max(MSE_i(:)),max(MSE_r(:))]);
% trial=trial+1;
% else
% end

casxismin=0;
%Max RRMSE was found across all runs and then set
casxismax=44; 

%custom color map
Cus_map=[(ones(2,3).*[0 0.5 0]);(ones(8,3).*[0 1 0]);(ones(40,3).*[1 1 0]);(ones(50,3).*[1 0.75 0]);...
    (ones(100,3).*[1 0.4 0]);(ones(240,3).*[1 0 0])];


figure()
%title("RRMSE for parameter Ranges with 1 pn step per 1 ode unit time", 'FontSize', 14)
tiledlayout(5,3);
nexttile
s=surf(spacer,spacer,MSE_s(:,:,1));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Susceptible \delta=0','fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h, 'RRMSE (%)');
%h.Label.String = "RRMSE (%)";h.Label.Rotation = 0;
view(2)

nexttile
s=surf(spacer,spacer,MSE_i(:,:,1));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Infected \delta=0','fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h, 'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r(:,:,1));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Recovered \delta=0','fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h, 'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_s(:,:,2));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Susceptible \delta=0.001','fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i(:,:,2));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Infected \delta=0.001','fontsize',12')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r(:,:,2));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Recovered \delta=0.001','fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_s(:,:,3));
ylabel('\gamma (rate of recovery per unit time)','fontsize',24); 
%xlabel('beta'); zlabel('RRMSE (%)');
title('Susceptible \delta=0.01','fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i(:,:,3));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Infected \delta=0.01','fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r(:,:,3));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Recovered \delta=0.01','fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_s(:,:,4));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Susceptible \delta=0.1','fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i(:,:,4));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Infected \delta=0.1','fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r(:,:,4));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Recovered \delta=0.1','fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_s(:,:,5));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Susceptible \delta=1','fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i(:,:,5));
%ylabel('gamma'); zlabel('RRMSE');
xlabel('\beta (rate of infection per unit of time)','FontSize',24); 
title('Infected \delta=1','fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r(:,:,5));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Recovered \delta=1','fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
view(2)
h=colorbar('Position',[0.931120331950207,0.11,0.014914910340541,0.816052332195677]);
title(h, 'RRMSE (%)');
%h=axes('Position', [0.85 0.15 0.05 0.7]);% axes('Position', [0.85 0.15 0.05 0.7]),   %'Location', 'westoutside',
%h.Label.String = "RRMSE (%)";h.Label.Rotation = 0;

MSE_s
MSE_i
MSE_r
%save all variables
save(strcat(num2str(timidivi),'paramsb.mat'))

end %end of section 2



















%% Run Section 3
%not currently used in paper so code not updated
if Run_Section==3
timidivi_col_S=[];
timidivi_col_I=[];
timidivi_col_R=[];
%timidivi_timestep=1:6:79;
%for timidivi=1:6:79
ode_run_counter=1;
%timidivi=1;
%for timidivi=[2]
timidivi=20;
% The following section is used for finding the RRMSE for various
%parameter values in comparison to the respective ODE, this is limited two
%comparison of 3 parameters in this configuration, for visual display.
  bglength=1;
% tau_S=[];
% tau_I=[];
% tau_R=[];
% tau_spread=linspace(0,1,10);

% for tau=linspace(0,1,10)
%spacer=[0 logspace(-5,0,bglength-1)];
%spacer=[0 logspace(-5,-1,bglength) linspace(0.2,0.8,bglength) (ones(1,bglength)-logspace(-1,-5,bglength)) 1];
% spacer=linspace(0,tau,bglength);
spacer=linspace(0,1,bglength);
% gridsize=size(spacer,2);
% 
% MSE_s=[zeros(gridsize,gridsize)];
% MSE_i=[zeros(gridsize,gridsize)];
% MSE_r=[zeros(gridsize,gridsize)];
MSE_s=[zeros(bglength,bglength,5)];
MSE_i=[zeros(bglength,bglength,5)];
MSE_r=[zeros(bglength,bglength,5)];

% MSE_S_c=zeros(1,5);
% MSE_S_r=zeros(1,5);
% MSE_S_v=zeros(1,5);
% MSE_S_C=zeros(bglength*bglength,5);
% %0.1
%deltar=[0,logspace(-3,0,4)];
 % % % % % deltactr=1;
 % % % % %  for delta=[0,logspace(-3,0,4)]
 % % % % %     betactr=1;
 % % % % %     %for beta=linspace(beta_ODE*0.9,beta_ODE*1.1,bglength) %beta_ODE = .0008;     gamma_ODE = 0.08;
 % % % % %     for beta=spacer   
 % % % % %         %beta=spacer(betaindex);
 % % % % %         gammactr=1;
 % % % % %         %for gamma=linspace(gamma_ODE*0.9,gamma_ODE*1.1,bglength)
 % % % % %         for gamma=spacer
             %gamma=spacer(gammaindex);
delta=1; beta=1; gamma=1;
% delta1=1/timidivi; beta1=1/timidivi; gamma1=1/timidivi;
for Rounder_test=1:1:4
    if Rounder_test==1
try
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;
%delta=1; beta=1; gamma=1;

%% Max iterations for Petri Net and ODE, petri net does 1 step per iteration
%in base settings
MAX_ITERATIONS = 100*timidivi;

%%%%% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 


%% Parameter values  %%%%%%%%%%%%%%%%%%%%%
%Used for setting specific parameter values (vs loop above for grid of
%param values)

%beta = 0.001;     gamma = 0.05;    delta = 0.001;
%beta = 0.1;     gamma = 0.2;    delta = 0.1;
%beta = 0.2;     gamma = 0.5;    delta = 0.4;
%beta = 3.1000e-04;     gamma =  0.0809;    
%delta = 0.00;


%% Initial Population
S_0=1000;
I_0=10;
R_0=10;

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

%% the loop   %%%%%%%%%%%%%%%%%%%%%%%%%
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;

while le(Number_of_iterations, MAX_ITERATIONS)
    %disp(['Iteration: ', int2str(Number_of_iterations)]);
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 
    %  Arc weights
    %This below is for splitting time
    % if Number_of_iterations==1
    %     beta=beta1;
    %     delta=delta1;
    %     gamma=gamma1;
    % else
    %     % beta=1;
    %     % delta=1;
    %     % gamma=1;
    % end
    if beta_T==0 || m0Infected==0
        global_info.pSus_tInf = zero;
        pSus_tInf_resEr=(beta_T * m0Susceptible *  m0Infected + pSus_tInf_resEr)-round(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
        global_info.tInf_pInf = zero;
        tInf_pInf_resEr=(beta_T * m0Susceptible *  m0Infected+tInf_pInf_resEr)-round(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
    elseif m0Infected<=(1/(beta_T))
        global_info.pSus_tInf = round(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr);
        pSus_tInf_resEr=(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr)-round(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
        global_info.tInf_pInf = round(beta_T * m0Susceptible *   m0Infected + tInf_pInf_resEr);
        tInf_pInf_resEr=(beta_T * m0Susceptible *   m0Infected +tInf_pInf_resEr)-round(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
    else 

        %beta2=(1/(m0Infected+??????));
        beta2=1/(m0Infected+1);
        global_info.pSus_tInf = round(beta2 * m0Susceptible * m0Infected + pSus_tInf_resEr);
        pSus_tInf_resEr=(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr)-round(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr);
        global_info.tInf_pInf = round(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);
        tInf_pInf_resEr=(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr)-round(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);
    
    end

    global_info.pInf_tRec = round(gamma_T * (m0Infected) + pInf_tRec_resEr);
    pInf_tRec_resEr=(gamma_T * (m0Infected) + pInf_tRec_resEr)-round(gamma_T * (m0Infected) + pInf_tRec_resEr);
    global_info.tRec_pRec = round(gamma_T * (m0Infected) + tRec_pRec_resEr);
    tRec_pRec_resEr=(gamma_T * (m0Infected) + tRec_pRec_resEr)-round(gamma_T * (m0Infected) + tRec_pRec_resEr);
    global_info.pRec_tSus = round(delta_T * (m0Recovered) + pRec_tSus_resEr);
    pRec_tSus_resEr=(delta_T * (m0Recovered) + pRec_tSus_resEr)-round(delta_T * (m0Recovered) + pRec_tSus_resEr);
    global_info.tSus_pSus = round(delta_T * (m0Recovered) + tSus_pSus_resEr);
    tSus_pSus_resEr=(delta_T * (m0Recovered)+ tSus_pSus_resEr)-round(delta_T * (m0Recovered) + tSus_pSus_resEr);
    %{ 
%One way to handle rounding error
    if m0Infected<=(1/beta)

    if round(beta * m0Susceptible * m0Infected)==0 && (beta * m0Susceptible * m0Infected)~=0 && MinfiringS==0 
        timestepstofireS=round((1/(beta * m0Susceptible * m0Infected)));
        tstf_counterS=tstf_counterS+1;
        MinfiringS=1;
        global_info.pSus_tInf = round(beta * m0Susceptible * m0Infected);
        global_info.tInf_pInf = round(beta * m0Susceptible * m0Infected);
        (beta * m0Susceptible * m0Infected)
    elseif MinfiringS==0
        global_info.pSus_tInf = round(beta * m0Susceptible * m0Infected);
        global_info.tInf_pInf = round(beta * m0Susceptible * m0Infected);
        (beta * m0Susceptible * m0Infected)
    elseif MinfiringS==1 && tstf_counterS<timestepstofireS
        tstf_counterS=tstf_counterS+1;
        global_info.pSus_tInf = 0;
        global_info.tInf_pInf = 0;
        (beta * m0Susceptible * m0Infected)
    elseif MinfiringS==1 && tstf_counterS==timestepstofireS
        global_info.pSus_tInf = 1;
        global_info.tInf_pInf = 1;
        (beta * m0Susceptible * m0Infected)
        MinfiringS=0;
        tstf_counterS=0;
    else
        error("logic loop broken")
    end
 
    
    else

        beta2=1/m0Infected;
    if round(beta2 * m0Susceptible * m0Infected)==0 && (beta2 * m0Susceptible * m0Infected)~=0 && MinfiringS==0 
        timestepstofireS=round(1/(beta2 * m0Susceptible * m0Infected));
        tstf_counterS=tstf_counterS+1;
        MinfiringS=1;
        global_info.pSus_tInf = round(beta2 * m0Susceptible * m0Infected);
        global_info.tInf_pInf = round(beta2 * m0Susceptible * m0Infected);
    elseif MinfiringS==0
        global_info.pSus_tInf = round(beta2 * m0Susceptible * m0Infected);
        global_info.tInf_pInf = round(beta2 * m0Susceptible * m0Infected);
    elseif MinfiringS==1 && tstf_counterS<timestepstofireS
        tstf_counterS=tstf_counterS+1;
        global_info.pSus_tInf = 0;
        global_info.tInf_pInf = 0;
    elseif MinfiringS==1 && tstf_counterS==timestepstofireS
        global_info.pSus_tInf = 1;
        global_info.tInf_pInf = 1;
        MinfiringS=0;
        tstf_counterS=0;
    else
        error("logic loop broken")
    end

    end
   

    % global_info.pInf_tInf = ceil(beta * m0Susceptible * m0Infected);
    if round(gamma * m0Infected)==0 && (gamma * m0Infected)~=0 && MinfiringI==0 
        timestepstofireI=round(1/(gamma * m0Infected));
        tstf_counterI=tstf_counterI+1;
        MinfiringI=1;
        global_info.pInf_tRec = round(gamma * m0Infected);
        global_info.tRec_pRec = round(gamma * m0Infected);
    elseif MinfiringI==0
        global_info.pInf_tRec = round(gamma * m0Infected);
        global_info.tRec_pRec = round(gamma * m0Infected);
    elseif MinfiringI==1 && tstf_counterI<timestepstofireI
        tstf_counterI=tstf_counterI+1;
        global_info.pInf_tRec = 0;
        global_info.tRec_pRec = 0;
    elseif MinfiringI==1 && tstf_counterI==timestepstofireI
        global_info.pInf_tRec =1;
        global_info.tRec_pRec = 1;
        MinfiringI=0;
        tstf_counterI=0;
    else
        error("logic loop broken")
    end


    if round(delta * m0Recovered)==0 && (delta * m0Recovered)~=0 && MinfiringR==0 
        timestepstofireR=round(1/(delta * m0Recovered));
        tstf_counterR=tstf_counterR+1;
        MinfiringR=1;
        global_info.pRec_tSus = round(delta * m0Recovered);
        global_info.tSus_pSus = round(delta * m0Recovered);
    elseif MinfiringR==0
        global_info.pRec_tSus = round(delta * m0Recovered);
        global_info.tSus_pSus = round(delta * m0Recovered);
    elseif MinfiringR==1 && tstf_counterR<timestepstofireR
        tstf_counterR=tstf_counterR+1;
        global_info.pRec_tSus = 0;
        global_info.tSus_pSus = 0;
    elseif MinfiringR==1 && tstf_counterI==timestepstofireI
        global_info.pRec_tSus = 1;
        global_info.tSus_pSus = 1;
        MinfiringR=0;
        tstf_counterR=0;
    else
        error("logic loop broken")
    end
    %}


    % if floor(beta * m0Susceptible * m0Infected)>=0
    % global_info.pSus_tInf = floor(beta * m0Susceptible * m0Infected);
    % global_info.tInf_pInf = floor(beta * m0Susceptible * m0Infected);
    % else
    % global_info.pSus_tInf = 0;
    % global_info.tInf_pInf = 0;
    % end
    % if floor(gamma * m0Infected)>=0
    % global_info.pInf_tRec = floor(gamma * m0Infected);
    % global_info.tRec_pRec = floor(gamma * m0Infected);
    % else
    % global_info.pInf_tRec = 0;
    % global_info.tRec_pRec = 0;
    % end
    % if floor(delta * m0Recovered)>=0
    % global_info.pRec_tSus = floor(delta * m0Recovered);
    % global_info.tSus_pSus = floor(delta * m0Recovered);
    % else
    % global_info.pRec_tSus = 0;
    % global_info.tSus_pSus = 0;    
    % end
    

    pns = pnstruct('simple_pn_SIRv1_pdf');
    dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    %figure()
    %plotp(Sim_Results, {'pSusceptible','pInfected','pRecovered'});

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
catch

    disp(['Error at delta=', int2str(delta), ' , gamma=',int2str(gamma), ' , beta=',int2str(beta), ' , timesteps=',int2str(timidivi)]);
    
    if delta==0
    delta=delta; 
    else
    delta=delta-10^-7; 
    end
   
    if beta==0
    beta=beta; 
    else
    beta=beta-10^-7;
    end
    
    if gamma==0
    gamma=gamma;
    else
    gamma=gamma-10^-7;
    end
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

%% Max iterations for Petri Net and ODE, petri net does 1 step per iteration
%in base settings
MAX_ITERATIONS = 100*timidivi;

%%%%% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

%% Initial Population
S_0=1000;
I_0=10;
R_0=10;

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

% the loop   %%%%%%%%%%%%%%%%%%%%%%%%%
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;
while le(Number_of_iterations, MAX_ITERATIONS)
    %disp(['Iteration: ', int2str(Number_of_iterations)]);
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 

    if beta_T==0 || m0Infected==0
        global_info.pSus_tInf = zero;
        pSus_tInf_resEr=(beta_T * m0Susceptible *  m0Infected + pSus_tInf_resEr)-round(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
        global_info.tInf_pInf = zero;
        tInf_pInf_resEr=(beta_T * m0Susceptible *  m0Infected+tInf_pInf_resEr)-round(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
    elseif m0Infected<=(1/(beta_T))
        global_info.pSus_tInf = round(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr);
        pSus_tInf_resEr=(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr)-round(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
        global_info.tInf_pInf = round(beta_T * m0Susceptible *   m0Infected + tInf_pInf_resEr);
        tInf_pInf_resEr=(beta_T * m0Susceptible *   m0Infected +tInf_pInf_resEr)-round(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
    else 

        %beta2=(1/(m0Infected+??????));
        beta2=1/(m0Infected+1);
        global_info.pSus_tInf = round(beta2 * m0Susceptible * m0Infected + pSus_tInf_resEr);
        pSus_tInf_resEr=(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr)-round(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr);
        global_info.tInf_pInf = round(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);
        tInf_pInf_resEr=(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr)-round(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);
    
    end

    global_info.pInf_tRec = round(gamma_T * (m0Infected) + pInf_tRec_resEr);
    pInf_tRec_resEr=(gamma_T * (m0Infected) + pInf_tRec_resEr)-round(gamma_T * (m0Infected) + pInf_tRec_resEr);
    global_info.tRec_pRec = round(gamma_T * (m0Infected) + tRec_pRec_resEr);
    tRec_pRec_resEr=(gamma_T * (m0Infected) + tRec_pRec_resEr)-round(gamma_T * (m0Infected) + tRec_pRec_resEr);
    global_info.pRec_tSus = round(delta_T * (m0Recovered) + pRec_tSus_resEr);
    pRec_tSus_resEr=(delta_T * (m0Recovered) + pRec_tSus_resEr)-round(delta_T * (m0Recovered) + pRec_tSus_resEr);
    global_info.tSus_pSus = round(delta_T * (m0Recovered) + tSus_pSus_resEr);
    tSus_pSus_resEr=(delta_T * (m0Recovered)+ tSus_pSus_resEr)-round(delta_T * (m0Recovered) + tSus_pSus_resEr);
    

    pns = pnstruct('simple_pn_SIRv1_pdf');
    %dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    %figure()
    %plotp(Sim_Results, {'pSusceptible','pInfected','pRecovered'});

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
end


MS_ER=MSusceptible;
MI_ER=MInfected;
MR_ER=MRecovered;


    elseif Rounder_test==2 %Floor
try
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;
%delta=1; beta=1; gamma=1;

%% Max iterations for Petri Net and ODE, petri net does 1 step per iteration
%in base settings
MAX_ITERATIONS = 100*timidivi;

%%%%% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 


%% Parameter values  %%%%%%%%%%%%%%%%%%%%%
%Used for setting specific parameter values (vs loop above for grid of
%param values)

%beta = 0.001;     gamma = 0.05;    delta = 0.001;
%beta = 0.1;     gamma = 0.2;    delta = 0.1;
%beta = 0.2;     gamma = 0.5;    delta = 0.4;
%beta = 3.1000e-04;     gamma =  0.0809;    
%delta = 0.00;


%% Initial Population
S_0=1000;
I_0=10;
R_0=10;

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

%% the loop   %%%%%%%%%%%%%%%%%%%%%%%%%
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;

while le(Number_of_iterations, MAX_ITERATIONS)
    %disp(['Iteration: ', int2str(Number_of_iterations)]);
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 


    if floor(beta * m0Susceptible * m0Infected)>=0
    global_info.pSus_tInf = floor(beta * m0Susceptible * m0Infected);
    global_info.tInf_pInf = floor(beta * m0Susceptible * m0Infected);
    else
    global_info.pSus_tInf = 0;
    global_info.tInf_pInf = 0;
    end
    if floor(gamma * m0Infected)>=0
    global_info.pInf_tRec = floor(gamma * m0Infected);
    global_info.tRec_pRec = floor(gamma * m0Infected);
    else
    global_info.pInf_tRec = 0;
    global_info.tRec_pRec = 0;
    end
    if floor(delta * m0Recovered)>=0
    global_info.pRec_tSus = floor(delta * m0Recovered);
    global_info.tSus_pSus = floor(delta * m0Recovered);
    else
    global_info.pRec_tSus = 0;
    global_info.tSus_pSus = 0;    
    end
    

    pns = pnstruct('simple_pn_SIRv1_pdf');
    dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    %figure()
    %plotp(Sim_Results, {'pSusceptible','pInfected','pRecovered'});

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
catch  % Floor Catch

    disp(['Error at delta=', int2str(delta), ' , gamma=',int2str(gamma), ' , beta=',int2str(beta), ' , timesteps=',int2str(timidivi)]);
    
    if delta==0
    delta=delta; 
    else
    delta=delta-10^-7; 
    end
   
    if beta==0
    beta=beta; 
    else
    beta=beta-10^-7;
    end
    
    if gamma==0
    gamma=gamma;
    else
    gamma=gamma-10^-7;
    end
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

% Max iterations for Petri Net and ODE, petri net does 1 step per iteration
%in base settings
MAX_ITERATIONS = 100*timidivi;

%%%%% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

%Initial Population
S_0=1000;
I_0=10;
R_0=10;
m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

%% the loop   %%%%%%%%%%%%%%%%%%%%%%%%%
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;
while le(Number_of_iterations, MAX_ITERATIONS)
    %disp(['Iteration: ', int2str(Number_of_iterations)]);
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 

    if floor(beta * m0Susceptible * m0Infected)>=0
    global_info.pSus_tInf = floor(beta * m0Susceptible * m0Infected);
    global_info.tInf_pInf = floor(beta * m0Susceptible * m0Infected);
    else
    global_info.pSus_tInf = 0;
    global_info.tInf_pInf = 0;
    end
    if floor(gamma * m0Infected)>=0
    global_info.pInf_tRec = floor(gamma * m0Infected);
    global_info.tRec_pRec = floor(gamma * m0Infected);
    else
    global_info.pInf_tRec = 0;
    global_info.tRec_pRec = 0;
    end
    if floor(delta * m0Recovered)>=0
    global_info.pRec_tSus = floor(delta * m0Recovered);
    global_info.tSus_pSus = floor(delta * m0Recovered);
    else
    global_info.pRec_tSus = 0;
    global_info.tSus_pSus = 0;    
    end

    pns = pnstruct('simple_pn_SIRv1_pdf');
    %dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    %figure()
    %plotp(Sim_Results, {'pSusceptible','pInfected','pRecovered'});

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
end
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MSusceptible, 'r--o','LineWidth',2);
MS_Floor=MSusceptible;
MI_Floor=MInfected;
MR_Floor=MRecovered;
    


%%CEILING
elseif Rounder_test==3 %Ceil
try
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;
%delta=1; beta=1; gamma=1;

%% Max iterations for Petri Net and ODE, petri net does 1 step per iteration
%in base settings
MAX_ITERATIONS = 100*timidivi;

%%%%% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 


%% Parameter values  %%%%%%%%%%%%%%%%%%%%%
%Used for setting specific parameter values (vs loop above for grid of
%param values)

%beta = 0.001;     gamma = 0.05;    delta = 0.001;

%% Initial Population
S_0=1000;
I_0=10;
R_0=10;

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

%% the loop   %%%%%%%%%%%%%%%%%%%%%%%%%
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;

while le(Number_of_iterations, MAX_ITERATIONS)
    %disp(['Iteration: ', int2str(Number_of_iterations)]);
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 
   


    global_info.pSus_tInf = ceil(beta * m0Susceptible * m0Infected);
    global_info.tInf_pInf = ceil(beta * m0Susceptible * m0Infected);
   
   
    global_info.pInf_tRec = ceil(gamma * m0Infected);
    global_info.tRec_pRec = ceil(gamma * m0Infected);
   
   
    global_info.pRec_tSus = ceil(delta * m0Recovered);
    global_info.tSus_pSus = ceil(delta * m0Recovered);
    

    pns = pnstruct('simple_pn_SIRv1_pdf');
    dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    %figure()
    %plotp(Sim_Results, {'pSusceptible','pInfected','pRecovered'});

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
catch

    disp(['Error at delta=', int2str(delta), ' , gamma=',int2str(gamma), ' , beta=',int2str(beta), ' , timesteps=',int2str(timidivi)]);
    
    if delta==0
    delta=delta; 
    else
    delta=delta-10^-7; 
    end
   
    if beta==0
    beta=beta; 
    else
    beta=beta-10^-7;
    end
    
    if gamma==0
    gamma=gamma;
    else
    gamma=gamma-10^-7;
    end
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

%% Max iterations for Petri Net and ODE, petri net does 1 step per iteration
%in base settings
MAX_ITERATIONS = 100*timidivi;

%%%%% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

%% Initial Population
S_0=1000;
I_0=10;
R_0=10;

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

%% the loop   %%%%%%%%%%%%%%%%%%%%%%%%%
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;
while le(Number_of_iterations, MAX_ITERATIONS)
    %disp(['Iteration: ', int2str(Number_of_iterations)]);
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 

   global_info.pSus_tInf = ceil(beta * m0Susceptible * m0Infected);
    global_info.tInf_pInf = ceil(beta * m0Susceptible * m0Infected);
   
   
    global_info.pInf_tRec = ceil(gamma * m0Infected);
    global_info.tRec_pRec = ceil(gamma * m0Infected);
   
   
    global_info.pRec_tSus = ceil(delta * m0Recovered);
    global_info.tSus_pSus = ceil(delta * m0Recovered);

    pns = pnstruct('simple_pn_SIRv1_pdf');
    %dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    %figure()
    %plotp(Sim_Results, {'pSusceptible','pInfected','pRecovered'});

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
end
MS_Ceil=MSusceptible;
MI_Ceil=MInfected;
MR_Ceil=MRecovered;




%%standar round
elseif Rounder_test==4 %Round standard
try
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;
%delta=1; beta=1; gamma=1;

%% Max iterations for Petri Net and ODE, petri net does 1 step per iteration
%in base settings
MAX_ITERATIONS = 100*timidivi;

%%%%% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 


%% Parameter values  %%%%%%%%%%%%%%%%%%%%%
%Used for setting specific parameter values (vs loop above for grid of
%param values)

%beta = 0.001;     gamma = 0.05;    delta = 0.001;

%% Initial Population
S_0=1000;
I_0=10;
R_0=10;

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

%% the loop   %%%%%%%%%%%%%%%%%%%%%%%%%
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;

while le(Number_of_iterations, MAX_ITERATIONS)
    %disp(['Iteration: ', int2str(Number_of_iterations)]);
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered;

    global_info.pSus_tInf = round(beta * m0Susceptible * m0Infected);
    global_info.tInf_pInf = round(beta * m0Susceptible * m0Infected);
   
   
    global_info.pInf_tRec = round(gamma * m0Infected);
    global_info.tRec_pRec = ceil(gamma * m0Infected);
   
   
    global_info.pRec_tSus = round(delta * m0Recovered);
    global_info.tSus_pSus = round(delta * m0Recovered);
    

    pns = pnstruct('simple_pn_SIRv1_pdf');
    dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    %figure()
    %plotp(Sim_Results, {'pSusceptible','pInfected','pRecovered'});

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
catch

    disp(['Error at delta=', int2str(delta), ' , gamma=',int2str(gamma), ' , beta=',int2str(beta), ' , timesteps=',int2str(timidivi)]);
    
    if delta==0
    delta=delta; 
    else
    delta=delta-10^-7; 
    end
   
    if beta==0
    beta=beta; 
    else
    beta=beta-10^-7;
    end
    
    if gamma==0
    gamma=gamma;
    else
    gamma=gamma-10^-7;
    end
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

%% Max iterations for Petri Net and ODE, petri net does 1 step per iteration
%in base settings
MAX_ITERATIONS = 100*timidivi;

%%%%% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

%% Initial Population
S_0=1000;
I_0=10;
R_0=10;

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

%% the loop   %%%%%%%%%%%%%%%%%%%%%%%%%
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;
while le(Number_of_iterations, MAX_ITERATIONS)
    %disp(['Iteration: ', int2str(Number_of_iterations)]);
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 

   global_info.pSus_tInf = round(beta * m0Susceptible * m0Infected);
    global_info.tInf_pInf = round(beta * m0Susceptible * m0Infected);
   
   
    global_info.pInf_tRec = round(gamma * m0Infected);
    global_info.tRec_pRec = round(gamma * m0Infected);
   
   
    global_info.pRec_tSus = round(delta * m0Recovered);
    global_info.tSus_pSus = round(delta * m0Recovered);

    pns = pnstruct('simple_pn_SIRv1_pdf');
    %dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    %figure()
    %plotp(Sim_Results, {'pSusceptible','pInfected','pRecovered'});

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
end
MS_SR=MSusceptible;
MI_SR=MInfected;
MR_SR=MRecovered;
    end
end
if ode_run_counter==1
y0=[S_0;I_0;R_0];
delta_ODE=delta; 
beta_ODE=beta; 
gamma_ODE=gamma;
params=[delta_ODE;beta_ODE;gamma_ODE];
[t,y]=ode15s(@SIR,[1:1:(MAX_ITERATIONS/timidivi)],y0,[],params);
else
end

figure()
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MSusceptible, 'b--o','LineWidth',2);
hold on
plot(t,y(:,1), 'y','LineWidth',2)
hold off
%% RRMSE
 err_S = (rmse(MSusceptible(1:timidivi:MAX_ITERATIONS) , y(:,1)')/sqrt(sumsqr(y(:,1)')))*100;
 err_I = (rmse(MInfected(1:timidivi:MAX_ITERATIONS) , y(:,2)')/sqrt(sumsqr( y(:,2)')))*100;
 err_R = (rmse(MRecovered(1:timidivi:MAX_ITERATIONS) , y(:,3)')/sqrt(sumsqr(y(:,3)')))*100;
 % % 
 % % timidivi_col_S=[timidivi_col_S err_S];
 % % timidivi_col_I=[timidivi_col_I err_I];
 % % timidivi_col_R=[timidivi_col_R err_R];
% %end
 %% Plot times steps vs rrmse %%%%%%%%%%%%%%%%%%
% figure()
% plot(timidivi_timestep, timidivi_col_S, 'b--o','LineWidth',2);
% hold on
% plot(timidivi_timestep, timidivi_col_I,'c--o','LineWidth',2);
% plot(timidivi_timestep, timidivi_col_R, 'r--o','LineWidth',2);
% legend('Susceptible_{RRMSE}','Infected_{RRMS}','Recovered_{RRMSE}');
% title("SIR Petri Net vs ODE RRMSE based on times steps", 'FontSize', 24)
% subtitle("beta="+beta*timidivi+', gamma='+gamma*timidivi+' delta='+delta*timidivi)
% xlabel('PN Time steps (per 1 ODE time interval)'); ylabel('RRMSE');
% hold off

figure()
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MS_Ceil, 'color', [.15 .15 .15],'LineWidth',12);
hold on
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MS_Floor,'color', [.65 .65 .65],'LineWidth',8);
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MS_SR, 'color', [.90 .90 .90],'LineWidth',4);
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MS_ER, 'b','LineWidth',4);
plot(t,y(:,1), 'r','LineWidth',2); 
lgd=legend('Ceiling','Floor','Standard','Standard+Residuals','ODE');legend('Location','best');
lgd.Title.String = 'Arc Weight Rounding Method';
title("SIR PN Rounding Methods vs ODE Susceptible", 'FontSize', 24)
subtitle("gamma (\gamma )="+gamma+", beta (\beta)="+beta+", delta (\delta)="+delta)
xlabel('Time')
ylabel('Population') 
xlim([0,105])
hold off

figure()
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MI_Ceil,'color', [.15 .15 .15],'LineWidth',12);
hold on
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MI_Floor,'color', [.65 .65 .65],'LineWidth',8);
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MI_SR, 'color', [.90 .90 .90],'LineWidth',4);
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MI_ER, 'b','LineWidth',4);
plot(t,y(:,2), 'r','LineWidth',2);
lgd=legend('Ceiling','Floor','Standard','Standard+Residuals','ODE');legend('Location','best');
lgd.Title.String = 'Arc Weight Rounding Method';
title("SIR PN Rounding Methods vs ODE Infected", 'FontSize', 24)
subtitle("gamma (\gamma )="+gamma+", beta (\beta)="+beta+", delta (\delta)="+delta)
xlabel('Time')
ylabel('Population') 
xlim([0,105])
hold off

figure()
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MR_Ceil,'color', [.15 .15 .15],'LineWidth',12);
hold on
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MR_Floor,'color', [.65 .65 .65],'LineWidth',8);
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MR_SR, 'color', [.90 .90 .90],'LineWidth',4);
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MR_ER, 'b','LineWidth',4);
plot(t,y(:,3), 'r','LineWidth',2);  
lgd=legend('Ceiling','Floor','Standard','Standard+Residuals','ODE');legend('Location','best');
lgd.Title.String = 'Arc Weight Rounding Method';
title("SIR PN Rounding Methods vs ODE Recovered", 'FontSize', 24)
subtitle("gamma (\gamma )="+gamma+", beta (\beta)="+beta+", delta (\delta)="+delta)
xlabel('Time')
ylabel('Population') 
xlim([0,105])
hold off
% % disp(['still going: ', int2str(deltactr)]);

%% Plot the results of single param value run vs ODE %%%%%%%%%%%%%%%%%%
 figure()
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MSusceptible, 'b--o','LineWidth',2);
hold on
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MInfected,'c--o','LineWidth',2);
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MRecovered, 'r--o','LineWidth',2);

plot(t,y(:,1), 'y','LineWidth',2)
plot(t,y(:,2),'m','LineWidth',2)
plot(t,y(:,3), 'k','LineWidth',2);
legend('Susceptible_{PN}','Infected_{PN}','Recovered_{PN}',...
    'Susceptible_{ODE}','Infected_{ODE}','Recovered_{ODE}');
title("SIR Petri Net vs ODE", 'FontSize', 24)
subtitle("beta="+beta*timidivi+', gamma='+gamma*timidivi+' delta='+delta*timidivi+", Susc. MSE="+err_S...
   +", Infe. MSE="+err_I+", Reco. MSE="+err_R, 'FontSize', 14)
hold off
% 
% figure()
% plot([1,1+1/timidivi,2:1:100], MSusceptible, 'b--o','LineWidth',2);
% hold on
% plot([1,1+1/timidivi,2:1:100], MInfected,'c--o','LineWidth',2);
% plot([1,1+1/timidivi,2:1:100], MRecovered, 'r--o','LineWidth',2);
% 
% plot(t,y(:,1), 'y','LineWidth',2)
% plot(t,y(:,2),'m','LineWidth',2)
% plot(t,y(:,3), 'k','LineWidth',2);
% legend('Susceptible_{PN}','Infected_{PN}','Recovered_{PN}',...
%     'Susceptible_{ODE}','Infected_{ODE}','Recovered_{ODE}');
% title("SIR Petri Net vs ODE", 'FontSize', 24)
% subtitle("beta="+beta*timidivi+', gamma='+gamma*timidivi+' delta='+delta*timidivi+", Susc. MSE="+err_S...
%    +", Infe. MSE="+err_I+", Reco. MSE="+err_R, 'FontSize', 14)
% hold off


        % MSE_s(gammactr,betactr)=err_S;
        % MSE_i(gammactr,betactr)=err_I;
        % MSE_iBeta(gammactr,betactr)=beta;
        % MSE_iGamma(gammactr,betactr)=gamma;
        % MSE_r(gammactr,betactr)=err_R;

        MSE_s(gammactr,betactr,deltactr)=err_S;
        MSE_i(gammactr,betactr,deltactr)=err_I;
        MSE_r(gammactr,betactr,deltactr)=err_R;
        % 
% MSE_S_c(deltactr)=zeros(1,5);
% MSE_S_r(deltactr)=zeros(1,5);
% MSE_S_v(deltactr)=zeros(1,5);
% MSE_S_C(dr)=[MSE_S_C(dr) err_S];
% % 
  % % % % %          gammactr=gammactr+1;
  % % % % %        end
  % % % % %       betactr=betactr+1;
  % % % % %    end
  % % % % %    deltactr=deltactr+1;
  % % % % % end
  
  
%   tau_S=[tau_S (sum(MSE_s,"all")/bglength^3)];
%   tau_I=[tau_I (sum(MSE_i,"all")/bglength^3)];
%   tau_R=[tau_R (sum(MSE_r,"all")/bglength^3)];
% end
% figure()
% plot(tau_spread, tau_S, 'b--o','LineWidth',2);
% hold on
% plot(tau_spread, tau_I, 'c--o','LineWidth',2);
% plot(tau_spread, tau_R, 'r--o','LineWidth',2);
% legend('Susceptible Mean RRMSE','Infected Mean RRMSE','Recovered Mean RRMSE');
% title("Mean Error for different Parameter ranges (0,\tau)", 'FontSize', 24)
% xlabel("\tau (parameters \gamma,\beta,\delta range [0,\tau])");ylabel('RRMSE Mean')
% subtitle("beta="+beta+', gamma='+gamma+' delta='+delta+", Susc. MSE="+err_S...
%     +", Infe. MSE="+err_I+", Reco. MSE="+err_R, 'FontSize', 14)
%  hold off


% minValue = min(MSE_i(:));
% % Find all (row, column) pairs where M = the min value.
% [rows, columns] = find(MSE_i == minValue)
% % Print them out:
% for k = 1 : length(rows)
% 	fprintf('M equals the min value of %f at row = %d, column = %d.\n', ...
% 		minValue, rows(k), columns(k));
% end
% MSE_iBeta(rows(k),columns(k))
%  MSE_iGamma(rows(k),columns(k))

%% Plot the results of grid param value RRMSE vs ODE
% casxismin=min([min(MSE_s(:)),min(MSE_i(:)),min(MSE_r(:))]);
% trial=1;
% if trial==1
% casxismax=max([max(MSE_s(:)),max(MSE_i(:)),max(MSE_r(:))]);
% trial=trial+1;
% else
% end

casxismin=0;
casxismax=44;
Cus_map=[(ones(2,3).*[0 0.5 0]);(ones(8,3).*[0 1 0]);(ones(40,3).*[1 1 0]);(ones(50,3).*[1 0.75 0]);...
    (ones(100,3).*[1 0.4 0]);(ones(240,3).*[1 0 0])];


figure()
%title("RRMSE for parameter Ranges with 1 pn step per 1 ode unit time", 'FontSize', 14)
tiledlayout(5,3);
nexttile
s=surf(spacer,spacer,MSE_s(:,:,1));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Susceptible Delta=0')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h, 'RRMSE (%)');
%h.Label.String = "RRMSE (%)";h.Label.Rotation = 0;
view(2)

nexttile
s=surf(spacer,spacer,MSE_i(:,:,1));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Infected Delta=0')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h, 'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r(:,:,1));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Recovered Delta=0')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h, 'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_s(:,:,2));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Susceptible Delta=0.001')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i(:,:,2));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Infected Delta=0.001')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r(:,:,2));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Recovered Delta=0.001')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_s(:,:,3));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Susceptible Delta=0.01')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i(:,:,3));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Infected Delta=0.01')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r(:,:,3));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Recovered Delta=0.01')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_s(:,:,4));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Susceptible Delta=0.1')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i(:,:,4));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Infected Delta=0.1')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r(:,:,4));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Recovered Delta=0.1')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_s(:,:,5));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Susceptible Delta=1')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i(:,:,5));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Infected Delta=1')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r(:,:,5));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Recovered Delta=1')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
view(2)
h=colorbar('Position',[0.931120331950207,0.11,0.014914910340541,0.816052332195677]);
title(h, 'RRMSE (%)');
%h=axes('Position', [0.85 0.15 0.05 0.7]);% axes('Position', [0.85 0.15 0.05 0.7]),   %'Location', 'westoutside',
%h.Label.String = "RRMSE (%)";h.Label.Rotation = 0;

MSE_s
MSE_i
MSE_r
ode_run_counter=ode_run_counter+1;
%end

end %end of Run Section 3














%% Run Section 4, comparing different types of rounding methods for single param value
if Run_Section==4

timidivi=20;                %tau in paper

  bglength=1;

spacer=linspace(0,1,bglength);

MSE_s=[zeros(bglength,bglength,5)];
MSE_i=[zeros(bglength,bglength,5)];
MSE_r=[zeros(bglength,bglength,5)];


delta=1; beta=1; gamma=1;   %parameter values

% Initial Population
S_0=1000;
I_0=10;
R_0=10;

MAX_ITERATIONS = 100*timidivi;  % Max iterations for Petri Net and ODE

for Rounder_test=1:1:4      %loop for each rounding methoded tested
    
    if Rounder_test==1      % Rounding method standard+residual
try
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

% the loop 
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;

while le(Number_of_iterations, MAX_ITERATIONS)
    
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 
    
    if beta_T==0 || m0Infected==0
        global_info.pSus_tInf = zero;
        pSus_tInf_resEr=(beta_T * m0Susceptible *  m0Infected + pSus_tInf_resEr)-round(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
        global_info.tInf_pInf = zero;
        tInf_pInf_resEr=(beta_T * m0Susceptible *  m0Infected+tInf_pInf_resEr)-round(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
    elseif m0Infected<=(1/(beta_T))
        global_info.pSus_tInf = round(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr);
        pSus_tInf_resEr=(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr)-round(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
        global_info.tInf_pInf = round(beta_T * m0Susceptible *   m0Infected + tInf_pInf_resEr);
        tInf_pInf_resEr=(beta_T * m0Susceptible *   m0Infected +tInf_pInf_resEr)-round(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
    else 

        beta2=1/(m0Infected+1);
        global_info.pSus_tInf = round(beta2 * m0Susceptible * m0Infected + pSus_tInf_resEr);
        pSus_tInf_resEr=(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr)-round(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr);
        global_info.tInf_pInf = round(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);
        tInf_pInf_resEr=(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr)-round(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);
    
    end

    global_info.pInf_tRec = round(gamma_T * (m0Infected) + pInf_tRec_resEr);
    pInf_tRec_resEr=(gamma_T * (m0Infected) + pInf_tRec_resEr)-round(gamma_T * (m0Infected) + pInf_tRec_resEr);
    global_info.tRec_pRec = round(gamma_T * (m0Infected) + tRec_pRec_resEr);
    tRec_pRec_resEr=(gamma_T * (m0Infected) + tRec_pRec_resEr)-round(gamma_T * (m0Infected) + tRec_pRec_resEr);
    global_info.pRec_tSus = round(delta_T * (m0Recovered) + pRec_tSus_resEr);
    pRec_tSus_resEr=(delta_T * (m0Recovered) + pRec_tSus_resEr)-round(delta_T * (m0Recovered) + pRec_tSus_resEr);
    global_info.tSus_pSus = round(delta_T * (m0Recovered) + tSus_pSus_resEr);
    tSus_pSus_resEr=(delta_T * (m0Recovered)+ tSus_pSus_resEr)-round(delta_T * (m0Recovered) + tSus_pSus_resEr);
    
    pns = pnstruct('simple_pn_SIRv1_pdf');
    dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
catch % Rounding method standard+residual

    disp(['Error at delta=', int2str(delta), ' , gamma=',int2str(gamma), ' , beta=',int2str(beta), ' , timesteps=',int2str(timidivi)]);
    
    if delta==0
    delta=delta; 
    else
    delta=delta-10^-7; 
    end
   
    if beta==0
    beta=beta; 
    else
    beta=beta-10^-7;
    end
    
    if gamma==0
    gamma=gamma;
    else
    gamma=gamma-10^-7;
    end
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

%Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

% the loop   
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;
while le(Number_of_iterations, MAX_ITERATIONS)
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 

    if beta_T==0 || m0Infected==0
        global_info.pSus_tInf = zero;
        pSus_tInf_resEr=(beta_T * m0Susceptible *  m0Infected + pSus_tInf_resEr)-round(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
        global_info.tInf_pInf = zero;
        tInf_pInf_resEr=(beta_T * m0Susceptible *  m0Infected+tInf_pInf_resEr)-round(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
    elseif m0Infected<=(1/(beta_T))
        global_info.pSus_tInf = round(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr);
        pSus_tInf_resEr=(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr)-round(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
        global_info.tInf_pInf = round(beta_T * m0Susceptible *   m0Infected + tInf_pInf_resEr);
        tInf_pInf_resEr=(beta_T * m0Susceptible *   m0Infected +tInf_pInf_resEr)-round(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
    else 

        beta2=1/(m0Infected+1);
        global_info.pSus_tInf = round(beta2 * m0Susceptible * m0Infected + pSus_tInf_resEr);
        pSus_tInf_resEr=(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr)-round(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr);
        global_info.tInf_pInf = round(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);
        tInf_pInf_resEr=(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr)-round(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);
    
    end

    global_info.pInf_tRec = round(gamma_T * (m0Infected) + pInf_tRec_resEr);
    pInf_tRec_resEr=(gamma_T * (m0Infected) + pInf_tRec_resEr)-round(gamma_T * (m0Infected) + pInf_tRec_resEr);
    global_info.tRec_pRec = round(gamma_T * (m0Infected) + tRec_pRec_resEr);
    tRec_pRec_resEr=(gamma_T * (m0Infected) + tRec_pRec_resEr)-round(gamma_T * (m0Infected) + tRec_pRec_resEr);
    global_info.pRec_tSus = round(delta_T * (m0Recovered) + pRec_tSus_resEr);
    pRec_tSus_resEr=(delta_T * (m0Recovered) + pRec_tSus_resEr)-round(delta_T * (m0Recovered) + pRec_tSus_resEr);
    global_info.tSus_pSus = round(delta_T * (m0Recovered) + tSus_pSus_resEr);
    tSus_pSus_resEr=(delta_T * (m0Recovered)+ tSus_pSus_resEr)-round(delta_T * (m0Recovered) + tSus_pSus_resEr);
    

    pns = pnstruct('simple_pn_SIRv1_pdf');
    %dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    
    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
end

%store values of standard+residual run
MS_ER=MSusceptible;
MI_ER=MInfected;
MR_ER=MRecovered;


    elseif Rounder_test==2 % Rounding method Floor
try
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

%the loop
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;

while le(Number_of_iterations, MAX_ITERATIONS)
    
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 

    if floor(beta * m0Susceptible * m0Infected)>=0
    global_info.pSus_tInf = floor(beta * m0Susceptible * m0Infected);
    global_info.tInf_pInf = floor(beta * m0Susceptible * m0Infected);
    else
    global_info.pSus_tInf = 0;
    global_info.tInf_pInf = 0;
    end
    if floor(gamma * m0Infected)>=0
    global_info.pInf_tRec = floor(gamma * m0Infected);
    global_info.tRec_pRec = floor(gamma * m0Infected);
    else
    global_info.pInf_tRec = 0;
    global_info.tRec_pRec = 0;
    end
    if floor(delta * m0Recovered)>=0
    global_info.pRec_tSus = floor(delta * m0Recovered);
    global_info.tSus_pSus = floor(delta * m0Recovered);
    else
    global_info.pRec_tSus = 0;
    global_info.tSus_pSus = 0;    
    end  

    pns = pnstruct('simple_pn_SIRv1_pdf');
    dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
catch  % Rounding method Floor Catch

    disp(['Error at delta=', int2str(delta), ' , gamma=',int2str(gamma), ' , beta=',int2str(beta), ' , timesteps=',int2str(timidivi)]);
    
    if delta==0
    delta=delta; 
    else
    delta=delta-10^-7; 
    end
   
    if beta==0
    beta=beta; 
    else
    beta=beta-10^-7;
    end
    
    if gamma==0
    gamma=gamma;
    else
    gamma=gamma-10^-7;
    end
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

% the loop   
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;
while le(Number_of_iterations, MAX_ITERATIONS)
   
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 

    if floor(beta * m0Susceptible * m0Infected)>=0
    global_info.pSus_tInf = floor(beta * m0Susceptible * m0Infected);
    global_info.tInf_pInf = floor(beta * m0Susceptible * m0Infected);
    else
    global_info.pSus_tInf = 0;
    global_info.tInf_pInf = 0;
    end
    if floor(gamma * m0Infected)>=0
    global_info.pInf_tRec = floor(gamma * m0Infected);
    global_info.tRec_pRec = floor(gamma * m0Infected);
    else
    global_info.pInf_tRec = 0;
    global_info.tRec_pRec = 0;
    end
    if floor(delta * m0Recovered)>=0
    global_info.pRec_tSus = floor(delta * m0Recovered);
    global_info.tSus_pSus = floor(delta * m0Recovered);
    else
    global_info.pRec_tSus = 0;
    global_info.tSus_pSus = 0;    
    end

    pns = pnstruct('simple_pn_SIRv1_pdf');
    %dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
end
%plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MSusceptible, 'r--o','LineWidth',2);
%store values of floor function
MS_Floor=MSusceptible;
MI_Floor=MInfected;
MR_Floor=MRecovered;
    


%%CEILING
elseif Rounder_test==3 % Rounding method Ceiling

    try
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

%the loop   %%%%%%%%%%%%%%%%%%%%%%%%%
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;

while le(Number_of_iterations, MAX_ITERATIONS)
    
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 

    global_info.pSus_tInf = ceil(beta * m0Susceptible * m0Infected);
    global_info.tInf_pInf = ceil(beta * m0Susceptible * m0Infected);
   
    global_info.pInf_tRec = ceil(gamma * m0Infected);
    global_info.tRec_pRec = ceil(gamma * m0Infected);
   
    global_info.pRec_tSus = ceil(delta * m0Recovered);
    global_info.tSus_pSus = ceil(delta * m0Recovered);
    
    pns = pnstruct('simple_pn_SIRv1_pdf');
    dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    %figure()
    %plotp(Sim_Results, {'pSusceptible','pInfected','pRecovered'});

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
    catch % Rounding method Ceiling

    disp(['Error at delta=', int2str(delta), ' , gamma=',int2str(gamma), ' , beta=',int2str(beta), ' , timesteps=',int2str(timidivi)]);
    
    if delta==0
    delta=delta; 
    else
    delta=delta-10^-7; 
    end
   
    if beta==0
    beta=beta; 
    else
    beta=beta-10^-7;
    end
    
    if gamma==0
    gamma=gamma;
    else
    gamma=gamma-10^-7;
    end
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

%Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

%the loop 
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;
while le(Number_of_iterations, MAX_ITERATIONS)
    %disp(['Iteration: ', int2str(Number_of_iterations)]);
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 

   global_info.pSus_tInf = ceil(beta * m0Susceptible * m0Infected);
    global_info.tInf_pInf = ceil(beta * m0Susceptible * m0Infected);
  
    global_info.pInf_tRec = ceil(gamma * m0Infected);
    global_info.tRec_pRec = ceil(gamma * m0Infected);
   
    global_info.pRec_tSus = ceil(delta * m0Recovered);
    global_info.tSus_pSus = ceil(delta * m0Recovered);

    pns = pnstruct('simple_pn_SIRv1_pdf');
    %dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    %figure()
    %plotp(Sim_Results, {'pSusceptible','pInfected','pRecovered'});

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
    end
    %store value from ceiling function
MS_Ceil=MSusceptible;
MI_Ceil=MInfected;
MR_Ceil=MRecovered;




%standard round
elseif Rounder_test==4 % Rounding method standard

try
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

% the loop  
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;

while le(Number_of_iterations, MAX_ITERATIONS)
    
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered;

    global_info.pSus_tInf = round(beta * m0Susceptible * m0Infected);
    global_info.tInf_pInf = round(beta * m0Susceptible * m0Infected);
   
    global_info.pInf_tRec = round(gamma * m0Infected);
    global_info.tRec_pRec = ceil(gamma * m0Infected);
   
    global_info.pRec_tSus = round(delta * m0Recovered);
    global_info.tSus_pSus = round(delta * m0Recovered);

    pns = pnstruct('simple_pn_SIRv1_pdf');
    dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    %figure()
    %plotp(Sim_Results, {'pSusceptible','pInfected','pRecovered'});

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
catch % Rounding method standard catch

    disp(['Error at delta=', int2str(delta), ' , gamma=',int2str(gamma), ' , beta=',int2str(beta), ' , timesteps=',int2str(timidivi)]);
    
    if delta==0
    delta=delta; 
    else
    delta=delta-10^-7; 
    end
   
    if beta==0
    beta=beta; 
    else
    beta=beta-10^-7;
    end
    
    if gamma==0
    gamma=gamma;
    else
    gamma=gamma-10^-7;
    end
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

% the loop 
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;
while le(Number_of_iterations, MAX_ITERATIONS)
    %disp(['Iteration: ', int2str(Number_of_iterations)]);
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 

    global_info.pSus_tInf = round(beta * m0Susceptible * m0Infected);
    global_info.tInf_pInf = round(beta * m0Susceptible * m0Infected);

    global_info.pInf_tRec = round(gamma * m0Infected);
    global_info.tRec_pRec = round(gamma * m0Infected);

    global_info.pRec_tSus = round(delta * m0Recovered);
    global_info.tSus_pSus = round(delta * m0Recovered);

    pns = pnstruct('simple_pn_SIRv1_pdf');
    %dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    %figure()
    %plotp(Sim_Results, {'pSusceptible','pInfected','pRecovered'});

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
end
%store standard rounding values
MS_SR=MSusceptible;
MI_SR=MInfected;
MR_SR=MRecovered;
    end
end

%Run ode model

y0=[S_0;I_0;R_0];
delta_ODE=delta; 
beta_ODE=beta; 
gamma_ODE=gamma;
params=[delta_ODE;beta_ODE;gamma_ODE];
[t,y]=ode15s(@SIR,[1:1:(MAX_ITERATIONS/timidivi)],y0,[],params);

% %% RRMSE
%  err_S = (rmse(MSusceptible(1:timidivi:MAX_ITERATIONS) , y(:,1)')/sqrt(sumsqr(y(:,1)')))*100;
%  err_I = (rmse(MInfected(1:timidivi:MAX_ITERATIONS) , y(:,2)')/sqrt(sumsqr( y(:,2)')))*100;
%  err_R = (rmse(MRecovered(1:timidivi:MAX_ITERATIONS) , y(:,3)')/sqrt(sumsqr(y(:,3)')))*100;

figure()
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MS_Ceil, 'color', [.15 .15 .15],'LineWidth',12);
hold on
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MS_Floor,'color', [.65 .65 .65],'LineWidth',8);
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MS_SR, 'color', [.90 .90 .90],'LineWidth',4);
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MS_ER, 'b','LineWidth',4);
plot(t,y(:,1), 'r','LineWidth',2); 
lgd=legend('Ceiling','Floor','Standard','Standard+Residuals','ODE');legend('Location','best');
lgd.Title.String = 'Arc Weight Rounding Method';
title("SIR PN Rounding Methods vs ODE Susceptible", 'FontSize', 24)
subtitle("gamma (\gamma )="+gamma+", beta (\beta)="+beta+", delta (\delta)="+delta)
xlabel('Time')
ylabel('Population') 
xlim([0,105])
hold off

figure()
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MI_Ceil,'color', [.15 .15 .15],'LineWidth',12);
hold on
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MI_Floor,'color', [.65 .65 .65],'LineWidth',8);
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MI_SR, 'color', [.90 .90 .90],'LineWidth',4);
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MI_ER, 'b','LineWidth',4);
plot(t,y(:,2), 'r','LineWidth',2);
lgd=legend('Ceiling','Floor','Standard','Standard+Residuals','ODE');legend('Location','best');
lgd.Title.String = 'Arc Weight Rounding Method';
title("SIR PN Rounding Methods vs ODE Infected", 'FontSize', 24)
subtitle("gamma (\gamma )="+gamma+", beta (\beta)="+beta+", delta (\delta)="+delta)
xlabel('Time')
ylabel('Population') 
xlim([0,105])
hold off

figure()
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MR_Ceil,'color', [.15 .15 .15],'LineWidth',12);
hold on
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MR_Floor,'color', [.65 .65 .65],'LineWidth',8);
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MR_SR, 'color', [.90 .90 .90],'LineWidth',4);
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MR_ER, 'b','LineWidth',4);
plot(t,y(:,3), 'r','LineWidth',2);  
lgd=legend('Ceiling','Floor','Standard','Standard+Residuals','ODE');legend('Location','best');
lgd.Title.String = 'Arc Weight Rounding Method';
title("SIR PN Rounding Methods vs ODE Recovered", 'FontSize', 24)
subtitle("gamma (\gamma )="+gamma+", beta (\beta)="+beta+", delta (\delta)="+delta)
xlabel('Time')
ylabel('Population') 
xlim([0,105])
hold off

end %end of Run Section 4



















%% Run Section 5, run SIRS PN model over range of parameters and Rounding methods and compares RRMSE
% not currently used in paper so code not updated
if Run_Section==5
timidivi_col_S=[];
timidivi_col_I=[];
timidivi_col_R=[];
%timidivi_timestep=1:6:79;
%for timidivi=1:6:79
ode_run_counter=1;
%timidivi=1;
%for timidivi=[2]
timidivi=20;
%% The following section is used for finding the RRMSE for various
%parameter values in comparison to the respective ODE, this is limited two
%comparison of 3 parameters in this configuration, for visual display.
  bglength=1;
% tau_S=[];
% tau_I=[];
% tau_R=[];
% tau_spread=linspace(0,1,10);

% for tau=linspace(0,1,10)
%spacer=[0 logspace(-5,0,bglength-1)];
%spacer=[0 logspace(-5,-1,bglength) linspace(0.2,0.8,bglength) (ones(1,bglength)-logspace(-1,-5,bglength)) 1];
% spacer=linspace(0,tau,bglength);
spacer=linspace(0,1,bglength);
% gridsize=size(spacer,2);
% 
% MSE_s=[zeros(gridsize,gridsize)];
% MSE_i=[zeros(gridsize,gridsize)];
% MSE_r=[zeros(gridsize,gridsize)];
MSE_s=[zeros(bglength,bglength,5)];
MSE_i=[zeros(bglength,bglength,5)];
MSE_r=[zeros(bglength,bglength,5)];

% MSE_S_c=zeros(1,5);
% MSE_S_r=zeros(1,5);
% MSE_S_v=zeros(1,5);
% MSE_S_C=zeros(bglength*bglength,5);
% %0.1
%deltar=[0,logspace(-3,0,4)];
 % % % % % deltactr=1;
 % % % % %  for delta=[0,logspace(-3,0,4)]
 % % % % %     betactr=1;
 % % % % %     %for beta=linspace(beta_ODE*0.9,beta_ODE*1.1,bglength) %beta_ODE = .0008;     gamma_ODE = 0.08;
 % % % % %     for beta=spacer   
 % % % % %         %beta=spacer(betaindex);
 % % % % %         gammactr=1;
 % % % % %         %for gamma=linspace(gamma_ODE*0.9,gamma_ODE*1.1,bglength)
 % % % % %         for gamma=spacer
             %gamma=spacer(gammaindex);
delta=1; beta=1; gamma=1;
% delta1=1/timidivi; beta1=1/timidivi; gamma1=1/timidivi;
for Rounder_test=1:1:4
    if Rounder_test==1
try
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;
%delta=1; beta=1; gamma=1;

%% Max iterations for Petri Net and ODE, petri net does 1 step per iteration
%in base settings
MAX_ITERATIONS = 100*timidivi;

%%%%% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 


%% Parameter values  %%%%%%%%%%%%%%%%%%%%%
%Used for setting specific parameter values (vs loop above for grid of
%param values)

%beta = 0.001;     gamma = 0.05;    delta = 0.001;
%beta = 0.1;     gamma = 0.2;    delta = 0.1;
%beta = 0.2;     gamma = 0.5;    delta = 0.4;
%beta = 3.1000e-04;     gamma =  0.0809;    
%delta = 0.00;


%% Initial Population
S_0=1000;
I_0=10;
R_0=10;

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

%% the loop   %%%%%%%%%%%%%%%%%%%%%%%%%
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;

while le(Number_of_iterations, MAX_ITERATIONS)
    %disp(['Iteration: ', int2str(Number_of_iterations)]);
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 
    %  Arc weights
    %This below is for splitting time
    % if Number_of_iterations==1
    %     beta=beta1;
    %     delta=delta1;
    %     gamma=gamma1;
    % else
    %     % beta=1;
    %     % delta=1;
    %     % gamma=1;
    % end
    if beta_T==0 || m0Infected==0
        global_info.pSus_tInf = zero;
        pSus_tInf_resEr=(beta_T * m0Susceptible *  m0Infected + pSus_tInf_resEr)-round(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
        global_info.tInf_pInf = zero;
        tInf_pInf_resEr=(beta_T * m0Susceptible *  m0Infected+tInf_pInf_resEr)-round(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
    elseif m0Infected<=(1/(beta_T))
        global_info.pSus_tInf = round(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr);
        pSus_tInf_resEr=(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr)-round(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
        global_info.tInf_pInf = round(beta_T * m0Susceptible *   m0Infected + tInf_pInf_resEr);
        tInf_pInf_resEr=(beta_T * m0Susceptible *   m0Infected +tInf_pInf_resEr)-round(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
    else 

        %beta2=(1/(m0Infected+??????));
        beta2=1/(m0Infected+1);
        global_info.pSus_tInf = round(beta2 * m0Susceptible * m0Infected + pSus_tInf_resEr);
        pSus_tInf_resEr=(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr)-round(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr);
        global_info.tInf_pInf = round(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);
        tInf_pInf_resEr=(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr)-round(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);
    
    end

    global_info.pInf_tRec = round(gamma_T * (m0Infected) + pInf_tRec_resEr);
    pInf_tRec_resEr=(gamma_T * (m0Infected) + pInf_tRec_resEr)-round(gamma_T * (m0Infected) + pInf_tRec_resEr);
    global_info.tRec_pRec = round(gamma_T * (m0Infected) + tRec_pRec_resEr);
    tRec_pRec_resEr=(gamma_T * (m0Infected) + tRec_pRec_resEr)-round(gamma_T * (m0Infected) + tRec_pRec_resEr);
    global_info.pRec_tSus = round(delta_T * (m0Recovered) + pRec_tSus_resEr);
    pRec_tSus_resEr=(delta_T * (m0Recovered) + pRec_tSus_resEr)-round(delta_T * (m0Recovered) + pRec_tSus_resEr);
    global_info.tSus_pSus = round(delta_T * (m0Recovered) + tSus_pSus_resEr);
    tSus_pSus_resEr=(delta_T * (m0Recovered)+ tSus_pSus_resEr)-round(delta_T * (m0Recovered) + tSus_pSus_resEr);
    %{ 
%One way to handle rounding error
    if m0Infected<=(1/beta)

    if round(beta * m0Susceptible * m0Infected)==0 && (beta * m0Susceptible * m0Infected)~=0 && MinfiringS==0 
        timestepstofireS=round((1/(beta * m0Susceptible * m0Infected)));
        tstf_counterS=tstf_counterS+1;
        MinfiringS=1;
        global_info.pSus_tInf = round(beta * m0Susceptible * m0Infected);
        global_info.tInf_pInf = round(beta * m0Susceptible * m0Infected);
        (beta * m0Susceptible * m0Infected)
    elseif MinfiringS==0
        global_info.pSus_tInf = round(beta * m0Susceptible * m0Infected);
        global_info.tInf_pInf = round(beta * m0Susceptible * m0Infected);
        (beta * m0Susceptible * m0Infected)
    elseif MinfiringS==1 && tstf_counterS<timestepstofireS
        tstf_counterS=tstf_counterS+1;
        global_info.pSus_tInf = 0;
        global_info.tInf_pInf = 0;
        (beta * m0Susceptible * m0Infected)
    elseif MinfiringS==1 && tstf_counterS==timestepstofireS
        global_info.pSus_tInf = 1;
        global_info.tInf_pInf = 1;
        (beta * m0Susceptible * m0Infected)
        MinfiringS=0;
        tstf_counterS=0;
    else
        error("logic loop broken")
    end
 
    
    else

        beta2=1/m0Infected;
    if round(beta2 * m0Susceptible * m0Infected)==0 && (beta2 * m0Susceptible * m0Infected)~=0 && MinfiringS==0 
        timestepstofireS=round(1/(beta2 * m0Susceptible * m0Infected));
        tstf_counterS=tstf_counterS+1;
        MinfiringS=1;
        global_info.pSus_tInf = round(beta2 * m0Susceptible * m0Infected);
        global_info.tInf_pInf = round(beta2 * m0Susceptible * m0Infected);
    elseif MinfiringS==0
        global_info.pSus_tInf = round(beta2 * m0Susceptible * m0Infected);
        global_info.tInf_pInf = round(beta2 * m0Susceptible * m0Infected);
    elseif MinfiringS==1 && tstf_counterS<timestepstofireS
        tstf_counterS=tstf_counterS+1;
        global_info.pSus_tInf = 0;
        global_info.tInf_pInf = 0;
    elseif MinfiringS==1 && tstf_counterS==timestepstofireS
        global_info.pSus_tInf = 1;
        global_info.tInf_pInf = 1;
        MinfiringS=0;
        tstf_counterS=0;
    else
        error("logic loop broken")
    end

    end
   

    % global_info.pInf_tInf = ceil(beta * m0Susceptible * m0Infected);
    if round(gamma * m0Infected)==0 && (gamma * m0Infected)~=0 && MinfiringI==0 
        timestepstofireI=round(1/(gamma * m0Infected));
        tstf_counterI=tstf_counterI+1;
        MinfiringI=1;
        global_info.pInf_tRec = round(gamma * m0Infected);
        global_info.tRec_pRec = round(gamma * m0Infected);
    elseif MinfiringI==0
        global_info.pInf_tRec = round(gamma * m0Infected);
        global_info.tRec_pRec = round(gamma * m0Infected);
    elseif MinfiringI==1 && tstf_counterI<timestepstofireI
        tstf_counterI=tstf_counterI+1;
        global_info.pInf_tRec = 0;
        global_info.tRec_pRec = 0;
    elseif MinfiringI==1 && tstf_counterI==timestepstofireI
        global_info.pInf_tRec =1;
        global_info.tRec_pRec = 1;
        MinfiringI=0;
        tstf_counterI=0;
    else
        error("logic loop broken")
    end


    if round(delta * m0Recovered)==0 && (delta * m0Recovered)~=0 && MinfiringR==0 
        timestepstofireR=round(1/(delta * m0Recovered));
        tstf_counterR=tstf_counterR+1;
        MinfiringR=1;
        global_info.pRec_tSus = round(delta * m0Recovered);
        global_info.tSus_pSus = round(delta * m0Recovered);
    elseif MinfiringR==0
        global_info.pRec_tSus = round(delta * m0Recovered);
        global_info.tSus_pSus = round(delta * m0Recovered);
    elseif MinfiringR==1 && tstf_counterR<timestepstofireR
        tstf_counterR=tstf_counterR+1;
        global_info.pRec_tSus = 0;
        global_info.tSus_pSus = 0;
    elseif MinfiringR==1 && tstf_counterI==timestepstofireI
        global_info.pRec_tSus = 1;
        global_info.tSus_pSus = 1;
        MinfiringR=0;
        tstf_counterR=0;
    else
        error("logic loop broken")
    end
    %}


    % if floor(beta * m0Susceptible * m0Infected)>=0
    % global_info.pSus_tInf = floor(beta * m0Susceptible * m0Infected);
    % global_info.tInf_pInf = floor(beta * m0Susceptible * m0Infected);
    % else
    % global_info.pSus_tInf = 0;
    % global_info.tInf_pInf = 0;
    % end
    % if floor(gamma * m0Infected)>=0
    % global_info.pInf_tRec = floor(gamma * m0Infected);
    % global_info.tRec_pRec = floor(gamma * m0Infected);
    % else
    % global_info.pInf_tRec = 0;
    % global_info.tRec_pRec = 0;
    % end
    % if floor(delta * m0Recovered)>=0
    % global_info.pRec_tSus = floor(delta * m0Recovered);
    % global_info.tSus_pSus = floor(delta * m0Recovered);
    % else
    % global_info.pRec_tSus = 0;
    % global_info.tSus_pSus = 0;    
    % end
    

    pns = pnstruct('simple_pn_SIRv1_pdf');
    dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    %figure()
    %plotp(Sim_Results, {'pSusceptible','pInfected','pRecovered'});

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
catch

    disp(['Error at delta=', int2str(delta), ' , gamma=',int2str(gamma), ' , beta=',int2str(beta), ' , timesteps=',int2str(timidivi)]);
    
    if delta==0
    delta=delta; 
    else
    delta=delta-10^-7; 
    end
   
    if beta==0
    beta=beta; 
    else
    beta=beta-10^-7;
    end
    
    if gamma==0
    gamma=gamma;
    else
    gamma=gamma-10^-7;
    end
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

%% Max iterations for Petri Net and ODE, petri net does 1 step per iteration
%in base settings
MAX_ITERATIONS = 100*timidivi;

%%%%% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

%% Initial Population
S_0=1000;
I_0=10;
R_0=10;

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

% the loop   %%%%%%%%%%%%%%%%%%%%%%%%%
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;
while le(Number_of_iterations, MAX_ITERATIONS)
    %disp(['Iteration: ', int2str(Number_of_iterations)]);
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 

    if beta_T==0 || m0Infected==0
        global_info.pSus_tInf = zero;
        pSus_tInf_resEr=(beta_T * m0Susceptible *  m0Infected + pSus_tInf_resEr)-round(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
        global_info.tInf_pInf = zero;
        tInf_pInf_resEr=(beta_T * m0Susceptible *  m0Infected+tInf_pInf_resEr)-round(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
    elseif m0Infected<=(1/(beta_T))
        global_info.pSus_tInf = round(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr);
        pSus_tInf_resEr=(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr)-round(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
        global_info.tInf_pInf = round(beta_T * m0Susceptible *   m0Infected + tInf_pInf_resEr);
        tInf_pInf_resEr=(beta_T * m0Susceptible *   m0Infected +tInf_pInf_resEr)-round(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
    else 

        %beta2=(1/(m0Infected+??????));
        beta2=1/(m0Infected+1);
        global_info.pSus_tInf = round(beta2 * m0Susceptible * m0Infected + pSus_tInf_resEr);
        pSus_tInf_resEr=(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr)-round(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr);
        global_info.tInf_pInf = round(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);
        tInf_pInf_resEr=(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr)-round(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);
    
    end

    global_info.pInf_tRec = round(gamma_T * (m0Infected) + pInf_tRec_resEr);
    pInf_tRec_resEr=(gamma_T * (m0Infected) + pInf_tRec_resEr)-round(gamma_T * (m0Infected) + pInf_tRec_resEr);
    global_info.tRec_pRec = round(gamma_T * (m0Infected) + tRec_pRec_resEr);
    tRec_pRec_resEr=(gamma_T * (m0Infected) + tRec_pRec_resEr)-round(gamma_T * (m0Infected) + tRec_pRec_resEr);
    global_info.pRec_tSus = round(delta_T * (m0Recovered) + pRec_tSus_resEr);
    pRec_tSus_resEr=(delta_T * (m0Recovered) + pRec_tSus_resEr)-round(delta_T * (m0Recovered) + pRec_tSus_resEr);
    global_info.tSus_pSus = round(delta_T * (m0Recovered) + tSus_pSus_resEr);
    tSus_pSus_resEr=(delta_T * (m0Recovered)+ tSus_pSus_resEr)-round(delta_T * (m0Recovered) + tSus_pSus_resEr);
    

    pns = pnstruct('simple_pn_SIRv1_pdf');
    %dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    %figure()
    %plotp(Sim_Results, {'pSusceptible','pInfected','pRecovered'});

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
end


MS_ER=MSusceptible;
MI_ER=MInfected;
MR_ER=MRecovered;


    elseif Rounder_test==2 %Floor
try
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;
%delta=1; beta=1; gamma=1;

%% Max iterations for Petri Net and ODE, petri net does 1 step per iteration
%in base settings
MAX_ITERATIONS = 100*timidivi;

%%%%% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 


%% Parameter values  %%%%%%%%%%%%%%%%%%%%%
%Used for setting specific parameter values (vs loop above for grid of
%param values)

%beta = 0.001;     gamma = 0.05;    delta = 0.001;
%beta = 0.1;     gamma = 0.2;    delta = 0.1;
%beta = 0.2;     gamma = 0.5;    delta = 0.4;
%beta = 3.1000e-04;     gamma =  0.0809;    
%delta = 0.00;


%% Initial Population
S_0=1000;
I_0=10;
R_0=10;

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

%% the loop   %%%%%%%%%%%%%%%%%%%%%%%%%
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;

while le(Number_of_iterations, MAX_ITERATIONS)
    %disp(['Iteration: ', int2str(Number_of_iterations)]);
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 


    if floor(beta * m0Susceptible * m0Infected)>=0
    global_info.pSus_tInf = floor(beta * m0Susceptible * m0Infected);
    global_info.tInf_pInf = floor(beta * m0Susceptible * m0Infected);
    else
    global_info.pSus_tInf = 0;
    global_info.tInf_pInf = 0;
    end
    if floor(gamma * m0Infected)>=0
    global_info.pInf_tRec = floor(gamma * m0Infected);
    global_info.tRec_pRec = floor(gamma * m0Infected);
    else
    global_info.pInf_tRec = 0;
    global_info.tRec_pRec = 0;
    end
    if floor(delta * m0Recovered)>=0
    global_info.pRec_tSus = floor(delta * m0Recovered);
    global_info.tSus_pSus = floor(delta * m0Recovered);
    else
    global_info.pRec_tSus = 0;
    global_info.tSus_pSus = 0;    
    end
    

    pns = pnstruct('simple_pn_SIRv1_pdf');
    dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    %figure()
    %plotp(Sim_Results, {'pSusceptible','pInfected','pRecovered'});

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
catch  % Floor Catch

    disp(['Error at delta=', int2str(delta), ' , gamma=',int2str(gamma), ' , beta=',int2str(beta), ' , timesteps=',int2str(timidivi)]);
    
    if delta==0
    delta=delta; 
    else
    delta=delta-10^-7; 
    end
   
    if beta==0
    beta=beta; 
    else
    beta=beta-10^-7;
    end
    
    if gamma==0
    gamma=gamma;
    else
    gamma=gamma-10^-7;
    end
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

%% Max iterations for Petri Net and ODE, petri net does 1 step per iteration
%in base settings
MAX_ITERATIONS = 100*timidivi;

%%%%% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

%% Initial Population
S_0=1000;
I_0=10;
R_0=10;
m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

%% the loop   %%%%%%%%%%%%%%%%%%%%%%%%%
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;
while le(Number_of_iterations, MAX_ITERATIONS)
    %disp(['Iteration: ', int2str(Number_of_iterations)]);
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 

    if floor(beta * m0Susceptible * m0Infected)>=0
    global_info.pSus_tInf = floor(beta * m0Susceptible * m0Infected);
    global_info.tInf_pInf = floor(beta * m0Susceptible * m0Infected);
    else
    global_info.pSus_tInf = 0;
    global_info.tInf_pInf = 0;
    end
    if floor(gamma * m0Infected)>=0
    global_info.pInf_tRec = floor(gamma * m0Infected);
    global_info.tRec_pRec = floor(gamma * m0Infected);
    else
    global_info.pInf_tRec = 0;
    global_info.tRec_pRec = 0;
    end
    if floor(delta * m0Recovered)>=0
    global_info.pRec_tSus = floor(delta * m0Recovered);
    global_info.tSus_pSus = floor(delta * m0Recovered);
    else
    global_info.pRec_tSus = 0;
    global_info.tSus_pSus = 0;    
    end

    pns = pnstruct('simple_pn_SIRv1_pdf');
    %dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    %figure()
    %plotp(Sim_Results, {'pSusceptible','pInfected','pRecovered'});

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
end
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MSusceptible, 'r--o','LineWidth',2);
MS_Floor=MSusceptible;
MI_Floor=MInfected;
MR_Floor=MRecovered;
    


%%CEILING
elseif Rounder_test==3 %Ceil
try
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;
%delta=1; beta=1; gamma=1;

%% Max iterations for Petri Net and ODE, petri net does 1 step per iteration
%in base settings
MAX_ITERATIONS = 100*timidivi;

%%%%% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 


%% Parameter values  %%%%%%%%%%%%%%%%%%%%%
%Used for setting specific parameter values (vs loop above for grid of
%param values)

%beta = 0.001;     gamma = 0.05;    delta = 0.001;

%% Initial Population
S_0=1000;
I_0=10;
R_0=10;

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

%% the loop   %%%%%%%%%%%%%%%%%%%%%%%%%
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;

while le(Number_of_iterations, MAX_ITERATIONS)
    %disp(['Iteration: ', int2str(Number_of_iterations)]);
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 
   


    global_info.pSus_tInf = ceil(beta * m0Susceptible * m0Infected);
    global_info.tInf_pInf = ceil(beta * m0Susceptible * m0Infected);
   
   
    global_info.pInf_tRec = ceil(gamma * m0Infected);
    global_info.tRec_pRec = ceil(gamma * m0Infected);
   
   
    global_info.pRec_tSus = ceil(delta * m0Recovered);
    global_info.tSus_pSus = ceil(delta * m0Recovered);
    

    pns = pnstruct('simple_pn_SIRv1_pdf');
    dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    %figure()
    %plotp(Sim_Results, {'pSusceptible','pInfected','pRecovered'});

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
catch

    disp(['Error at delta=', int2str(delta), ' , gamma=',int2str(gamma), ' , beta=',int2str(beta), ' , timesteps=',int2str(timidivi)]);
    
    if delta==0
    delta=delta; 
    else
    delta=delta-10^-7; 
    end
   
    if beta==0
    beta=beta; 
    else
    beta=beta-10^-7;
    end
    
    if gamma==0
    gamma=gamma;
    else
    gamma=gamma-10^-7;
    end
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

%% Max iterations for Petri Net and ODE, petri net does 1 step per iteration
%in base settings
MAX_ITERATIONS = 100*timidivi;

%%%%% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

%% Initial Population
S_0=1000;
I_0=10;
R_0=10;

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

%% the loop   %%%%%%%%%%%%%%%%%%%%%%%%%
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;
while le(Number_of_iterations, MAX_ITERATIONS)
    %disp(['Iteration: ', int2str(Number_of_iterations)]);
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 

   global_info.pSus_tInf = ceil(beta * m0Susceptible * m0Infected);
    global_info.tInf_pInf = ceil(beta * m0Susceptible * m0Infected);
   
   
    global_info.pInf_tRec = ceil(gamma * m0Infected);
    global_info.tRec_pRec = ceil(gamma * m0Infected);
   
   
    global_info.pRec_tSus = ceil(delta * m0Recovered);
    global_info.tSus_pSus = ceil(delta * m0Recovered);

    pns = pnstruct('simple_pn_SIRv1_pdf');
    %dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    %figure()
    %plotp(Sim_Results, {'pSusceptible','pInfected','pRecovered'});

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
end
MS_Ceil=MSusceptible;
MI_Ceil=MInfected;
MR_Ceil=MRecovered;




%%standar round
elseif Rounder_test==4 %Round standard
try
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;
%delta=1; beta=1; gamma=1;

%% Max iterations for Petri Net and ODE, petri net does 1 step per iteration
%in base settings
MAX_ITERATIONS = 100*timidivi;

%%%%% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 


%% Parameter values  %%%%%%%%%%%%%%%%%%%%%
%Used for setting specific parameter values (vs loop above for grid of
%param values)

%beta = 0.001;     gamma = 0.05;    delta = 0.001;

%% Initial Population
S_0=1000;
I_0=10;
R_0=10;

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

%% the loop   %%%%%%%%%%%%%%%%%%%%%%%%%
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;

while le(Number_of_iterations, MAX_ITERATIONS)
    %disp(['Iteration: ', int2str(Number_of_iterations)]);
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered;

    global_info.pSus_tInf = round(beta * m0Susceptible * m0Infected);
    global_info.tInf_pInf = round(beta * m0Susceptible * m0Infected);
   
   
    global_info.pInf_tRec = round(gamma * m0Infected);
    global_info.tRec_pRec = ceil(gamma * m0Infected);
   
   
    global_info.pRec_tSus = round(delta * m0Recovered);
    global_info.tSus_pSus = round(delta * m0Recovered);
    

    pns = pnstruct('simple_pn_SIRv1_pdf');
    dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    %figure()
    %plotp(Sim_Results, {'pSusceptible','pInfected','pRecovered'});

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
catch

    disp(['Error at delta=', int2str(delta), ' , gamma=',int2str(gamma), ' , beta=',int2str(beta), ' , timesteps=',int2str(timidivi)]);
    
    if delta==0
    delta=delta; 
    else
    delta=delta-10^-7; 
    end
   
    if beta==0
    beta=beta; 
    else
    beta=beta-10^-7;
    end
    
    if gamma==0
    gamma=gamma;
    else
    gamma=gamma-10^-7;
    end
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

%% Max iterations for Petri Net and ODE, petri net does 1 step per iteration
%in base settings
MAX_ITERATIONS = 100*timidivi;

%%%%% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

%% Initial Population
S_0=1000;
I_0=10;
R_0=10;

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

%% the loop   %%%%%%%%%%%%%%%%%%%%%%%%%
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;
while le(Number_of_iterations, MAX_ITERATIONS)
    %disp(['Iteration: ', int2str(Number_of_iterations)]);
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 

   global_info.pSus_tInf = round(beta * m0Susceptible * m0Infected);
    global_info.tInf_pInf = round(beta * m0Susceptible * m0Infected);
   
   
    global_info.pInf_tRec = round(gamma * m0Infected);
    global_info.tRec_pRec = round(gamma * m0Infected);
   
   
    global_info.pRec_tSus = round(delta * m0Recovered);
    global_info.tSus_pSus = round(delta * m0Recovered);

    pns = pnstruct('simple_pn_SIRv1_pdf');
    %dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    %figure()
    %plotp(Sim_Results, {'pSusceptible','pInfected','pRecovered'});

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
end
MS_SR=MSusceptible;
MI_SR=MInfected;
MR_SR=MRecovered;
    end
end
if ode_run_counter==1
y0=[S_0;I_0;R_0];
delta_ODE=delta; 
beta_ODE=beta; 
gamma_ODE=gamma;
params=[delta_ODE;beta_ODE;gamma_ODE];
[t,y]=ode15s(@SIR,[1:1:(MAX_ITERATIONS/timidivi)],y0,[],params);
else
end

figure()
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MSusceptible, 'b--o','LineWidth',2);
hold on
plot(t,y(:,1), 'y','LineWidth',2)
hold off
%% RRMSE
 err_S = (rmse(MSusceptible(1:timidivi:MAX_ITERATIONS) , y(:,1)')/sqrt(sumsqr(y(:,1)')))*100;
 err_I = (rmse(MInfected(1:timidivi:MAX_ITERATIONS) , y(:,2)')/sqrt(sumsqr( y(:,2)')))*100;
 err_R = (rmse(MRecovered(1:timidivi:MAX_ITERATIONS) , y(:,3)')/sqrt(sumsqr(y(:,3)')))*100;
 % % 
 % % timidivi_col_S=[timidivi_col_S err_S];
 % % timidivi_col_I=[timidivi_col_I err_I];
 % % timidivi_col_R=[timidivi_col_R err_R];
% %end
 %% Plot times steps vs rrmse %%%%%%%%%%%%%%%%%%
% figure()
% plot(timidivi_timestep, timidivi_col_S, 'b--o','LineWidth',2);
% hold on
% plot(timidivi_timestep, timidivi_col_I,'c--o','LineWidth',2);
% plot(timidivi_timestep, timidivi_col_R, 'r--o','LineWidth',2);
% legend('Susceptible_{RRMSE}','Infected_{RRMS}','Recovered_{RRMSE}');
% title("SIR Petri Net vs ODE RRMSE based on times steps", 'FontSize', 24)
% subtitle("beta="+beta*timidivi+', gamma='+gamma*timidivi+' delta='+delta*timidivi)
% xlabel('PN Time steps (per 1 ODE time interval)'); ylabel('RRMSE');
% hold off

figure()
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MS_Ceil, 'color', [.15 .15 .15],'LineWidth',12);
hold on
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MS_Floor,'color', [.65 .65 .65],'LineWidth',8);
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MS_SR, 'color', [.90 .90 .90],'LineWidth',4);
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MS_ER, 'b','LineWidth',4);
plot(t,y(:,1), 'r','LineWidth',2); 
lgd=legend('Ceiling','Floor','Standard','Standard+Residuals','ODE');legend('Location','best');
lgd.Title.String = 'Arc Weight Rounding Method';
title("SIR PN Rounding Methods vs ODE Susceptible", 'FontSize', 24)
subtitle("gamma (\gamma )="+gamma+", beta (\beta)="+beta+", delta (\delta)="+delta)
xlabel('Time')
ylabel('Population') 
xlim([0,105])
hold off

figure()
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MI_Ceil,'color', [.15 .15 .15],'LineWidth',12);
hold on
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MI_Floor,'color', [.65 .65 .65],'LineWidth',8);
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MI_SR, 'color', [.90 .90 .90],'LineWidth',4);
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MI_ER, 'b','LineWidth',4);
plot(t,y(:,2), 'r','LineWidth',2);
lgd=legend('Ceiling','Floor','Standard','Standard+Residuals','ODE');legend('Location','best');
lgd.Title.String = 'Arc Weight Rounding Method';
title("SIR PN Rounding Methods vs ODE Infected", 'FontSize', 24)
subtitle("gamma (\gamma )="+gamma+", beta (\beta)="+beta+", delta (\delta)="+delta)
xlabel('Time')
ylabel('Population') 
xlim([0,105])
hold off

figure()
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MR_Ceil,'color', [.15 .15 .15],'LineWidth',12);
hold on
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MR_Floor,'color', [.65 .65 .65],'LineWidth',8);
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MR_SR, 'color', [.90 .90 .90],'LineWidth',4);
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MR_ER, 'b','LineWidth',4);
plot(t,y(:,3), 'r','LineWidth',2);  
lgd=legend('Ceiling','Floor','Standard','Standard+Residuals','ODE');legend('Location','best');
lgd.Title.String = 'Arc Weight Rounding Method';
title("SIR PN Rounding Methods vs ODE Recovered", 'FontSize', 24)
subtitle("gamma (\gamma )="+gamma+", beta (\beta)="+beta+", delta (\delta)="+delta)
xlabel('Time')
ylabel('Population') 
xlim([0,105])
hold off
% % disp(['still going: ', int2str(deltactr)]);

%% Plot the results of single param value run vs ODE %%%%%%%%%%%%%%%%%%
 figure()
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MSusceptible, 'b--o','LineWidth',2);
hold on
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MInfected,'c--o','LineWidth',2);
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MRecovered, 'r--o','LineWidth',2);

plot(t,y(:,1), 'y','LineWidth',2)
plot(t,y(:,2),'m','LineWidth',2)
plot(t,y(:,3), 'k','LineWidth',2);
legend('Susceptible_{PN}','Infected_{PN}','Recovered_{PN}',...
    'Susceptible_{ODE}','Infected_{ODE}','Recovered_{ODE}');
title("SIR Petri Net vs ODE", 'FontSize', 24)
subtitle("beta="+beta*timidivi+', gamma='+gamma*timidivi+' delta='+delta*timidivi+", Susc. MSE="+err_S...
   +", Infe. MSE="+err_I+", Reco. MSE="+err_R, 'FontSize', 14)
hold off
% 
% figure()
% plot([1,1+1/timidivi,2:1:100], MSusceptible, 'b--o','LineWidth',2);
% hold on
% plot([1,1+1/timidivi,2:1:100], MInfected,'c--o','LineWidth',2);
% plot([1,1+1/timidivi,2:1:100], MRecovered, 'r--o','LineWidth',2);
% 
% plot(t,y(:,1), 'y','LineWidth',2)
% plot(t,y(:,2),'m','LineWidth',2)
% plot(t,y(:,3), 'k','LineWidth',2);
% legend('Susceptible_{PN}','Infected_{PN}','Recovered_{PN}',...
%     'Susceptible_{ODE}','Infected_{ODE}','Recovered_{ODE}');
% title("SIR Petri Net vs ODE", 'FontSize', 24)
% subtitle("beta="+beta*timidivi+', gamma='+gamma*timidivi+' delta='+delta*timidivi+", Susc. MSE="+err_S...
%    +", Infe. MSE="+err_I+", Reco. MSE="+err_R, 'FontSize', 14)
% hold off


        % MSE_s(gammactr,betactr)=err_S;
        % MSE_i(gammactr,betactr)=err_I;
        % MSE_iBeta(gammactr,betactr)=beta;
        % MSE_iGamma(gammactr,betactr)=gamma;
        % MSE_r(gammactr,betactr)=err_R;

        MSE_s(gammactr,betactr,deltactr)=err_S;
        MSE_i(gammactr,betactr,deltactr)=err_I;
        MSE_r(gammactr,betactr,deltactr)=err_R;
        % 
% MSE_S_c(deltactr)=zeros(1,5);
% MSE_S_r(deltactr)=zeros(1,5);
% MSE_S_v(deltactr)=zeros(1,5);
% MSE_S_C(dr)=[MSE_S_C(dr) err_S];
% % 
  % % % % %          gammactr=gammactr+1;
  % % % % %        end
  % % % % %       betactr=betactr+1;
  % % % % %    end
  % % % % %    deltactr=deltactr+1;
  % % % % % end
  
  
%   tau_S=[tau_S (sum(MSE_s,"all")/bglength^3)];
%   tau_I=[tau_I (sum(MSE_i,"all")/bglength^3)];
%   tau_R=[tau_R (sum(MSE_r,"all")/bglength^3)];
% end
% figure()
% plot(tau_spread, tau_S, 'b--o','LineWidth',2);
% hold on
% plot(tau_spread, tau_I, 'c--o','LineWidth',2);
% plot(tau_spread, tau_R, 'r--o','LineWidth',2);
% legend('Susceptible Mean RRMSE','Infected Mean RRMSE','Recovered Mean RRMSE');
% title("Mean Error for different Parameter ranges (0,\tau)", 'FontSize', 24)
% xlabel("\tau (parameters \gamma,\beta,\delta range [0,\tau])");ylabel('RRMSE Mean')
% subtitle("beta="+beta+', gamma='+gamma+' delta='+delta+", Susc. MSE="+err_S...
%     +", Infe. MSE="+err_I+", Reco. MSE="+err_R, 'FontSize', 14)
%  hold off


% minValue = min(MSE_i(:));
% % Find all (row, column) pairs where M = the min value.
% [rows, columns] = find(MSE_i == minValue)
% % Print them out:
% for k = 1 : length(rows)
% 	fprintf('M equals the min value of %f at row = %d, column = %d.\n', ...
% 		minValue, rows(k), columns(k));
% end
% MSE_iBeta(rows(k),columns(k))
%  MSE_iGamma(rows(k),columns(k))

%% Plot the results of grid param value RRMSE vs ODE
% casxismin=min([min(MSE_s(:)),min(MSE_i(:)),min(MSE_r(:))]);
% trial=1;
% if trial==1
% casxismax=max([max(MSE_s(:)),max(MSE_i(:)),max(MSE_r(:))]);
% trial=trial+1;
% else
% end

casxismin=0;
casxismax=44;
Cus_map=[(ones(2,3).*[0 0.5 0]);(ones(8,3).*[0 1 0]);(ones(40,3).*[1 1 0]);(ones(50,3).*[1 0.75 0]);...
    (ones(100,3).*[1 0.4 0]);(ones(240,3).*[1 0 0])];


figure()
%title("RRMSE for parameter Ranges with 1 pn step per 1 ode unit time", 'FontSize', 14)
tiledlayout(5,3);
nexttile
s=surf(spacer,spacer,MSE_s(:,:,1));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Susceptible Delta=0')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h, 'RRMSE (%)');
%h.Label.String = "RRMSE (%)";h.Label.Rotation = 0;
view(2)

nexttile
s=surf(spacer,spacer,MSE_i(:,:,1));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Infected Delta=0')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h, 'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r(:,:,1));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Recovered Delta=0')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h, 'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_s(:,:,2));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Susceptible Delta=0.001')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i(:,:,2));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Infected Delta=0.001')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r(:,:,2));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Recovered Delta=0.001')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_s(:,:,3));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Susceptible Delta=0.01')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i(:,:,3));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Infected Delta=0.01')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r(:,:,3));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Recovered Delta=0.01')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_s(:,:,4));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Susceptible Delta=0.1')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i(:,:,4));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Infected Delta=0.1')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r(:,:,4));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Recovered Delta=0.1')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_s(:,:,5));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Susceptible Delta=1')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i(:,:,5));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Infected Delta=1')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r(:,:,5));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Recovered Delta=1')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
view(2)
h=colorbar('Position',[0.931120331950207,0.11,0.014914910340541,0.816052332195677]);
title(h, 'RRMSE (%)');
%h=axes('Position', [0.85 0.15 0.05 0.7]);% axes('Position', [0.85 0.15 0.05 0.7]),   %'Location', 'westoutside',
%h.Label.String = "RRMSE (%)";h.Label.Rotation = 0;

MSE_s
MSE_i
MSE_r
ode_run_counter=ode_run_counter+1;
%end
end %end of Run Section 5












%% Run_Section 6, loads params from different time steps, showing RRMSE with dif. time steps

if Run_Section==6

load("1stepvariables.mat")
MSE_S_1_M=mean(MSE_s,"all");
MSE_I_1_M=mean(MSE_i,"all");
MSE_R_1_M=mean(MSE_r,"all");

load("20params.mat")
MSE_S_20_M=mean(MSE_s,"all");
MSE_I_20_M=mean(MSE_i,"all");
MSE_R_20_M=mean(MSE_r,"all");

load("40params.mat")
MSE_S_40_M=mean(MSE_s,"all");
MSE_I_40_M=mean(MSE_i,"all");
MSE_R_40_M=mean(MSE_r,"all");

load("60params.mat")
MSE_S_60_M=mean(MSE_s,"all");
MSE_I_60_M=mean(MSE_i,"all");
MSE_R_60_M=mean(MSE_r,"all");

load("80values.mat")
MSE_S_80_M=mean(MSE_s,"all");
MSE_I_80_M=mean(MSE_i,"all");
MSE_R_80_M=mean(MSE_r,"all");

timidivi_col_S=[MSE_S_1_M,MSE_S_20_M,MSE_S_40_M,MSE_S_60_M,MSE_S_80_M];
timidivi_col_I=[MSE_I_1_M,MSE_I_20_M,MSE_I_40_M,MSE_I_60_M,MSE_I_80_M];
timidivi_col_R=[MSE_R_1_M,MSE_R_20_M,MSE_R_40_M,MSE_R_60_M,MSE_R_80_M];
figure()
plot([1,20,40,60,80], timidivi_col_S, 'b--o','LineWidth',2);
hold on
plot([1,20,40,60,80], timidivi_col_I,'c--o','LineWidth',2);
plot([1,20,40,60,80], timidivi_col_R, 'r--o','LineWidth',2);
legend('Susceptible_{RRMSE}','Infected_{RRMS}','Recovered_{RRMSE}');
title("SIR Petri Net vs ODE RRMSE based on times steps", 'FontSize', 24)
%subtitle("beta="+beta*timidivi+', gamma='+gamma*timidivi+' delta='+delta*timidivi)
xlabel('PN Time steps (per 1 ODE time interval)'); ylabel('RRMSE');
hold off

end %End of Run Section 6













%% Run_Section 7, plots comutation time of single run of dif. time steps
if Run_Section==7

%values obtained from finding mean run time of runs both on laptop and SOL
%supercomputer

figure()
plot([1,20,40,60,80], [10.64,160.51,254.77,458.17,611.09], 'r--o','LineWidth',2);
hold on
legend('Computation Time');
title("Computation time for Petri Net depending on times steps", 'FontSize', 24)
subtitle("This is mean run time for 1 PN model")
xlabel('PN Time steps (per 1 ODE time interval)'); ylabel('Mean Computation Time (seconds)');
hold off

end %end of Run section 7

toc


%% ODE Function
function dSIR = SIR(t,a,params)
S = a(1); I = a(2); R = a(3);
delta=params(1); beta=params(2); gamma=params(3);
dS=delta*R-beta*S*I;
dI=beta*S*I-gamma*I;
dR=gamma*I-delta*R;
dSIR=[dS;dI;dR];

end



% function y = optimization_func(x, params)
% 
%     y = params(1) * (params(2) * x);
% end
