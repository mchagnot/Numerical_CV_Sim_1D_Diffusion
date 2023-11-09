%%Numerical Cyclic Voltammetry Simulation from 1D Diffusion 
%%by Matthew Chagnot

% We solve the one-dimensional diffusion equation using the Crank-Nicolson
% method. The film thickness is set by the parameter Lx, and the number of subdivisions is set by
% the parameter Nx. The number of grid points is nx=Nx+1. The spatial grid
% size is dx and the time-step is dt.

%Before running the script, you must define a vector "scanrates" for all of
%the potential scan rates you wish to use in mV/s

%% Looping parameters
%This section is for establishing "lists" to hold measured data from
%multiple loops of the code. e.g. peak currents as a function of scan
%rates, peak potentials, capacities

jlist1 = [];
jlist2 = [];
klist1 = [];
klist2 = [];
llist1 = [];
Eioutmatrix = []; %Matrix which couples each "Eout" with each associated "iout"
iVmatrix = []; %Matrix which couples the full potential datat with with each associated current value
loopcount = 0;

testlist = -0.001*scanrates; %converting mV/s to units of V/s
multiple_scanrates = 1; %if you want to collect b-value information, set multiple_scanrates = 1
%if you are only using one scanrate, set this to 0

%% Main Loop
%this is the workhorse loop of the code. It creates one CV per value of
%"jval" and calls the other functions

for jval = testlist
loopcount = loopcount +1;
Dcm = 5*10^-11; %Diffusion constant in cm^2/s 
D = Dcm*(10^14); %Converting diffusion coefficient to nm
Lx = 100; %film thickness in nm
nsteps = 10000000; %number of time steps (note: this can be as large as you want. the simulation will terminate when the CV reaches its starting potential again and delete all extra timesteps
Eout = 0.2; %plot every Eout potential steps (volts)
Nx = 10*Lx; %number of intervals
nx = Nx+1;%number of gridpoints in x direction including boundaries
dx = Lx/Nx; %grid size in x
x = (0:Nx)*dx; %x values on the grid
dt = 0.1; %time step in seconds
Estart = 0.5 ; %Initial potential value, volts
Eto = -0.5; %Turnover potential, Volts
nu = jval; %Sweep rate, V/s. MUST be a negative number
f = 38.9; %F/RT, V
Eo = 0; %Equilibrium potential, V
R = 0; %Resistance, ohm
F = 96485; %Faraday Constant
Area_cm = 1; %Area in square cm
Area = Area_cm*(10^14); %Area in square nm
CStar = (1.3*(10^-2))*(10^-21); %Concentration of electroactive species, moles per nm^-3
EoutNet = 0; 
Cxmatrix = []; %Matrix of concentration vs. distance, written inside the big for loop so that it resets on every CV

% Next, the nx-by-nx matrix is constructed using spdiags. The matrix is
% altered to represent the condition of (dC/dx) = 0 at the appropriate
% interface
% Here we are setting up our matrix A to use the Crank-Nicolson method

alpha=dt*D/dx^2;
diagonals = [2*(1+alpha)*ones(nx,1), -alpha*ones(nx,2)];
A=spdiags(diagonals,[0 -1 1], nx, nx);
I=speye(nx);
A([1 nx],:)=I([1 nx],:); 
A(Nx,Nx+1) = 0;
A(Nx,Nx) = (1+alpha); 
A(Nx,Nx-1) = -alpha; %original boundaries, modeling diffusion at the film/electrolyte interface

% We next initialize the solution and plot the initial conditions

u_old = CStar*ones(Nx+1,1);
y=u_old./CStar;
figure;
plot(x,y); hold on; %initializing the concentration plot
xlabel('Distance from Electrode/Electrolyte Interface (nm)','Interpreter','latex','FontSize',14);
ylabel('$C_O / C^*(x, t)$','Interpreter','latex','FontSize',14);
title('Concentration Profiles for Selected Time-steps','Interpreter','latex','FontSize',16);
redcounter = 1; % we start counting the number of concentration plots on the reductive and oxidative sweeps for the figure legend
oxcounter = 0;


%this section is basically a grouping of different variables that are
%established before the main for-loop
iguess = 0; %initializing the guess for the first timestep
ioutlist = double.empty(nsteps,0); %this list holds the current values for each Eout
EEffList = double.empty(nsteps,0); 
EAppOutList = double.empty(nsteps,0);
tic
linenum = 1;
EEffoutstring = [string(Estart)];
EEffoutvals = [];
ioutvals = [];
iout = [];
iterationslist = zeros(0,nsteps);
stopcode = 0;
i1list = [];
i2list = [];
i3list = [];
cathodic_sweep = 1; %track if I am in a cathodic sweep or an anodic sweep
mto = 0; %tracking the step where I cross my turnover potential

% This for loop calculates current and potential for every time-step in
% nsteps
for m=1:nsteps
iterations = 0;
if m >= 3
    iguess = (ioutlist(m-1)+(ioutlist(m-1)-ioutlist(m-2)));
end

[EEffEst,iout,u_new,y,EAppOut,iguessout] = crank_nicolson(cathodic_sweep,m,Estart,nu,dt,iguess,R,CStar,f,Eo,alpha,u_old,nx,A,dx,F,Area,D,mto);

u_old = u_new;
iguess = iout;
EEffList(m,1) = EEffEst; %making a list of the effective potentials 
ioutlist(m,1) = iout;%making a list of the output currents
EAppOutList(m,1) = EAppOut;
if round(EAppOut-EoutNet,2) == 0 && cathodic_sweep == 0
    plot(x(2:Nx),y(2:Nx));
    Cxmatrix = [Cxmatrix,x(:),y(:)];
    EEffoutstring = [EEffoutstring,string(round(EEffEst,2))];
    EEffoutvals = [EEffoutvals,EEffEst];
    ioutvals = [ioutvals,iout];
    EoutNet = EoutNet + Eout;
    oxcounter = oxcounter + 1;
elseif round(EAppOut-EoutNet,2) == 0 && cathodic_sweep == 1
    plot(x(2:Nx),y(2:Nx));
    Cxmatrix = [Cxmatrix,x(:),y(:)];  
    EEffoutstring = [EEffoutstring,string(round(EEffEst,2))];
    EEffoutvals = [EEffoutvals,EEffEst];
    ioutvals = [ioutvals,iout];
    EoutNet = EoutNet - Eout;
    redcounter = redcounter + 1;
end
if EAppOut < Eto && cathodic_sweep == 1
    cathodic_sweep = 0; %switch from cathodic to anodic sweep when I go below my turnover potential
    mto = m;
    EoutNet = EoutNet + Eout;
end
%deleting all of the excess timesteps once I return to my starting
%potential
if cathodic_sweep == 0
    if EAppOut > Estart
        stopcode = 1;
        EAppOutList(m:nsteps) = [];
        ioutlist(m:nsteps) = [];
        EEffList(m:nsteps) = [];
    end
end
if stopcode == 1
    writematrix(Cxmatrix,'concdata.xlsx','Sheet',loopcount,'Range','A2')  
    break
end
end

toc
colors = [];
redvals = linspace(0.4,0.9,redcounter);
oxvals = linspace(0.4,0.9,oxcounter);
for n = 1:redcounter
    colors = [colors;0 0.4 redvals(n)];
end
for n = 1:oxcounter
    colors = [colors;oxvals(n) 0.1 0.2];
end
colororder(colors);
lgd = legend(EEffoutstring,'location','eastoutside');
lgd.Title.String = 'Potential (V)';
box on
hold off
jlist1 = [jlist1;max(ioutlist)];
klist1 = [klist1;EAppOutList(find(ioutlist==jlist1(length(jlist1))))];
jlist2 = [jlist2;min(ioutlist)];
klist2 = [klist2;EAppOutList(find(ioutlist==jlist2(length(jlist2))))];

%Capacity calculations
tvals = 0:dt:(dt*(length(ioutlist)-1)); %generate a vector of the time values in seconds
cumulative_capacity = cumtrapz(tvals,ioutlist); %generate a vector of the cumulative capacity to a given step
cathodic_capacity_Coulombs = min(cumulative_capacity); %cathodic capacity is the minimum of the cumulative values
llist1 = [llist1;cathodic_capacity_Coulombs];
%%

legvals = []; %initializing the vector to hold the legend values
figure %making the voltammogram for applied (un-perturbed) potential
hold on
plot(EAppOutList,ioutlist);
iVmatrix = [EAppOutList(:),ioutlist(:)];
writematrix(iVmatrix,'cvdata.xlsx','Sheet',loopcount,'Range','A2')
scatter(EEffoutvals,ioutvals,'O');
%writing this data to a csv
%nustr = strcat('\nu =',string(abs(nu*1000)),' mV/s');
Eioutmatrix = [Eioutmatrix,EEffoutvals(:),ioutvals(:)];
xline(0);
legvals = [legvals,'E = 0 V'];
%plot formatting
axis([1.2*min(EEffList),1.2*max(EEffList),1.1*min(ioutlist),1.1*max(ioutlist)])
xlabel('E (V)')
ylabel('i (A)')
%Adding variable descriptors to the plot
Dstr = strcat('D_{Li} = ',string(Dcm),' cm^2/s');
nustr = strcat('\nu =',string(abs(nu*1000)),' mV/s');
rstr = strcat('R_u = ',string(R),' \Omega');
vars = {nustr,rstr,Dstr};
dim = [.65 0 .3 .3];
annotation('textbox',dim,'String',vars,'FitBoxToText','on','Interpreter','tex');
title('Cyclic Voltammogram with E_{Applied}','Interpreter','tex','FontSize',16);
box on
hold off
writematrix(Eioutmatrix,'Eiout.xlsx','Range','A2')
% close all
end

%% post-processing

%b-value analysis
%requires multiple scanrates
if multiple_scanrates == 1
    logrates = [ones(length(testlist),1),transpose(log10(-1*testlist))];
    logpeakspos = log10(jlist1);
    logpeaksneg = log10(-1*jlist2);
    cathodic_bval_linreg = logrates\logpeaksneg;
    anodic_bval_linreg = logrates\logpeakspos;
    cathodic_bval = cathodic_bval_linreg(2)
    anodic_bval = anodic_bval_linreg(2)
end

%% functions
%The function that takes the initial conditions and uses it to solve the 1D
%diffusion equation. This function outputs our "effective" potential,
%output current, new concentration gradient, and the difference between the
%input and output currents
function [EEffEst,iout,u_new,y,EAppOut,iguessout] = crank_nicolson(cathodic_sweep,m,Estart,nu,dt,iguess,R,CStar,f,Eo,alpha,u_old,nx,A,dx,F,Area,D,mto)
      if cathodic_sweep == 1
        EEffEst = Estart + nu*(dt*m) + iguess*R; %Initializing potential "guess"
        EAppOut = Estart + nu*(dt*m);
      else
       EEffEst = Estart + nu*(dt*(mto)) -  nu*(dt*(m-mto))  + iguess*R; %Initializing potential "guess"
       EAppOut = Estart + nu*(dt*(mto)) -  nu*(dt*(m-mto));
      end
        Co = CStar*exp(f*(EEffEst - Eo))/(1+exp(f*(EEffEst - Eo))); %calculating concentration at the interface
       b=[Co; [alpha*u_old(1:nx-3) + 2*(1-alpha)*u_old(2:nx-2) + alpha*u_old(3:nx-1)]; [alpha*u_old(nx-2) + 1*(1-alpha)*u_old(nx-1) + 0*alpha*u_old(nx)]; 1];
        u_new=A\b;
        y=u_new./CStar;
        CoOut1 = u_new(1,1); %taking the first two data points inside the film to approximate dC/dx at the film/electrolyte interface
        CoOut2 = u_new(2,1);
        SLSlope = (CoOut2-CoOut1)/dx;
        iout = -1*F*Area*D*SLSlope;
        iguessout = iguess;
end
