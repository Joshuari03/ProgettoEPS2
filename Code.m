clear
close all
clc
%% copyright AV
f = 50; % nominal frequency in Hz
w = 2*pi*f; % omega in [rad/s]
wpu = 1;
t = 0:0.001:5; % time
delta = 0:pi/1536:pi; % delta initi.

% data of generator
deltar0 = 30.0373; % input the power angle
deltar0 = deg2rad(deltar0); % in radians
deltamax = pi-deltar0; % max angle
H = 3; % inertia [s]
Pm = 80; % P generated [MW]
Sb = 90; % S power [MVA]
Pmpu = Pm/Sb; % in p.u.
Pmax = Pmpu/sin(deltar0); % maximum power
%% Calculating delta2 with the equal area criterion

% % tc = 0.12
disp('Stable case')
tc = 0.12; % Clearing time
delta1r = ((w*Pmpu*tc^2)/(4*H))+deltar0; % delta at t_clearing
Pele_tc = Pmax*sin(delta1r); % electrical power at t_clearing
EA1 = Pmpu*(delta1r-deltar0); % energy A1

funA2 = @(x) Pmax*(cos(delta1r)-cos(x))-Pmpu*(x-delta1r)-EA1; % my function
x0 = delta1r+tc*delta1r; % i.c.

delta2r = fsolve(funA2,x0) % solve
Pele_ws = Pmax*sin(delta2r); % electrical power at ws
EA2_check = Pmax*(cos(delta1r)-cos(delta2r))-Pmpu*(delta2r-delta1r); % energy A2 (check)

% % tc = 0.2
disp('Unstable case')
tc2 = 0.2; % Clearing time
delta1r2 = (w*Pmpu*(tc2^2))/(4*H)+deltar0; % delta
Pele_tc_2 = Pmax*sin(delta1r2); % electrical power at t_clearing
EA12 = Pmpu*(delta1r2-deltar0); % energy A1_2

funA22 = @(x) Pmax*(cos(delta1r2)-cos(x))-Pmpu*(x-delta1r2)-EA12; % my function
x02 = delta1r2+tc2*delta1r2; % i.c.

delta2rtc2 = fsolve(funA22,x02) % solve --> it stops without finding a solutions,
% so this number is meaningless

% critical clearing time
deltacr = acos((pi-2*deltar0)*sin(deltar0)-cos(deltar0)) % critical angle
tcr = sqrt(4*H*(deltacr-deltar0)/w/Pmpu) % critical clearing time
Pele_tcr = Pmax*sin(deltacr); % electrical power at tcr
% check
EA1_cr = Pmpu*(deltacr-deltar0)
EA2_cr = Pmax*(cos(deltacr)-cos(deltamax))-Pmpu*(deltamax-deltacr)
%% Plotting the graphs
% % Pe(delta) % initial conditions
figure
func_pe = Pmax*sin(delta); % power function
plot(delta,func_pe,'LineWidth',2);
ylabel('Pe [p.u.]')
xlabel('\delta [rad]')
line([0,3.5],[Pmpu,Pmpu],'Color','red','LineStyle','--') % mech. power
text(0.01,0.95,num2str(Pmpu),'Color','red')
line([0,3.5],[Pmax,Pmax],'Color','red','LineStyle','--') % max power
text(0.01,1.73,num2str(Pmax),'Color','red')
line([deltar0,deltar0],[0,1.8],'Color','red','LineStyle','--') % delta0
text(0.1,0.1,num2str(deltar0),'Color','red')
line([deltamax,deltamax],[0,1.8],'Color','red','LineStyle','--') % delta max
text(2.2,0.1,num2str(deltamax),'Color','red')

% % Pe(delta) % stable
figure
plot(delta,func_pe,'LineWidth',2);
ylabel('Pe [p.u.]')
xlabel('\delta [rad]')
line([0,3.5],[Pmpu,Pmpu],'Color','red','LineStyle','--') % mech. power
text(0.01,0.95,num2str(Pmpu),'Color','red')
line([0,3.5],[Pele_tc,Pele_tc],'Color','red','LineStyle','--')  % power at tc
text(0.01,1.39,num2str(Pele_tc),'Color','red')
line([0,3.5],[Pele_ws,Pele_ws],'Color','red','LineStyle','--') % power at ws
text(0.01,1.76,num2str(Pele_ws),'Color','red')
line([deltar0,deltar0],[0,1.8],'Color','red','LineStyle','--') % delta0
text(0.1,0.1,num2str(deltar0),'Color','red')
line([deltamax,deltamax],[0,1.8],'Color','red','LineStyle','--') % delta max
text(2.2,0.1,num2str(deltamax),'Color','red')
line([delta1r,delta1r],[0,1.8],'Color','red','LineStyle','--') % delta at tc
text(0.9,0.1,num2str(delta1r),'Color','red')
line([delta2r,delta2r],[0,1.8],'Color','red','LineStyle','--') % delta at ws
text(1.35,0.1,num2str(delta2r),'Color','red')
% % filling with color
% hold on
% x = [deltar0 delta1r delta1r deltar0];
% y = [0 0 Pmpu Pmpu];
% fill(x,y,'r')
% x = [delta1r delta2r delta2r 1.076 delta1r];
% y = [Pmpu Pmpu Pele_ws 1.563 Pele_tc];
% fill(x,y,'r')

% % Pe(delta) % unstable
figure
plot(delta,func_pe,'LineWidth',2);
ylabel('Pe [p.u.]')
xlabel('\delta [rad]')
line([0,3.5],[Pmpu,Pmpu],'Color','red','LineStyle','--') % mech. power
text(0.01,0.95,num2str(Pmpu),'Color','red')
line([0,3.5],[Pele_tc_2,Pele_tc_2],'Color','red','LineStyle','--') % power at tc
text(0.01,1.70,num2str(Pele_tc_2),'Color','red')
line([deltar0,deltar0],[0,1.8],'Color','red','LineStyle','--') % delta0
text(0.1,0.1,num2str(deltar0),'Color','red')
line([deltamax,deltamax],[0,1.8],'Color','red','LineStyle','--') % delta max
text(2.2,0.1,num2str(deltamax),'Color','red')
line([delta1r2,delta1r2],[0,1.8],'Color','red','LineStyle','--') % delta at tc
text(1.55,0.1,num2str(delta1r2),'Color','red')
% % filling the space
% hold on
% x = [deltar0 delta1r2 delta1r2 deltar0];
% y = [0 0 Pmpu Pmpu];
% fill(x,y,'r')

% % Pe(delta) % critical angle
figure
plot(delta,func_pe,'LineWidth',2);
ylabel('Pe [p.u.]')
xlabel('\delta [rad]')
line([0,3.5],[Pmpu,Pmpu],'Color','red','LineStyle','--') % mech. power
text(0.01,0.95,num2str(Pmpu),'Color','red')
line([0,3.5],[Pele_tcr,Pele_tcr],'Color','red','LineStyle','--') % power at tcr
text(0.01,1.70,num2str(Pele_tcr),'Color','red')
line([deltar0,deltar0],[0,1.8],'Color','red','LineStyle','--') % delta0
text(0.1,0.1,num2str(deltar0),'Color','red')
line([deltamax,deltamax],[0,1.8],'Color','red','LineStyle','--') % delta max
text(2.2,0.1,num2str(deltamax),'Color','red')
line([deltacr,deltacr],[0,1.8],'Color','red','LineStyle','--') % delta cr.
text(0.9,0.1,num2str(deltacr),'Color','red')

% % delta(t)
figure
func_delta = (w*Pmpu*(t.^2))/(4*H)+deltar0;
plot(t,func_delta,'LineWidth',2);
ylim([0 3.5])
ylabel('\delta [rad]')
xlabel('t [s]')
line([0,0.4],[deltamax,deltamax],'Color','red','LineStyle','--')
text(0.01,deltamax+0.08,'2.61 - max angle','Color','red')

disp('')