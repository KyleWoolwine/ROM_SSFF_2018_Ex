function [ext_dyn] = DisturbanceTransientResponse(ext_mag,dis)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Author: Kyle Woolwine
%       Date Created: 1-20-17
% Date Last Modified: 1-20-17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function Description:
%
% This function models the dynamics of the supersonic flow field in 
% response to freestream perturbations. The function 
% DisturbanceMagnitude.m is used to determeine the 3 
% new tempoary steady state solutions caused by the applied freestream 
% disturbance which moves through the flow field as three disturbance 
% waves. That function also determined the propagation speeds as a function
% x,y location. This function uses that solution to model the disturbance 
% as a sine wave or a step input in  1-D model, depending on the input 
% parameters.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Inputs:
%
%                  ext_mag = Dynamic external flow field. Similar to ext in
%                            structure except it has the temporary steady 
%                            state solution associated with each 
%                            disturbance wave (fast acoustic, entropy, slow
%                            acoustic) added to the rows below the inital
%                            solution.
%
%                 dis.type = Flag to determine type of disturbance. If = 1,
%                            apply step disturbance, else apply sinusoidal
%                            disturbance.
%
%                 dis.freq = Frequency of sinusoidal disturbance (Hz).
%
%                   dis.TF = Total time (s).
%
% Outputs:
%
%                  ext_dyn = External flow field solution and 
%                            geometry. For each flow variable the columns 
%                            correspond to each x-location and the rows 
%                            correspond to time step.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
[~,c] = size(ext_mag.x);                                       % Number of x locations
[r,~] = size(ext_mag.y);                                       % Number of streamlines
ext.x = ext_mag.x;                                             % Pass along x coords of the external geometry
ext.A = ext_mag.A;                                             % Pass along area of the external geometry

dt = 1e-5;                                                     % Time step (s)
t_out = (0:dt:dis.TF)';                                        % Output time vector of disturbance
ext.t = t_out; 
[Nt,~] = size(t_out);                                          % Total number of time steps
ext.U1(1:Nt,:) = ones(Nt,1)*ext_mag.U1(1,:);
ext.U2(1:Nt,:) = ones(Nt,1)*ext_mag.U2(1,:);
ext.U3(1:Nt,:) = ones(Nt,1)*ext_mag.U3(1,:);
PtsY = 65;                                                     % Number of interpolation points, needs to be in the form (2^n)+1, where n>=3

%%%External Inlet
LY = 0.80767;
Xcowl = ext_mag.x(end);                                        % X-coordinate of cowl lip (m)
Ycowl = LY;                                                    % Y-coordinate of cowl lip (m)
Ycb_cl = interp1([0.83069 1.08385], [0.26852 0.38672], Xcowl); % Y location at cowl lip on centerbody
Ystart = Ycb_cl;
A = pi*(Ycowl^2-Ycb_cl^2);                                     % Area at cowl lip, going to need to do this at all x locations eventually
CrossType = 1;                                                 % Type of cross section = axisymmetric
Y = linspace(Ystart,LY,PtsY);                                    % Interpolation points;

for i=1:Nt % Loop through time steps
    for j=1:c % Loop through points in x
        dU1 = zeros(r,1);
        dU2 = zeros(r,1);
        dU3 = zeros(r,1);
        if dis.type==1
            % Determine which disturbance waves have arrived per streamline at each x caused by a freestream step disturbance
            U1_Interp = [ext_mag.U1(1,j), ext_mag.U1(1,j), ext_mag.U1(2,j), ext_mag.U1(2,j), ext_mag.U1(3,j), ext_mag.U1(3,j), ext_mag.U1(4,j), ext_mag.U1(4,j)];
            U2_Interp = [ext_mag.U2(1,j), ext_mag.U2(1,j), ext_mag.U2(2,j), ext_mag.U2(2,j), ext_mag.U2(3,j), ext_mag.U2(3,j), ext_mag.U2(4,j), ext_mag.U2(4,j)];
            U3_Interp = [ext_mag.U3(1,j), ext_mag.U3(1,j), ext_mag.U3(2,j), ext_mag.U3(2,j), ext_mag.U3(3,j), ext_mag.U3(3,j), ext_mag.U3(4,j), ext_mag.U3(4,j)];
            
            if t_out(i) > max(max(max(ext_mag.t))) % No need to keep interpolating if last disturbance time delay is past
                ext.U1(i:Nt,:) = ext_mag.U1(4,:);
                ext.U2(i:Nt,:) = ext_mag.U2(4,:);
                ext.U3(i:Nt,:) = ext_mag.U3(4,:);
                ext_dyn = ext;
                return
            end
        end
        for k=1:r % Loop through streamlines
            t_Interp = [0, ext_mag.t(k,j,1), 1.00001*ext_mag.t(k,j,1), ext_mag.t(k,j,2), 1.00001*ext_mag.t(k,j,2), ext_mag.t(k,j,3), 1.00001*ext_mag.t(k,j,3), dis.TF];
            if dis.type==1 % Step disturbance
                dU1(k,1) = interp1(t_Interp,U1_Interp,t_out(i)); %Interpolate U1 disturbance over new time vector
                dU2(k,1) = interp1(t_Interp,U2_Interp,t_out(i)); %Interpolate U2 disturbance over new time vector
                dU3(k,1) = interp1(t_Interp,U3_Interp,t_out(i)); %Interpolate U3 disturbance over new time vector
                    
            else % Sinusoidal disturbance
                dU1(k,1) = ext_mag.U1(1,j);
                dU2(k,1) = ext_mag.U2(1,j);
                dU3(k,1) = ext_mag.U3(1,j);
                
                if t_out(i) >= ext_mag.t(k,j,1)
                    dU1(k,1) = dU1(k,1) + (ext_mag.U1(2,j)-ext_mag.U1(1,j))*sin(2*pi*dis.freq*(t_out(i)-ext_mag.t(k,j,1))); %Effect from fast acoustic wave caused by sinisoidal disturbance
                    dU2(k,1) = dU2(k,1) + (ext_mag.U2(2,j)-ext_mag.U2(1,j))*sin(2*pi*dis.freq*(t_out(i)-ext_mag.t(k,j,1))); %Effect from fast acoustic wave caused by sinisoidal disturbance
                    dU3(k,1) = dU3(k,1) + (ext_mag.U3(2,j)-ext_mag.U3(1,j))*sin(2*pi*dis.freq*(t_out(i)-ext_mag.t(k,j,1))); %Effect from fast acoustic wave caused by sinisoidal disturbance
                    
                    if t_out(i) >= ext_mag.t(k,j,2)
                        dU1(k,1) = dU1(k,1) + (ext_mag.U1(3,j)-ext_mag.U1(2,j))*sin(2*pi*dis.freq*(t_out(i)-ext_mag.t(k,j,2))); %Effect from fast acoustic and entropy waves caused by sinisoidal disturbance
                        dU2(k,1) = dU2(k,1) + (ext_mag.U2(3,j)-ext_mag.U2(2,j))*sin(2*pi*dis.freq*(t_out(i)-ext_mag.t(k,j,2))); %Effect from fast acoustic and entropy waves caused by sinisoidal disturbance
                        dU3(k,1) = dU3(k,1) + (ext_mag.U3(3,j)-ext_mag.U3(2,j))*sin(2*pi*dis.freq*(t_out(i)-ext_mag.t(k,j,2))); %Effect from fast acoustic and entropy waves caused by sinisoidal disturbance
                    
                        if t_out(i) >= ext_mag.t(k,j,3)
                            dU1(k,1) = dU1(k,1) + (ext_mag.U1(4,j)-ext_mag.U1(3,j))*sin(2*pi*dis.freq*(t_out(i)-ext_mag.t(k,j,3))); %Effect from fast acoustic, entropy and slow acoustic waves caused by sinisoidal disturbance
                            dU2(k,1) = dU2(k,1) + (ext_mag.U2(4,j)-ext_mag.U2(3,j))*sin(2*pi*dis.freq*(t_out(i)-ext_mag.t(k,j,3))); %Effect from fast acoustic, entropy and slow acoustic waves caused by sinisoidal disturbance
                            dU3(k,1) = dU3(k,1) + (ext_mag.U3(4,j)-ext_mag.U3(3,j))*sin(2*pi*dis.freq*(t_out(i)-ext_mag.t(k,j,3))); %Effect from fast acoustic, entropy and slow acoustic waves caused by sinisoidal disturbance
                    
                        end
                    end
                end
            end
        end

        dU1=dU1';
        dU2=dU2';
        dU3=dU3';
        if CrossType == 1 % Axisymmetric cross section
            dU1 = interp1(ext_mag.y(:,j),dU1,Y,'linear','extrap'); % Interpolate rho*A over new points Y
            U1 = NumericalIntegrator(dU1.*Y,Y);                    % Find 1D equivalent
            ext.U1(i,j) = 2*pi*U1/A;                               % Find 1D equivalent
            dU2 = interp1(ext_mag.y(:,j),dU2,Y,'linear','extrap'); % Interpolate rho*A*V over new points Y
            U2 = NumericalIntegrator(dU2.*Y,Y);                    % Find 1D equivalent
            ext.U2(i,j) = 2*pi*U2/A;                               % Find 1D equivalent
            dU3 = interp1(ext_mag.y(:,j),dU3,Y,'linear','extrap'); % Interpolate rho*A*e_tot over new points Y
            U3 = NumericalIntegrator(dU3.*Y,Y);                    % Find 1D equivalent
            ext.U3(i,j) = 2*pi*U3/A;                               % Find 1D equivalent
        else % Rectangular cross section
            dU1 = interp1(ext_mag.y(:,j),dU1,Y,'linear','extrap'); % Interpolate rho*A over new points Y
            U1 = NumericalIntegrator(dU1,Y);                       % Find 1D equivalent
            ext.U1(i,j) = z*U1/A;                                  % Find 1D equivalent
            dU2 = interp1(ext_mag.y(:,j),dU2,Y,'linear','extrap'); % Interpolate rho*A*V over new points Y
            U2 = NumericalIntegrator(dU2,Y);                       % Find 1D equivalent
            ext.U2(i,j) = z*U2/A;                                  % Find 1D equivalent
            dU3 = interp1(ext_mag.y(:,j),dU3,Y,'linear','extrap'); % Interpolate rho*A*e_tot over new points Y
            U3 = NumericalIntegrator(dU3,Y);                       % Find 1D equivalent
            ext.U3(i,j) = z*U3/A;                                  % Find 1D equivalent
        end
    end

end
ext_dyn = ext;
end
function [I_RT] = NumericalIntegrator(f,x)
%%
% This function numerically integrates a function using the Trapezoidal Rule
% combined with Richardson Extrapolation.
% Inputs: 
% fn(x) = function values at discrete points = 2^n+1 and  (n is an integer and > 5)
% x = evaluation points
% a,b = integral bounds
% Outputs:
% I = Integrand
%%
a=x(1,1);
b=x(end);
[~,N]=size(x);
h=(b-a)/(N-1); %Grid size
h2=h*2;
h4=h*4;
Ih=0;
Ih2=0;
Ih4=0;
%[I_Bar]=MyPhatTrap(N,fn,a,b)
%Trapezoidal rule (h)
for i=1:N-1
    Ih=Ih+((f(i)+f(i+1))/2)*h; %Integrand
end
%Trapezoidal rule (2h)
for i=1:2:N-1
    Ih2=Ih2+((f(i)+f(i+2))/2)*h2; %Integrand
end
%Trapezoidal rule (4h)
for i=1:4:N-1
    Ih4=Ih4+((f(i)+f(i+4))/2)*h4; %Integrand
end
[I_RT]=RichTrap(3,[Ih4,Ih2,Ih]); %Richardson extrapolation

end
function [I_Bar]=RichTrap(n,I)
%%
k=0;
l=0;
b=zeros(n,1);
b(1,1)=1;
for i=1:n
    for j=1:n
        A(i,j)=1/(2^(l*k*2));
        k=k+1;
    end
    k=0;
    l=l+1;
end
x=A\b;
I_Bar=0;
for i=1:n
    I_Bar=I_Bar+I(i)*x(i,1);
end
end