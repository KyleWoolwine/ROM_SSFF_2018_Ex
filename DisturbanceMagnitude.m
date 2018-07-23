function [ext_mag] = DisturbanceMagnitude(Dir)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Author: Kyle Woolwine
%       Date Created: 1-17-17
% Date Last Modified: 1-17-17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function Description:
%
% This function models the dynamics of a supersonic flow field using the 4
% steady state solutions gathered in response to a freestream disturbance
% as detailed in 'Reduced Order Modeling of a Supersonic Flow Field' by
% Kyle Woolwine et. al.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Inputs:
%                      Dir = Directory which holds the steady solutions
%                            (ext) and the steamline solution (FlowVarsSL)
%                            at a single x location.
%
% Outputs:
%
%                  ext_mag = Temporary steady state solution associated 
%                            with each disturbance wave (fast acoustic, 
%                            entropy, slow acoustic) added to the rows 
%                            below the initial solution.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Initialization
SaveFlag = 0;                                         % Save Flag
gam = 1.4;                                            % Ratio of specific heats
gm1 = gam -1;
R = 287;                                              % Gas constant of air
load([Dir 'StreamlineSolution.mat'])                  % Load streamline solution
load([Dir 'ext.mat'])                                 % Load steady state solutions
%%
% Get propagation delay along streamlines
[r,~] = size(ySL);                                    % Number of streamlines
[~,c] = size(ext.x);                                  % Number of points in x
ext.y = zeros(r,c);

% Loop through streamlines
for j=1:r
    theSL = FlowVarsSL(j,:,4);                        % Flow angles along current streamline
    MxSL = FlowVarsSL(j,:,1).*cos(theSL);             % x-component of Mach number
    TSL = FlowVarsSL(j,:,3);                          % Temperature (K) along streamline
    ext.y(j,:) = interp1(xSL,ySL(j,:),ext.x,'linear',ySL(j,end));         % Y components of streamlines interpolated to 1D x locations
    % Loop through points in x to find the time delay
    for i=1:c
        % Streamline information

        if c==1
            PtsX = 65;                                % Needs to be in the form (2^n)+1, where n>=3
            L = ext.x(i);                             % End of interpolation points
            X = linspace(xSL(1,1),L,PtsX);            % Flow values interpolated along streamline from shock wave to exit
        else
            if i<c/2
                PtsX = 17;                                % Needs to be in the form (2^n)+1, where n>=3
                L = ext.x(i);                             % End of interpolation points
                X = linspace(0,L,PtsX);                   % Flow values interpolated along streamline from shock wave to exit
            else
                PtsX = 33;                                % Needs to be in the form (2^n)+1, where n>=3
                L = ext.x(i);                             % End of interpolation points
                X = linspace(0,L,PtsX);                   % Flow values interpolated along streamline from shock wave to exit
            end
        end
        Mx = interp1(xSL,MxSL,X,'linear',MxSL(end));       % Interpolate x-component of Mach number
        the = interp1(xSL,theSL,X,'linear',theSL(end));    % Interpolate flow angle
        T = interp1(xSL,TSL,X,'linear',TSL(end));          % Interpolate temperature (K) along streamline
        a = (gam*R*T).^0.5;                                % Speed of sound along current streamline
        dTau_Plus = 1./(a.*(Mx+cos(the)));
        dTau_Ent = 1./(a.*Mx);
        dTau_Minus = 1./(a.*(Mx-cos(the)));

        Tau(1,1,1) = NumericalIntegrator(dTau_Plus,X);  % Time delay of fast acoustic wave from shock to end of streamline
        Tau(1,1,2) = NumericalIntegrator(dTau_Ent,X);   % Time delay of entropy wave from shock to end of streamline
        Tau(1,1,3) = NumericalIntegrator(dTau_Minus,X); % Time delay of slow acoustic wave from shock to end of streamline
        ext.t(j,i,:) = abs(Tau);                             % Time delay for each disturbance wave at current x,y location
    end
end
ext_mag=ext;
% For a given disturbance magnitude, this information can be used to 
% determine any frequency disturbance so its useful to save and avoid 
% costly reruns.
if SaveFlag
    save([Dir 'DisturbanceMagnitude.mat'],'ext_mag');
end
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