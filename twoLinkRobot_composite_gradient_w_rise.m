function twoLinkRobotAdaptive
    close all
    %% Select gains for controller
    K        = 50; %Enter a number
    a1       = 2; %Enter a number
    a2       = 2;
    gamma    = [2.15,1.9,1.2,4.1,1]';
    B        = 3; 
    B2       = 3;
    k2       = 50;
    gains = [K;a1;a2;gamma;B;B2;k2]
    %% Stacked parameter vector
    %Set up parameters for sim
    p1       = 3.473;
    p2       = 0.196;
    p3       = 0.242;
    f1       = 5.3;
    f2       = 1.1;
    theta    = [p1;p2;p3;f1;f2];

    % Simulation final time
    tf   = 30;

    % Initial condition vector (X0 must be same size and "form" as X and XDot below)
    % (i.e., in this sim, X0 = [e0;r0;thetahat0])
    X0   = [ones(23,1);0;0];

    % Options for integration function
    opts = odeset('RelTol',1e-3,'AbsTol',1e-3);

    % Integrate (you can send the paramters theta to the dynamics as seen below)
    [t,STATES] = ode45(@(t,X) twoLinkdynamics(t,X,theta,gains),[0 tf],X0,opts);
    u = zeros(2,length(t));
        for i=1:length(t)
        u(:,i) =  getu(t(i),STATES(i,:),theta,gains);
        end
    % Set up desired trajectory data for plots (enter desired trajectory for your simulation)
    qd = [cos(0.5*t) 2*cos(t)]';

    % Parse integrated states (STATES is the same "form" as X0)
    % (i.e., in this sim, STATES = [e r thetahat] over all time);
    e  = STATES(:,1:2)';
    thetaHat = STATES(:,9:13)';

    % Compute x from e and xd for plotting purposes
    q  = -e + qd;
    thetatilda = repmat(theta,1,length(t))-thetaHat;
    
    figure(1)
    plot(t,qd,'-','LineWidth',2)
    hold on
    ax = gca;
    ax.ColorOrderIndex = 1;
    plot(t,q,':','LineWidth',2)
    title('Actual Vs Desred Trajectories')
    xlabel('time')
    ylabel('Trajectories')
    legend('y = q1','y = q2','y=qd1','y=qd2')
    hold off
   
    % Plot the filtered tracking error
     figure(2)
    plot(t,e,'--','LineWidth',2)
    title(' Tracking Error')
    xlabel('time')
    ylabel('Errors')
    legend('y = e1','y = e2')
    
    % Plot the link torques
    figure(3)
    plot(t,u(1,:),'b',t,u(2,:),'LineWidth',2)
    title('Control Input')
    xlabel('time')
    ylabel('Torques')
    legend('y = u1','y = u2')
    
    % Plot the adaptive estimates vs actual parameters
     figure(4)
    plot(t,repmat(theta,1,length(t)),'-','LineWidth',2)
    hold on
    ax = gca;
    ax.ColorOrderIndex = 1;
    plot(t,thetaHat,':','LineWidth',2)
    title('Adaptive estimate')
    xlabel('time') 
    ylabel('Adaptive estimate') 
    legend('y = p1','y = p2','y=p3','y=fd1','y = fd2','y = p1 hat','y = p2 hat','y=p3 hat','y=fd1 hat','y = fd2 hat')
    hold off
    
     %Plot adaptive error estimates 
    figure(5)
    plot(t,thetatilda,'--','LineWidth',2)
    title('Parameter estimates error ')
    xlabel('time')
    ylabel('Parameter estimates error')
    legend('y = p1tilde','y = p2tilde','y=p3tilde','y=fd1tilde','y = fd2tilde')
    
    
end
function [u] = getu(t,X,theta,gains)
global e2init
%% Parse for controller gains = [K;a1;a2;gamma;B;B2;k2]
K        = gains(1); %Enter a number
a1       = gains(2); %Enter a number
a2       = gains(3);
gamma    = gains(4:8);
gamma    = diag(gamma);
B        = gains(9);
B2       = gains(10);
k2       = gains(11);
%% Parse parameter vector
p1 = theta(1);
p2 = theta(2);
p3 = theta(3);
f1 = theta(4);
f2 = theta(5);


% Desired trajectory and needed derivatives
qd       = [cos(0.5*t);2*cos(t)];
qdDot    = [-0.5*sin(0.5*t);-2*sin(t)]; %Enter the expression
qdDotDot = [-0.25*cos(0.5*t);-2*cos(t)]; %Enter the expression
qdDotDotDot = [0.125*sin(0.5*t);2*sin(t)]; %Enter the expression

% Parse current states (X is the same size and "form" as X0)
% (i.e., in this sim, XDot= [eDot;e2Dot;TfDot;Mu2Dot;thetaHatDot;YdfDot(1,:)';YdfDot(2,:)';Usgn];
e        = [X(1);X(2)];
e2       = [X(3);X(4)];
Tf       = [X(5);X(6)];
Mu2      = [X(7);X(8)];
thetaHat = [X(9);X(10);X(11);X(12);X(13)];
Ydf      = [X(14),X(15),X(16),X(17),X(18);X(19),X(20),X(21),X(22),X(23)];
Usgn     = [X(24);X(25)];

% Compute current x and xDot for convenience
q        = qd-e;
qDot     = -e2 + a1*e + qdDot;

% Compute cos(x2) and sin(x2) for convenience
c2       = cos(q(2));
s2       = sin(q(2));
% Compute cos(x2) and sin(x2) for convenience
cd2      = cos(qd(2));
sd2      = sin(qd(2));

% Compute current matrices for the dynamics
M        = [p1 + 2*p3*c2 p2 + p3*c2;p2 + p3*c2 p2];
V       = [-p3*s2*qDot(2) -p3*s2*(qDot(1) + qDot(2));p3*s2*qDot(1) 0];
fd       = [f1 0;0 f2];
Td       = [0.5*cos(0.5*t);sin(t)];
% %% Compute current regression matrix
% y11      = qdDot(1); %Enter the expression
% y12      = qdDot(2); %Enter the expression
% y13      = 2*c2*qdDot(1) + c2*qdDot(2) - qDot(2)*s2*(qDot(1) + qDot(2)) - qDot(1)*qDot(2)*s2; %Enter the expression
% y14      = qDot(1); %Enter the expression
% y15      = 0; %Enter the expression
% y21      = 0; %Enter the expression
% y22      = qdDot(1) + qdDot(2); %Enter the expression
% y23      = s2*qDot(1)^2 + c2*qdDot(1); %Enter the expression
% y24      = 0; %Enter the expression
% y25      = qDot(2); %Enter the expression
% Y        = [y11 y12 y13 y14 y15;y21 y22 y23 y24 y25];

%% Compute current regression matrix
yd11      = qdDotDot(1); %Enter the expression
yd12      = qdDotDot(2); %Enter the expression
yd13      = 2*cd2*qdDotDot(1) + cd2*qdDotDot(2) - qdDot(2)*sd2*(qdDot(1) + qdDot(2)) - qdDot(2)*qdDot(1)*sd2; %Enter the expression
yd14      = qdDot(1); %Enter the expression
yd15      = 0; %Enter the expression
yd21      = 0; %Enter the expression
yd22      = qdDotDot(1) + qdDotDot(2); %Enter the expression
yd23      = cd2*qdDotDot(1) + qdDot(1)*qdDot(1)*sd2; %Enter the expression
yd24      = 0; %Enter the expression
yd25      = qdDot(2); %Enter the expression
Yd        = [yd11 yd12 yd13 yd14 yd15;yd21 yd22 yd23 yd24 yd25];

%% Compute desired Dot regression matrix
ydDot11      = qdDotDotDot(1); %Enter the expression
ydDot12      = qdDotDotDot(2); %Enter the expression
ydDot13      = 2*cd2*qdDotDotDot(1)-2*sd2*qdDotDot(1)*qdDot(2) + cd2*qdDotDotDot(2)...
                -sd2*qdDotDot(2)*qdDot(2) - qdDot(2)*sd2*(qdDotDot(1) + qdDotDot(2))...
                -qdDot(2)*cd2*qdDot(2)*(qdDot(1) + qdDot(2))...
                -qdDotDot(2)*sd2*(qdDot(1) + qdDot(2)) - qdDotDot(2)*qdDot(1)*sd2...
                - qdDot(2)*qdDotDot(1)*sd2- qdDot(2)*qdDot(1)*cd2*qdDot(2); %Enter the expression
ydDot14      = qdDotDot(1); %Enter the expression
ydDot15      = 0; %Enter the expression
ydDot21      = 0; %Enter the expression
ydDot22      = qdDotDotDot(1) + qdDotDotDot(2); %Enter the expression
ydDot23      = -sd2*qdDot(2)*qdDotDot(1) +cd2*qdDotDotDot(1)+ qdDot(1)*qdDot(1)*cd2*qdDot(2)+2*qdDot(1)*qdDotDot(1)*sd2; %Enter the expression
ydDot24      = 0; %Enter the expression
ydDot25      = qdDotDot(2); %Enter the expression
YdDot        = [ydDot11 ydDot12 ydDot13 ydDot14 ydDot15;ydDot21 ydDot22 ydDot23 ydDot24 ydDot25];
eDot        = e2 - a1*e;
S = V*qDot+fd*q+M*qdDot+M*(a1*eDot+a2*e2)-Yd*theta;%+a1*M*e2-a1^2*M*e2+M*a2*e2

%% Design controller
TfHat     = Ydf*thetaHat+Mu2;
    if (t==0)
    e2init    = e2;
    end
 Mu1       = (K+1)*(e2-e2init)+Usgn; %Enter the expression
% Mu1=Usgn;
ep        = Tf - TfHat;
u         = Yd*thetaHat+Mu1;

end
function [XDot] = twoLinkdynamics(t,X,theta,gains)
global e2init
%% Parse for controller gains = [K;a1;a2;gamma;B;B2;k2]
K        = gains(1); %Enter a number
a1       = gains(2); %Enter a number
a2       = gains(3);
gamma    = gains(4:8);
gamma    = diag(gamma);
B        = gains(9);
B2       = gains(10);
k2       = gains(11);
%% Parse parameter vector
p1 = theta(1);
p2 = theta(2);
p3 = theta(3);
f1 = theta(4);
f2 = theta(5);


% Desired trajectory and needed derivatives
qd       = [cos(0.5*t);2*cos(t)];
qdDot    = [-0.5*sin(0.5*t);-2*sin(t)]; %Enter the expression
qdDotDot = [-0.25*cos(0.5*t);-2*cos(t)]; %Enter the expression
qdDotDotDot = [0.125*sin(0.5*t);2*sin(t)]; %Enter the expression

% Parse current states (X is the same size and "form" as X0)
% (i.e., in this sim, XDot= [eDot;e2Dot;TfDot;Mu2Dot;thetaHatDot;YdfDot(1,:)';YdfDot(2,:)';Usgn];
e        = [X(1);X(2)];
e2       = [X(3);X(4)];
Tf       = [X(5);X(6)];
Mu2      = [X(7);X(8)];
thetaHat = [X(9);X(10);X(11);X(12);X(13)];
Ydf      = [X(14),X(15),X(16),X(17),X(18);X(19),X(20),X(21),X(22),X(23)];
Usgn     = [X(24);X(25)];

% Compute current x and xDot for convenience
q        = qd-e;
qDot     = -e2 + a1*e + qdDot;

% Compute cos(x2) and sin(x2) for convenience
c2       = cos(q(2));
s2       = sin(q(2));
% Compute cos(x2) and sin(x2) for convenience
cd2      = cos(qd(2));
sd2      = sin(qd(2));

% Compute current matrices for the dynamics
M        = [p1 + 2*p3*c2 p2 + p3*c2;p2 + p3*c2 p2];
V       = [-p3*s2*qDot(2) -p3*s2*(qDot(1) + qDot(2));p3*s2*qDot(1) 0];
fd       = [f1 0;0 f2];
Td       = [0.5*cos(0.5*t);sin(t)];

%% Compute current regression matrix
yd11      = qdDotDot(1); %Enter the expression
yd12      = qdDotDot(2); %Enter the expression
yd13      = 2*cd2*qdDotDot(1) + cd2*qdDotDot(2) - qdDot(2)*sd2*(qdDot(1) + qdDot(2)) - qdDot(2)*qdDot(1)*sd2; %Enter the expression
yd14      = qdDot(1); %Enter the expression
yd15      = 0; %Enter the expression
yd21      = 0; %Enter the expression
yd22      = qdDotDot(1) + qdDotDot(2); %Enter the expression
yd23      = cd2*qdDotDot(1) + qdDot(1)*qdDot(1)*sd2; %Enter the expression
yd24      = 0; %Enter the expression
yd25      = qdDot(2); %Enter the expression
Yd        = [yd11 yd12 yd13 yd14 yd15;yd21 yd22 yd23 yd24 yd25];

%% Compute desired Dot regression matrix
ydDot11      = qdDotDotDot(1); %Enter the expression
ydDot12      = qdDotDotDot(2); %Enter the expression
ydDot13      = 2*cd2*qdDotDotDot(1)-2*sd2*qdDotDot(1)*qdDot(2) + cd2*qdDotDotDot(2)...
                -sd2*qdDotDot(2)*qdDot(2) - qdDot(2)*sd2*(qdDotDot(1) + qdDotDot(2))...
                -qdDot(2)*cd2*qdDot(2)*(qdDot(1) + qdDot(2))...
                -qdDotDot(2)*sd2*(qdDot(1) + qdDot(2)) - qdDotDot(2)*qdDot(1)*sd2...
                - qdDot(2)*qdDotDot(1)*sd2- qdDot(2)*qdDot(1)*cd2*qdDot(2); %Enter the expression
ydDot14      = qdDotDot(1); %Enter the expression
ydDot15      = 0; %Enter the expression
ydDot21      = 0; %Enter the expression
ydDot22      = qdDotDotDot(1) + qdDotDotDot(2); %Enter the expression
ydDot23      = -sd2*qdDot(2)*qdDotDot(1) +cd2*qdDotDotDot(1)+ qdDot(1)*qdDot(1)*cd2*qdDot(2)+2*qdDot(1)*qdDotDot(1)*sd2; %Enter the expression
ydDot24      = 0; %Enter the expression
ydDot25      = qdDotDot(2); %Enter the expression
YdDot        = [ydDot11 ydDot12 ydDot13 ydDot14 ydDot15;ydDot21 ydDot22 ydDot23 ydDot24 ydDot25];

eDot        = e2 - a1*e;
S = V*qDot+fd*q+M*qdDot+M*(a1*eDot+a2*e2)-Yd*theta;%+a1*M*e2-a1^2*M*e2+M*a2*e2

%% Design controller
TfHat     = Ydf*thetaHat+Mu2;
    if (t==0)
    e2init    = e2;
    end
 Mu1       = (K+1)*(e2-e2init)+Usgn; %Enter the expression
% Mu1=Usgn;
ep        = Tf - TfHat;
u         = Yd*thetaHat+Mu1;
r         = M\(S-u+Td+Yd*theta);
%% Compute current closed-loop dynamics
 Usgn        = (K+1)*a2*e2+B*sign(e2); 
% Usgn = (K+1)*r+B*sign(e2);
eDot        = e2 - a1*e;
e2Dot       = r - a2*e2;
TfDot       = -B*Tf+B*u;
Mu2Dot      = k2*ep+B2*sign(ep);
YdfDot      = -Ydf*B+Yd;
thetaHatDot = gamma*YdDot'*r+gamma*YdfDot'*ep; %Enter the expression

% Stacked dynamics vector (XDot is the same size and "form" as X)
XDot        = [eDot;e2Dot;TfDot;Mu2Dot;thetaHatDot;YdfDot(1,:)';YdfDot(2,:)';Usgn];
end