function twoLinkRobotAdaptive
    close all
    %Set up parameters for sim
    p1       = 3.473;
    p2       = 0.196;
    p3       = 0.242;
    f1       = 5.3;
    f2       = 1.1;

    % Select gains for controller
    K        = .01; %Enter a number
    a        = .1875; %Enter a number
    gamma    = [30,20,40,90,20]';
    B        = 0.1875;
    gains = [K;a;gamma;B]
    % Stacked parameter vector
    theta    = [p1;p2;p3;f1;f2];

    % Simulation final time
    tf   = 60;

    % Initial condition vector (X0 must be same size and "form" as X and XDot below)
    % (i.e., in this sim, X0 = [e0;r0;thetahat0])
    X0   = [4;10;3;2;1;1;1;1;1;ones(12,1)];

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
    r  = STATES(:,3:4)';
    thetaHat = STATES(:,5:9)';

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
% Parse parameter vector
p1 = theta(1);
p2 = theta(2);
p3 = theta(3);
f1 = theta(4);
f2 = theta(5);



% Select gains for controller
K        = gains(1); %Enter a number
a        = gains(2); %Enter a number
gamma    = gains(3:7);
gamma    = diag(gamma);
B        = gains(8);
% Desired trajectory and needed derivatives
qd       = [cos(0.5*t);2*cos(t)];
qdDot    = [-0.5*sin(0.5*t);-2*sin(t)]; %Enter the expression
qdDotDot = [-0.25*cos(0.5*t);-2*cos(t)]; %Enter the expression

% Parse current states (X is the same size and "form" as X0)
% (i.e., in this sim, X = [e;r;thetahat])
e        = [X(1);X(2)];
r        = [X(3);X(4)];
thetaHat = [X(5);X(6);X(7);X(8);X(9)];
Ydf      = [X(10),X(11),X(12),X(13),X(14);X(15),X(16),X(17),X(18),X(19)];
uf       = [X(20);X(21)];
% Compute current x and xDot for convenience
q        = qd-e;
qDot     = -r + a*e + qdDot;

% Compute cos(x2) and sin(x2) for convenience
c2       = cos(q(2));
s2       = sin(q(2));

% Compute current matrices for the dynamics
M        = [p1 + 2*p3*c2 p2 + p3*c2;p2 + p3*c2 p2];
Vm       = [-p3*s2*qDot(2) -p3*s2*(qDot(1) + qDot(2));p3*s2*qDot(1) 0];
fd       = [f1 0;0 f2];

% Compute current regression matrix
y11      = qdDotDot(1) - a*qDot(1) + a*qdDot(1); %Enter the expression
y12      = qdDotDot(2) - a*qDot(2) + a*qdDot(2); %Enter the expression
y13      = 2*c2*qdDotDot(1) + c2*qdDotDot(2) - a*(qd(2)*s2*(qDot(1) + qDot(2)) + qDot(2)*qd(1)*s2) - qdDot(2)*s2*(qDot(1) + qDot(2)) - 2*a*c2*qDot(1) - a*c2*qDot(2) + 2*a*c2*qdDot(1) + a*c2*qdDot(2) - qDot(2)*qdDot(1)*s2 + a*q(1)*qDot(2)*s2 + a*q(2)*s2*(qDot(1) + qDot(2)); %Enter the expression
y14      = qDot(1); %Enter the expression
y15      = 0; %Enter the expression
y21      = 0; %Enter the expression
y22      = qdDotDot(1) + qdDotDot(2) - a*qDot(1) - a*qDot(2) + a*qdDot(1) + a*qdDot(2); %Enter the expression
y23      = c2*qdDotDot(1) - a*c2*qDot(1) + a*c2*qdDot(1) + qDot(1)*qdDot(1)*s2 - a*q(1)*qDot(1)*s2 + a*qDot(1)*qd(1)*s2; %Enter the expression
y24      = 0; %Enter the expression
y25      = qDot(2); %Enter the expression
Y        = [y11 y12 y13 y14 y15;y21 y22 y23 y24 y25];
% Compute cos(x2) and sin(x2) for convenience
cd2       = cos(qd(2));
sd2       = sin(qd(2));
% Compute current regression matrix
yd11      = qdDotDot(1); %Enter the expression
yd12      = qdDotDot(2); %Enter the expression
yd13      = 2*cd2*qdDotDot(1) + cd2*qdDotDot(2) - qdDot(2)*sd2*(qDot(1) + qDot(2)) - qDot(2)*qdDot(1)*sd2; %Enter the expression
yd14      = qDot(1); %Enter the expression
yd15      = 0; %Enter the expression
yd21      = 0; %Enter the expression
yd22      = qdDotDot(1) + qdDotDot(2); %Enter the expression
yd23      = cd2*qdDotDot(1) + qDot(1)*qdDot(1)*sd2; %Enter the expression
yd24      = 0; %Enter the expression
yd25      = qDot(2); %Enter the expression
Yd        = [yd11 yd12 yd13 yd14 yd15;yd21 yd22 yd23 yd24 yd25];


% Design controller

u         = K*r+Y*thetaHat+e; %Enter the expression
end
function [XDot] = twoLinkdynamics(t,X,theta,gains)

% Parse parameter vector
p1 = theta(1);
p2 = theta(2);
p3 = theta(3);
f1 = theta(4);
f2 = theta(5);

% Select gains for controller
K        = gains(1); %Enter a number
a        = gains(2); %Enter a number
gamma    = gains(3:7);
gamma    = diag(gamma);
B        = gains(8);
% Desired trajectory and needed derivatives
qd       = [cos(0.5*t);2*cos(t)];
qdDot    = [-0.5*sin(0.5*t);-2*sin(t)]; %Enter the expression
qdDotDot = [-0.25*cos(0.5*t);-2*cos(t)]; %Enter the expression

% Parse current states (X is the same size and "form" as X0)
% (i.e., in this sim, X = [e;r;thetahat])
e        = [X(1);X(2)];
r        = [X(3);X(4)];
thetaHat = [X(5);X(6);X(7);X(8);X(9)];
Ydf      = [X(10),X(11),X(12),X(13),X(14);X(15),X(16),X(17),X(18),X(19)];
uf       = [X(20);X(21)];
% Compute current x and xDot for convenience
q        = qd-e;
qDot     = -r + a*e + qdDot;

% Compute cos(x2) and sin(x2) for convenience
c2       = cos(q(2));
s2       = sin(q(2));

% Compute current matrices for the dynamics
M        = [p1 + 2*p3*c2 p2 + p3*c2;p2 + p3*c2 p2];
Vm       = [-p3*s2*qDot(2) -p3*s2*(qDot(1) + qDot(2));p3*s2*qDot(1) 0];
fd       = [f1 0;0 f2];
Td       = [0.5*cos(0.5*t);sin(t)];
% Compute current regression matrix
y11      = qdDotDot(1) - a*qDot(1) + a*qdDot(1); %Enter the expression
y12      = qdDotDot(2) - a*qDot(2) + a*qdDot(2); %Enter the expression
y13      = 2*c2*qdDotDot(1) + c2*qdDotDot(2) - a*(qd(2)*s2*(qDot(1) + qDot(2)) + qDot(2)*qd(1)*s2) - qdDot(2)*s2*(qDot(1) + qDot(2)) - 2*a*c2*qDot(1) - a*c2*qDot(2) + 2*a*c2*qdDot(1) + a*c2*qdDot(2) - qDot(2)*qdDot(1)*s2 + a*q(1)*qDot(2)*s2 + a*q(2)*s2*(qDot(1) + qDot(2)); %Enter the expression
y14      = qDot(1); %Enter the expression
y15      = 0; %Enter the expression
y21      = 0; %Enter the expression
y22      = qdDotDot(1) + qdDotDot(2) - a*qDot(1) - a*qDot(2) + a*qdDot(1) + a*qdDot(2); %Enter the expression
y23      = c2*qdDotDot(1) - a*c2*qDot(1) + a*c2*qdDot(1) + qDot(1)*qdDot(1)*s2 - a*q(1)*qDot(1)*s2 + a*qDot(1)*qd(1)*s2; %Enter the expression
y24      = 0; %Enter the expression
y25      = qDot(2); %Enter the expression
Y        = [y11 y12 y13 y14 y15;y21 y22 y23 y24 y25];
% Compute cos(x2) and sin(x2) for convenience
cd2       = cos(qd(2));
sd2       = sin(qd(2));
% Compute current regression matrix
yd11      = qdDotDot(1); %Enter the expression
yd12      = qdDotDot(2); %Enter the expression
yd13      = 2*cd2*qdDotDot(1) + cd2*qdDotDot(2) - qdDot(2)*sd2*(qDot(1) + qDot(2)) - qDot(2)*qdDot(1)*sd2; %Enter the expression
yd14      = qDot(1); %Enter the expression
yd15      = 0; %Enter the expression
yd21      = 0; %Enter the expression
yd22      = qdDotDot(1) + qdDotDot(2); %Enter the expression
yd23      = cd2*qdDotDot(1) + qDot(1)*qdDot(1)*sd2; %Enter the expression
yd24      = 0; %Enter the expression
yd25      = qDot(2); %Enter the expression
Yd        = [yd11 yd12 yd13 yd14 yd15;yd21 yd22 yd23 yd24 yd25];


% Design controller

u         = K*r+Y*thetaHat+e; %Enter the expression
ep        = uf- Ydf*thetaHat;
% Compute current closed-loop dynamics
eDot        = r - a*e;
rDot        = M\(Y*theta-u-Vm*r+Td); %Enter the expression
thetaHatDot = gamma*Y'*r+gamma*Ydf'*ep; %Enter the expression
YdfDot      = -Ydf*B+Yd;
ufDot       = -B*uf+u;
% Stacked dynamics vector (XDot is the same size and "form" as X)
XDot        = [eDot;rDot;thetaHatDot;YdfDot(1,:)';YdfDot(2,:)';ufDot];
end