function twoLinkRobotAdaptive
clc
%Set up parameters for sim
p1       = 3.473;
p2       = 0.196;
p3       = 0.242;
f1       = 5.3;
f2       = 1.1;

% Stacked parameter vector
theta    = [p1;p2;p3;f1;f2];

% Simulation final time
tf   = 30;

% Initial condition vector (X0 must be same size and "form" as X and XDot below)
% (i.e., in this sim, X0 = [e0;r0;thetahat0])
X0   = 2*[ones(6,1)];

% Options for integration function
opts = odeset('RelTol',1e-3,'AbsTol',1e-3);

% Integrate (you can send the paramters theta to the dynamics as seen below)
[t,STATES] = ode45(@(t,X) twoLinkdynamics(t,X,theta),[0 tf],X0,opts);
u = zeros(2,length(t));
        for i=1:length(t)
        u(:,i) =  getu(t(i),STATES(i,:),theta);
        end
% Set up desired trajectory data for plots (enter desired trajectory for your simulation)
qd = [cos(0.5*t) 2*cos(t)]';

% Parse integrated states (STATES is the same "form" as X0)
% (i.e., in this sim, STATES = [e r thetahat] over all time);
e  = STATES(:,1:2)';
r  = STATES(:,3:4)';
% thetaHat = STATES(:,5:9)';

% Compute x from e and xd for plotting purposes
q  = qd-e;

% Plot the actual vs desired trajectories
figure(1)
plot(t,qd,'-','LineWidth',1)
hold on
ax = gca;
ax.ColorOrderIndex = 1;
plot(t,q,':','LineWidth',1)
title('Actual Vs Desred Trajectories')
xlabel('time') % x-axis label
ylabel('Trajectories') % y-axis label
legend('y = qd1','y = qd2','y=q1','y=q2')
hold off

% Plot the tracking error
figure(2)
plot(t,e,'--','LineWidth',2)
title('Tracking Error')
xlabel('time') % x-axis label
ylabel('Errors') % y-axis label
legend('y = e1','y = e2')

figure(3)
plot(t,u(1,:),'-',t,u(2,:),':','LineWidth',1)
title('Control Input')
xlabel('time') % x-axis label
ylabel('Torques') % y-axis label
legend('y = u1','y = u2')

end

function [u] = getu(t,X,theta)
global e2init
% Parse parameter vector
p1 = theta(1);
p2 = theta(2);
p3 = theta(3);
f1 = theta(4);
f2 = theta(5);

% Select gains for controller
K        = 5; 
a1       = 2; 
a2       = 2;
% Desired trajectory and needed derivatives
qd       = [cos(0.5*t);2*cos(t)];
qdDot    = [-0.5*sin(0.5*t);-2*sin(t)]; 
qdDotDot = [-0.25*cos(0.5*t);-2*cos(t)]; 

% Parse current states (X is the same size and "form" as X0)
% (i.e., in this sim, X = [e;r;thetahat])
e        = [X(1);X(2)];
e2       = [X(3);X(4)];
Usgn     = [X(5);X(6)];
% Compute current x and xDot for convenience
q        = qd-e;
qDot     = -e2 + a1*e + qdDot;

% Compute cos(x2) and sin(x2) for convenience
c2       = cos(q(2));
s2       = sin(q(2));

% Compute current matrices for the dynamics
M        = [p1 + 2*p3*c2 p2 + p3*c2;p2 + p3*c2 p2];
V        = [-p3*s2*qDot(2) -p3*s2*(qDot(1) + qDot(2));p3*s2*qDot(1) 0];
fd       = [f1 0;0 f2];
Td       = [0.5*cos(0.5*t);sin(t)];

S = V*qDot+fd*q+M*qdDot+a1*M*e2-a1^2*M*e2+M*a2*e2;
if (t==0)
    e2init = e2;
end
u        = (K+1)*(e2-e2init)+Usgn; 
end
function [XDot] = twoLinkdynamics(t,X,theta)
global e2init
% Parse parameter vector
p1 = theta(1);
p2 = theta(2);
p3 = theta(3);
f1 = theta(4);
f2 = theta(5);

% Select gains for controller
K        = 50; %Enter a number
a1       = 2; %Enter a number
a2       = 2;
B        = 3;
% Desired trajectory and needed derivatives
qd       = [cos(0.5*t);2*cos(t)];
qdDot    = [-0.5*sin(0.5*t);-2*sin(t)]; %Enter the expression
qdDotDot = [-0.25*cos(0.5*t);-2*cos(t)]; %Enter the expression

% Parse current states (X is the same size and "form" as X0)
% (i.e., in this sim, X = [e;r;thetahat])
e        = [X(1);X(2)];
e2       = [X(3);X(4)];
Usgn     = [X(5);X(6)];

% Compute current x and xDot for convenience
q        = qd-e;
qDot     = -e2 + a1*e + qdDot;

% Compute cos(x2) and sin(x2) for convenience
c2       = cos(q(2));
s2       = sin(q(2));

% Compute current matrices for the dynamics
M        = [p1 + 2*p3*c2 p2 + p3*c2;p2 + p3*c2 p2];
V        = [-p3*s2*qDot(2) -p3*s2*(qDot(1) + qDot(2));p3*s2*qDot(1) 0];
fd       = [f1 0;0 f2];
Td       = [0.5*cos(0.5*t);sin(t)];

S = V*qDot+fd*q+M*qdDot+a1*M*e2-a1^2*M*e2+M*a2*e2;

% Design controller
if (t==0)
    e2init = e2;
end
u         =(K+1)*(e2-e2init)+Usgn; %Enter the expression

r         = M\(S-u+Td);
%udot=(K+1)*r+B*sign(e2); 
% Compute current closed-loop dynamics
eDot       = e2 - a1*e;
e2Dot      = r - a2*e2;
UsgnDot    = (K+1)*e2+B*sign(e2);

% Stacked dynamics vector (XDot is the same size and "form" as X)
XDot        = [eDot;e2Dot;UsgnDot];
end
