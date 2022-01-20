%% Variables
global K;
global obstacles;
global goal;
global world;
global start;
global maxSteps;
global epsilon;
global mag;

K = 3;
obstacles = [struct('center',[-8;5],'radius',5,'distanceInfluence',4); struct('center',[0;-9],'radius',4,'distanceInfluence',20); struct('center',[5;0],'radius',4,'distanceInfluence',20)];
goal = [15;0];
start = [-5,-15,-10,0,-5,-10;-5,5,12,3,-15,15];
world = struct('center',[0;0], 'radius', 20,'distanceInfluence',2);
maxSteps = 60000;
epsilon = 1;
mag = 0.1;



close all;

%% Plot Contour w/path
plot_contour(true);

%% Plot surface
plot_surface();

%% Plotting
function plot_contour(withPath)
global obstacles;
global goal;

x = linspace(-20, 20);
y =linspace(-20, 20);
[xx, yy] = meshgrid(x, y);
UNav = get_potential(xx, yy);

figure
contour(xx, yy, UNav)
colormap winter

hold on

for i=1:numel(obstacles)
    obstacle = obstacles(i);
    plot(obstacle.center(1,1), obstacle.center(2,1),'b');
    hold on
end

plot(goal(1),goal(2),"r*");
hold on

if withPath
    plotPath()
end
end

function plotPath()
global start;
pathColors = ['r', 'b', 'g', 'm','c','k','y'];
if size(start,2) > 0
    for i = 1:size(start,2)
        start_ = start(:,i);
        path = calculate_path(start_);
        color = ".-"+pathColors(i);
        plot(path(1,:), path(2,:),color,'linewidth',2);
        plot(start_(1), start_(2), "b*");
        hold on
    end
else
    path = calculate_path(start);
    plot(path(1,:), path(2,:), 'b.-');
    plot(start(1), start(2), "b*");
    plot(path(1,:), path(2,:), '-r','linewidth',2)
end
end

function plot_surface()
figure
x = linspace(-20, 20);
y =linspace(-20, 20);
[xx, yy] = meshgrid(x, y);
UNav = get_potential(xx, yy);

surf(xx, yy, UNav,'FaceColor','interp','EdgeColor','interp')

end

function [UNav] = get_potential(xx,yy)
UNav = [];
for i = 1:numel(xx)
    qx = xx(i);
    qy = yy(i);
    q = [qx; qy];
    phi = calculate_phi(q);
    UNav(i) = phi;
end
UNav =reshape(UNav, size(xx));

end

%% Calculate path
function [path] = calculate_path(start_)
global maxSteps;
global epsilon;
global mag;
global start;
global goal;

if isnan(start_)
    start_ = start;
end

startingDistance = distance(start_, goal);
path = zeros(2, maxSteps);
path(:,1) = start_;

for step = 2:maxSteps
    q = path(:,step - 1);
    phiGrad = calculate_phiGrad(q);
    dist_toGoal = min(distance(q, goal), startingDistance);
     
    if ~isnan(mag)
        eta = (dist_toGoal / startingDistance)^2;
        phiGrad = eta * mag * phiGrad / norm(phiGrad);
    end
    path(:,step) = path(:,step - 1) - epsilon * phiGrad;
    
    if dist_toGoal < 0.005
        path = path(:,1:step);
        break;
    end
end
end

%% Calculate distance between 2 points
function [dist] = distance(q1, q2)
dist = norm(q1 - q2);
end

%% Calculate gamma
function [gamma] = calculate_gamma(q)
global K;
global goal;
gamma = distance(q, goal)^(2 * K);
end

function [gammaGrad] = calculate_gammaGrad(q)
global K;
global goal;

gammaGrad = 2 * K * distance(q, goal)^(2 * K - 1) * (q - goal) / distance(q, goal);
end

%% Calculate beta and its gradient
function [beta] = calculate_beta(q)
global world;
global obstacles;

beta = world.radius^2 - distance(q, world.center)^2;
for i = 1:numel(obstacles)
    obstacle = obstacles(i);
    beta = beta * (distance(q, obstacle.center)^2 - obstacle.radius^2);
end
end

function [betaGrad] = calculate_betaGrad(q)
global world;
global obstacles;

num_ofObstacles = numel(obstacles);

beta_0 = world.radius^2 - distance(q, world.center)^2;
betaGrad_0 = -2 * (q - world.center);

betas = ones(1,num_ofObstacles + 1);
betasGrad = ones(2,num_ofObstacles + 1);

betas(1,1) = beta_0;
betasGrad(:,1) = betaGrad_0;

for i = 1:num_ofObstacles
    obstacle = obstacles(i);
    betas(:,i) = distance(q, obstacle.center)^2 - obstacle.radius^2;
    betasGrad(:,i) = 2 * (q - obstacle.center);
end

betaGrad = zeros(2,1);
for i = 1:size(betas,2)
    temp = betasGrad(:,i);
    for j = 1:size(betas,2)
        if j ~= i
            temp = temp * betas(1,j);
        end
    end
    
    betaGrad = betaGrad + temp;
end
end

%% Calculate alpha and its gradient
function [alpha] = calculate_alpha(q)
gamma = calculate_gamma(q);
beta = calculate_beta(q);
alpha = gamma / beta;
end

function [alphaGrad] = calculate_alphaGrad(q)
gamma = calculate_gamma(q);
beta = calculate_beta(q);
gammaGrad = calculate_gammaGrad(q);
betaGrad = calculate_betaGrad(q);
alphaGrad = (gammaGrad * beta - gamma * betaGrad) / beta^2;
end

%% Calculate phi and its gradient
function [phi] = calculate_phi(q)
global K

alpha = calculate_alpha(q);
if alpha < 0
    phi = 1;
else
    phi = (alpha / (1 + alpha))^(1/K);
end
end

function [phiGrad] = calculate_phiGrad(q)
global K;
alpha = calculate_alpha(q);
alphaGrad = calculate_alphaGrad(q);
phiGrad = (1 / K) * (alpha / (1 + alpha))^((1 - K) / K) * (1 / (1 + alpha)^2) * alphaGrad;

end
