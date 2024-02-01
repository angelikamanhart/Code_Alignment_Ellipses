clear all

% biological parameters
p.N = 125; %125; % number of initial ellipses
% cell sizes
p.rBar = 2; % equilibrium value of r
p.nu = 0.5; % non dimensional velocity parameter
p.gamma = 0.1; % non dimensional shape change parameter
p.kappa = 1; % non dimensional spring constant
p.mu = 5; % non dimensional bending stiffness
p.junction_range = 0.2; % non dimensional junction range (lambda)

p.sq = 20; % size of domain (square)

% which model (with or without cell-cell junctions)
%{
opt 0 .. No cell-cell junctions
opt 1 .. With cell-cell junctions (only at head and tail)
%}
p.optModel=0;

% what type of cell-cell junctions can form?
%{
opt 0 .. Any
opt 1 .. Only head-head or tail-tail
opt 2 .. Only head-tail
%}
p.optJunction=0;

% extra info recorded
%{
opt 0 .. No extra info
opt 1 .. Extra info on cell-cell junctions
%}
p.optTesting=0;

% varying or fixed cell shapes
%{
opt 0 .. Fixed cell shapes
opt 1 .. Varying cell shapes
%}
p.optCellShape=1;

% numerical parameters
p.n_t = 200; % number of theta points
p.search_range = 0.04; %0.9; % search range for overlapping ellipses

% define theta array
p.theta = linspace(0,2*pi-2*pi/p.n_t,p.n_t); % define theta points to discretise ellipse boundary

% simulation parameters
p.dt = 1/100; %timestep
p.T = 1; %max time

% INITIAL CONDITIONS
[xPos0,yPos0] = evenlySpacedPoints(p.N,p.sq,1);
dir0 = 2*pi*rand(1,p.N); % assign random ellipse directions
r0 = p.rBar*ones(1,p.N);

%% RUN SIMULATION

[x_centre, y_centre, alpha, r, junction_details] = fullSim(p,xPos0,yPos0,dir0,r0);

%% Visualisation

p.n_b = 8000; 
boundary = [zeros(p.n_b/4,1),linspace(0,p.sq,p.n_b/4)';linspace(0,p.sq,p.n_b/4)',p.sq*ones(p.n_b/4,1);
    p.sq*ones(p.n_b/4,1),flip(linspace(0,p.sq,p.n_b/4))';linspace(0,p.sq,p.n_b/4)',zeros(p.n_b/4,1)];

figure('Position',[200,200,650,600])
hold on

n_colors = 16;
for k = 1:10:p.T/p.dt
    p.a = sqrt(r(k,:));
    p.b = 1./sqrt(r(k,:));
    p.n_e = p.N;
    colormap(hsv(n_colors)); % map cells onto colourscale depending on orientation
    map = colormap;
    [N,edges,bin] = histcounts(mod(alpha(k,1:p.n_e),pi),linspace(0,pi,n_colors+1));
    clf
    for i = 1:p.n_e
        x_eb(:,i) = x_centre(k,i) + p.a(i)*cos(p.theta)*cos(alpha(k,i)) - p.b(i)*sin(p.theta)*sin(alpha(k,i));
        % y points on ellipses
        y_eb(:,i) = y_centre(k,i) + p.a(i)*cos(p.theta)*sin(alpha(k,i)) + p.b(i)*sin(p.theta)*cos(alpha(k,i));
        % plot ellipses
        plot(x_eb(:,i),y_eb(:,i),'color',map(bin(i),:),'linewidth',1.5);
        hold on
    end 
    switch p.optTesting
        case 1
            if ~isempty(junction_details{k,1})
                plot(junction_details{k,1}(1,:),junction_details{k,1}(2,:),'xr','linewidth',1)
%                 plot([x_eb(1,junction_details{k,3}(1,:))',x_eb(p.n_t/2+1,junction_details{k,3}(1,:))'].',[y_eb(1,junction_details{k,3}(1,:))',y_eb(p.n_t/2+1,junction_details{k,3}(1,:))'].','r','linewidth',1)
                hold on
                plot(junction_details{k,2}(1,:),junction_details{k,2}(2,:),'xb','linewidth',1)
%                 plot([x_eb(1,junction_details{k,3}(2,:))',x_eb(p.n_t/2+1,junction_details{k,3}(2,:))'].',[y_eb(1,junction_details{k,3}(2,:))',y_eb(p.n_t/2+1,junction_details{k,3}(2,:))'].','b','linewidth',1)
            end 
    end 
    plot(boundary(:,1),boundary(:,2))
    axis([-2,p.sq+2,-2,p.sq+2])
    set(gca,'fontsize',20)
    set(gcf,'color','w')
    c = colorbar('Ticks',[0,0.5,1],...
         'TickLabels',{'$0$', '$\pi/2$', '$\pi$'},...
         'TickLabelInterpreter', 'latex');
    c.Label.String = 'Orientation';
    pause(0.1)
end 

