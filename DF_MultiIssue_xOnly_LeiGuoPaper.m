clear all
close all
clc

%This script studies the original DF model under switching topologies,
%where each topology is strongly connected.

%% Simulation parameters
t_final = 99;   %Total number of time steps-1 (29 for 30 time steps)
t = 0:1:t_final; %Time-steps for a single issue
n = 4;   % Number of agents
s = 60;   % Number of issues
plot_y = 0; % Plot individual issues
Nst = 3;  %The number of topologies in the switching set

%% Generate topology switching signal
%This generates a random switching signal. The + X shifts the topologies in the case number
% swt_sig = randi(Nst,[1,s])+0;

%Generated from swt_sig = randi(Nst,[1,s])+0;
%Time varying 1 in JSSC paper
% swt_sig = [1,1,1,3,2,2,1,3,2,3,2,2,3,3,2,3,3,2,2,2,3,3,2,3,3,2,2,1,1,1,1,2,1,3,2,1,1,1,2,1,1,3,2,3,2,2,1,1,1,2,1,3,3,3,2,1,1,2,3,2];

%Time varying 2 in JSSC paper
% swt_sig = [3,2,2,3,1,3,2,1,3,1,1,2,1,1,2,2,2,2,1,1,1,1,3,1,3,3,3,3,1,1,2,1,3,3,2,1,3,1,2,2,2,3,3,3,1,2,1,3,2,3,3,3,2,1,2,1,1,1,1,3];


%% This generates a periodic switch, with an adjustable period of issues
period = 1; %The number of issues between switches
swt = [2,1,3,1];  %The order of the topology switching  % Periodic 1 in JSSC paper
% swt = [3,2,1,3];  %The order of the topology switching    % Periodic 2 in JSSC paper
per_swt = kron(swt,ones(1,period));    %Index of topologies for a single period
num_per = ceil(s/length(per_swt));    %Number of periods for the given number of issues s
swt_sig = kron(ones(1,num_per),per_swt);  %Generate periodic switching signal
swt_sig = swt_sig(1:s);  %This truncates the signal if it is too long

%% Scenario parameters

[C,c_vect] = influence_matrix_switch_LeiGuo(swt_sig(1));   %Initial

w_ii_0 = [0.9;0;0.05;0.05];   %Initial condition set 1 in JSSC paper
% w_ii_0 = [0.1;0;0.8;0.1];   %Initial condition set 2 in JSSC paper



%% Influence dynamics

w_ii = zeros(n,s);   %Generate storage for self-influence weights

w_ii(:,1) = w_ii_0;   %The initial self-influence (i.e. for issue 1)

for j = 1:s-1
    
    %Update the adjucency matrix according to switching rule
    [C,c_vect] = influence_matrix_switch_LeiGuo(swt_sig(j));
    
    % Update self-weights using dominant left eigenvector
    for i = 1:n
        u(i,j) = c_vect(i)/(1-w_ii(i,j));  %w_ii(j) is x(s)
    end
    alpha(j) = 1/(sum(u(:,j)));   %\alpha(x(s))
    for i = 1:n
        F(i,j) = alpha(j)*u(i,j);  %F(x(s))
    end
    w_ii(:,j+1) = F(:,j);   % x(s+1) = F(x(s))
    
    
end

%% Plots

figure
hold on
plot(w_ii(1,:),':r','LineWidth',2)
plot(w_ii(2,:),':b','LineWidth',2)
plot(w_ii(3,:),':c','LineWidth',2)
plot(w_ii(4,:),':k','LineWidth',2)
ylabel('Self-confidence, x_i(s)')
xlabel('Topic, s')
axis([1 20 0 1])

plot(w_ii(1,:),'r','LineWidth',2)
plot(w_ii(2,:),'b','LineWidth',2)
plot(w_ii(3,:),'c','LineWidth',2)
plot(w_ii(4,:),'k','LineWidth',2)
