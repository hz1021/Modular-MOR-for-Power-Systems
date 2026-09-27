function mpc = mydata2matpower(sys_data, branch_data)
%   Detailed explanation goes here
mpc.version = '2';
mpc.baseMVA = 1;

bus_data = sys_data(:,1:6);
bus_data(isnan(bus_data)) = 0;
gen_buses = sys_data(:,7)==10 | sys_data(:,7)==100 | sys_data(:,7)==200 | sys_data(:,7)==50;
gen_buses = sys_data(gen_buses,1);

slack_bus = sys_data(sys_data(:,2)==3,1);
load_buses = sys_data(sys_data(:,7)==12,1);


n = size(bus_data,1); % total number of buses
m = size(branch_data,1); % total number of branches
n_gen = sum(sys_data(:,7)==10) + sum(sys_data(:,7)==100) + sum(sys_data(:,7)==200) + sum(sys_data(:,7)==50);
% counting infinite bus, gfl, gfm, and sg (sg not included yet)
% 
mpc.bus = zeros(n,13); % initializing a bus matrix with zeros 
mpc.bus(:,1:2) = bus_data(:,[1,2]); % copying first two columns from bus data containing the bus number and it's type
mpc.bus(load_buses,3:4) = -bus_data(load_buses,[5,6]); % copying load bus data. since it is entered negative in excel sheet we need to multiply it by 1
mpc.bus(:,8) = 1; % making voltage magnitude at each bus equal to one
mpc.bus(slack_bus,[8,9]) = bus_data(slack_bus,[3,4]); % copying slack bus voltage and angle
mpc.bus(:,[7,10,11]) = 1; % making zone, baseKV and area equal to 1
mpc.bus(:,12:13) = [2*ones(n,1) 0.5*ones(n,1)]; % maxVm and minVm
% 
mpc.gen = zeros(n_gen,21);  % initializing a generator matrix with zeros
mpc.gen(:,[1,2,3,6]) = [bus_data(gen_buses,1) bus_data(gen_buses,5) bus_data(gen_buses,6) bus_data(gen_buses,3)]; % copying bus number, Pg, Qg, and Vg
mpc.gen(:,[4,5,7,8,9,10]) = [5*ones(n_gen,1) -5*ones(n_gen,1) ones(n_gen,2) 100*ones(n_gen,1) -100*ones(n_gen,1)]; % setting Qmax = 5, 
% Qmin = -5, mbase = 1, status = 1, Pmax = 5, Pmin = -5
% 
% 
% 
mpc.branch = zeros(m,13); % initializing branch matrix with zeros
mpc.branch(:,1:5) = branch_data(:,1:5); % copying from_bus, to_bus, r, x, bc
mpc.branch(:,[9,11]) = 1; % transformer tap and status = 1

end