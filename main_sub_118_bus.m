%% Setting up environment

startup_project
config_file

% Load subsystem models
load linSys1.mat
load linSys2.mat
load linSys3.mat

% Load equilibrium data
load nominalU.mat
load nominalY.mat
load nominalInter.mat
P_input = 0;

file_name = 'case118_system_s0';
sim_name = "subsystem_118_bus";

% Base case load flow
[sys_data, branch_data] = loaddata(file_name);

T = incidencematrix(branch_data);
mpc = mydata2matpower(sys_data, branch_data);
mpopt = mpoption('out.all', 1,'verbose',1);
results = runpf(mpc,mpopt);

%% Initialise simulink models

% Calculate shunt susceptances at each bus for areas
branch_data_reshaped = [branch_data(:,1) branch_data(:,5)/2;
                        branch_data(:,2) branch_data(:,5)/2];
[unique_keys, ~, key_indices] = unique(branch_data_reshaped(:, 1));
sys_data(sys_data(:,1) == unique_keys, size(sys_data,2)+1) ...
    = accumarray(key_indices, branch_data_reshaped(:, 2)); % shunt susceptances at each bus


%% Subsystem model
sub = struct;
bus_area = sys_data(:, end-1); % area in the second last column
branch_area = branch_data(:, end);
areas = unique(bus_area);
na = length(areas);
temp_branch = bus_area(branch_data(:,1:2));
for i = 1:na
    sub(i).branch_data = branch_data(temp_branch(:,1) == i & temp_branch(:,2) == i,:);
    temp_T = incidencematrix(sub(i).branch_data);
    sub(i).T = temp_T(sum(boolean(temp_T),2) ~= 0, :);
    sub(i).sys_data = sys_data(bus_area==i,:);
    sub(i).bus_data = results.bus(bus_area==i,:);
    [~, ~, igen] = intersect(sub(i).sys_data(:,1),results.gen(:,1));
    sub(i).gen_data = results.gen(igen,:);

    [sub(i).n_devices, sub(i).busnum, sub(i).reordering, sub(i).device_selector, ...
        sub(i).total_devices, sub(i).infbus, sub(i).float, sub(i).load1, sub(i).gfl, ...
        sub(i).gfm, sub(i).line] = sub_assignvariables_39bus(sub(i).sys_data, ...
        sub(i).branch_data, sub(i).gen_data, sub(i).bus_data, sub(i).T, w, i, sim_name);

    sub(i).n_gfl = sub(i).n_devices(4);
    sub(i).n_gfm = sub(i).n_devices(5);
    
    Tstep_gfl = 0.1;
    Tstop = 1;
    if sub(i).n_gfl > 0
        sub(i).gfl_gains = design_gains_gfl_simple(sub(i).sys_data, 20, 0.9, 20, 20, 0.5); % PLL BW  zeta fp fv xg %
        gfl_sel = [];  % Select from 1:n_gfl
        sub(i).V_before_gfl = ones(sub(i).n_gfl,1); sub(i).V_after_gfl = sub(i).V_before_gfl; % V_after_gfl(gfl_sel) = 1
        sub(i).P_before_gfl = ones(sub(i).n_gfl,1); sub(i).P_after_gfl = sub(i).P_before_gfl; % P_after_gfl(gfl_sel) = 0.99
    end
    
    Tstep_gfm = 0.1;
    if sub(i).n_gfm > 0
        sub(i).gfm_gains = design_gains_gfm_simple(sub(i).sys_data, sub(i).gfm, 50);  % fv
        gfm_sel = [];  % Select from 1:n_gfm
        sub(i).V_before_gfm = ones(sub(i).n_gfm,1); sub(i).V_after_gfm = sub(i).V_before_gfm;
        sub(i).P_before_gfm = ones(sub(i).n_gfm,1); sub(i).P_after_gfm = sub(i).P_before_gfm;
    end
end

%% Connection line model
con_size = nchoosek(length(areas),2);
con_matrix = [];
for i_area = 1:length(areas)
    con_matrix = [con_matrix; [i_area*ones(size(i_area+1:length(areas))); i_area+1:length(areas)]'];
end

for i_con = 1:con_size
    ii = con_matrix(i_con, 1);
    jj = con_matrix(i_con, 2);
    con.branch_data = branch_data((temp_branch(:,1) == ii & temp_branch(:,2) == jj) | (temp_branch(:,1) == jj & temp_branch(:,2) == ii),:);
    temp_T = incidencematrix(con.branch_data);
    con.T = temp_T(sum(boolean(temp_T),2) ~= 0, :);
    con.nbus = size(con.T,1);
    temp_sys = reshape(con.branch_data(:,1:2), [], 1);
    temp_branch1 = reshape(bus_area(con.branch_data(:,[1,2])), [], 1);
    
    temp_con_bus = [];
    for i = 1:na
        sub(i).con_bus{i_con} = temp_sys(temp_branch1 == i);
        temp_con_bus = [temp_con_bus; sub(i).con_bus{i_con}];
        [~, ~, sub(i).vout{i_con}] = intersect(sub(i).con_bus{i_con}, sub(i).sys_data(:,1),'stable');
        sub(i).I_in{i_con} = zeros(size(sub(i).sys_data,1), length(sub(i).vout{i_con}));
        for j = 1:length(sub(i).vout{i_con})
            sub(i).I_in{i_con}(sub(i).vout{i_con}(j), j) = 1;
        end
    end

    [con_bus, con.reorder_v] = sort(temp_con_bus); 

    for i = 1:na
        [~, ~, con.I_select{i}] = intersect(sub(i).con_bus{i_con}, con_bus, 'stable');
    end

    % Putting Loadflow Results in V|delta|P|Q form
    % for floating buses, I need to insert zero value for P and Q in the matrix below
    % define_constants % matpower constants are being named
    PQ_at_each_bus = zeros(size(sys_data,1), 2); % initializing a matrix with containing P and Q values for each bus
    PQ_at_each_bus(results.gen(:,1),:) = results.gen(:,[2,3]); % inserting P and Q values of generator buses. Floating buses remain zero
    V_Delta_P_Q = [results.bus(:, 8) results.bus(:, 9)*pi/180  PQ_at_each_bus]; % forming a matrix for each bus containing solved V, delta, P, Q

    r_line = con.branch_data(:,3);
    x_line = con.branch_data(:,4);

    v_bus_dq = V_Delta_P_Q(con_bus,1); % q component is zero in their local frames
    theta_bus = V_Delta_P_Q(con_bus,2); % theta of each bus in radians
    v_bus_DQ = v_bus_dq.*exp(1j*theta_bus);
    v_diff_DQ = con.T'*v_bus_DQ;
    i_line_DQ = v_diff_DQ./(r_line+1j*x_line);

    % LINE
    % constants
    % initial conditions
    con.line.i_D = real(i_line_DQ);
    con.line.i_Q = imag(i_line_DQ);
    % multipliers
    con.line.r = r_line;
    con.line.l = x_line/w;
    con_sub(i_con) = con;
    clear con;
end

% clear sys_data branch_data
  
%% obtaining system model

Gsys1 = linearize(sim_name);
ev_G = eig(Gsys1)/2/pi;
si1 = any(real(ev_G)>1e-6)
max(real(ev_G))