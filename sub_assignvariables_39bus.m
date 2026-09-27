function [n_devices, busnum, reordering, device_selector, total_devices, ...
    inf, float, load, gfl, gfm, line] = sub_assignvariables_39bus(sys_data, ...
    branch_data, gen_data, bus_data, T, w, sub_i, file1)

define_constants %matpower constants are being named
n = size(sys_data,1); % total number of buses
open(file1+".slx");

%% Getting Bus Organizing Variables: busnum, reordering, device_selector

device_data = sys_data(:,[1,2,7]);
temp_bus = (1:n)';

busnum.inf = temp_bus(device_data(:,3)==10,1); % getting the infinite bus number
busnum.float = temp_bus(device_data(:,3)==11,1); % getting the floating bus number
busnum.load = temp_bus(device_data(:,3)==12,1); % getting the load bus number
busnum.gfl = temp_bus(device_data(:,3)==100,1); % getting the gfl pv number
busnum.gfm = temp_bus(device_data(:,3)==200,1); % getting the gfm number

% the vector below is very important. It corresponds to the muxing I do in Simulink file named toolbox_analysis.slx
incoming = [busnum.inf; busnum.float; busnum.load; busnum.gfl; busnum.gfm;];
to_sort = [incoming, (1:n)'];
sss = sortrows(to_sort);
reordering = sss(:,2);

% now getting device availability
device_availability = zeros(5,1); % inf, float, load, gfl, gfm
fields = fieldnames(busnum);
for i = 1:length(fields)
    if isempty(busnum.(fields{i}))
        busnum.(fields{i}) = 0;
    else
        device_availability(i) = numel(busnum.(fields{i}));
    end
end
n_inf = device_availability(1);
n_float = device_availability(2);
n_load = device_availability(3);
n_gfl = device_availability(4);
n_gfm = device_availability(5);

n_devices = [n_inf, n_float, n_load, n_gfl, n_gfm];

device_nonibr = boolean(device_availability(1:3)');
device_ibr = boolean(device_availability(4:end)');

%%%%%%%%%%%%%%%%%%%%%%%%%%% non-ibr buses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isequal(device_nonibr, [0 0 1]) % gfm available only
   ndg = 1:(1+1+n_load); % ndg stands for number of devices gross
   device_selector.nonibr = ndg(3:end);
   total_devices.nonibr = numel(ndg);
end
if isequal(device_nonibr, [0 1 0]) % gfl-pq available only
   ndg = 1:(1+n_float+1); % ndg stands for number of devices gross
   device_selector.nonibr = ndg(2:end-1);
   total_devices.nonibr = numel(ndg);
end
if isequal(device_nonibr, [0 1 1]) % gfl-pq and gfm available
   ndg = 1:(1+n_float+n_load); % ndg stands for number of devices gross
   device_selector.nonibr = ndg(2:end);
   total_devices.nonibr = numel(ndg);
end
if isequal(device_nonibr, [1 0 0]) % gfl-pv available only
   ndg = 1:(n_inf+1+1); % ndg stands for number of devices gross
   device_selector.nonibr = ndg(1:end-2);
   total_devices.nonibr = numel(ndg);
end
if isequal(device_nonibr, [1 0 1]) % gfl-pv and gfm available
   ndg = 1:(n_inf+1+n_load); % ndg stands for number of devices gross
   device_selector.nonibr = [ndg(1:n_inf), ndg(end-n_load+1:end)];
   total_devices.nonibr = numel(ndg);
end
if isequal(device_nonibr, [1 1 0]) % gfl-pv and gfl-pq available
   ndg = 1:(n_inf+n_float+1); % ndg stands for number of devices gross
   device_selector.nonibr = ndg(1:end-1);
   total_devices.nonibr = numel(ndg);
end
if isequal(device_nonibr, [1 1 1]) % all available
   ndg = 1:(n_inf+n_float+n_load); % ndg stands for number of devices gross
   device_selector.nonibr = ndg;
   total_devices.nonibr = numel(ndg);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% ibr buses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isequal(device_ibr, [0 1]) % gfm available only
   ndg = 1:(1+n_gfm); % ndg stands for number of devices gross
   device_selector.ibr = ndg(2:end);
   total_devices.ibr = numel(ndg);
end

if isequal(device_ibr, [1 0]) % gfl-pv available only
   ndg = 1:(n_gfl+1); % ndg stands for number of devices gross
   device_selector.ibr = ndg(1:end-1);
   total_devices.ibr = numel(ndg);
end

if isequal(device_ibr, [1 1]) % gfl-pv and gfm available
   ndg = 1:(n_gfl+n_gfm); % ndg stands for number of devices gross
   device_selector.ibr = ndg;
   total_devices.ibr = numel(ndg);
end
%% Calculating all bus quantities v_bus_DQ, i_bus_DQ, i_inj_DQ, c_bus, theta_bus, p_ref, q_ref, v_ref and all line quantities, i_line_DQ, r_line, x_line, l_line

% Putting Loadflow Results in V|delta|P|Q form
% for floating buses, I need to insert zero value for P and Q in the matrix below
PQ_at_each_bus = zeros(n,2); % initializing a matrix with containing P and Q values for each bus
[~, igen, ~] = intersect(bus_data(:,1), gen_data(:,1));
PQ_at_each_bus(igen,:) = gen_data(:,[2,3]); % inserting P and Q values of generator buses. Floating buses remain zero
V_Delta_P_Q = [bus_data(:, VM) bus_data(:, VA)*pi/180  PQ_at_each_bus]; % forming a matrix for each bus containing solved V, delta, P, Q

r_line = branch_data(:,3);
x_line = branch_data(:,4);

p_ref = V_Delta_P_Q(:, 3);
q_ref = V_Delta_P_Q(:, 4);
v_ref = V_Delta_P_Q(:, 1);

v_bus_dq = V_Delta_P_Q(:,1); % q component is zero in their local frames
theta_bus = V_Delta_P_Q(:,2); % theta of each bus in radians
v_bus_DQ = v_bus_dq.*exp(1j*theta_bus);
v_diff_DQ = T'*v_bus_DQ;
i_line_DQ = v_diff_DQ./(r_line+1j*x_line);
i_bus_DQ = T*i_line_DQ;

%%%%%%%%%%%%%% determining susceptances at all buses c_bus %%%%%%%%%%%%%%%%%%%
b_bus = zeros(n,1); % initializing Bs at all buses
% for ibr buses, we simply get the data from sys_data's 8th column
if n_gfl > 0
    b_bus(busnum.gfl)=sys_data(busnum.gfl,8);
end
if n_gfm > 0
    b_bus(busnum.gfm)=sys_data(busnum.gfm,8);
end

% for non IBR buses shunt susceptatnce will be obtained from branch data
b_ss = sys_data(:, end); % check if shunt susceptances are in the last column
b_bus = b_bus + b_ss;
b_bus(b_bus==0) = 0.1;
c_bus = b_bus/w;

% Now I can find injected current at each bus %
i_inj_DQ = 1j.*b_bus.*v_bus_DQ + i_bus_DQ;
i_inj_dq = i_inj_DQ.*exp(-1j*theta_bus);



%% string names
str_inf = "/Subsystem" + num2str(sub_i) + "/area" + num2str(sub_i) + "/non_ibr/infinite bus";
str_float = "/Subsystem" + num2str(sub_i) + "/area" + num2str(sub_i) + "/non_ibr/floating buses";
str_load = "/Subsystem" + num2str(sub_i) + "/area" + num2str(sub_i) + "/non_ibr/load buses";
str_gfl = "/Subsystem" + num2str(sub_i) + "/area" + num2str(sub_i) + "/ibr/gfl";
str_gfm = "/Subsystem" + num2str(sub_i) + "/area" + num2str(sub_i) + "/ibr/gfm";

%% Device 1: INFINITE BUS %%
if n_inf > 0
    set_param(file1+str_inf,'commented','off');
    % constants
    inf.v_D = real(v_bus_DQ(busnum.inf));
    inf.v_Q = imag(v_bus_DQ(busnum.inf));
    % initial conditions
    % multipliers
end
if n_inf == 0
    inf = 0;
    set_param(file1+str_inf,'commented','on');
end
%% Device 2: FLOATING BUS %%
if n_float > 0
    set_param(file1+str_float,'commented','off');
    % initial conditions
    float.v_D = real(v_bus_DQ(busnum.float));
    float.v_Q = imag(v_bus_DQ(busnum.float));
    % multipliers
    float.c = c_bus(busnum.float);
end
if n_float == 0
    float = 0;
    set_param(file1+str_float,'commented','on');
end
%% Device 3: LOAD BUS %%
if n_load > 0
    set_param(file1+str_load,'commented','off');
    % constants
    load.p = sys_data(busnum.load,5);
    load.q = sys_data(busnum.load,6);
    load.ptype = sys_data(busnum.load,22);
    load.qtype = sys_data(busnum.load,23);
    % initial conditions
    load.v_D = real(v_bus_DQ(busnum.load));
    load.v_Q = imag(v_bus_DQ(busnum.load));
    % multipliers
    load.c = c_bus(busnum.load);
end
if n_load == 0
    load = 0;
    set_param(file1+str_load,'commented','on');
end

%%
% This is how sys_data IBR values are tabulated, starting from 8th index:
% 8 |    9    |    10   |    11    |    12    |    13    |    14      |  15  |  16  |    17    | 18 | 19 |  20 | 21 |
% Bf|i-Loop-BW|P-meas-BW|QV-meas-BW|Dpw_or_Dwp|Dqv_or_Dvq|Phase_Margin|PLL-kp|PLL-ki|PLL-LPF-BW|P-kp|P-ki|QV-kp|QV-ki
%% Device 4: GFL-PV %%
if n_gfl > 0
    set_param(file1+str_gfl,'commented','off');

    % constants
    gfl.p_ref = p_ref(busnum.gfl);
    gfl.v_ref = v_ref(busnum.gfl);
    gfl.q_ref = q_ref(busnum.gfl);

    % initial conditions
    gfl.v_d = real(v_bus_dq(busnum.gfl));
    gfl.v_q = imag(v_bus_dq(busnum.gfl));
    gfl.i_d = real(i_inj_dq(busnum.gfl));
    gfl.i_q = imag(i_inj_dq(busnum.gfl));
    gfl.theta = theta_bus(busnum.gfl);

    % multipliers
    gfl.c = c_bus(busnum.gfl);
    gfl.i_loop_w = sys_data(busnum.gfl,9)*2*pi;
    gfl.p_meas_w = sys_data(busnum.gfl,10)*2*pi;
    gfl.v_meas_w = sys_data(busnum.gfl,11)*2*pi;
    gfl.dwp = sys_data(busnum.gfl,12); % zero means droop disabled
    gfl.dqv = sys_data(busnum.gfl,13); % zero means droop disabled
    gfl.pll_lpf_w = sys_data(busnum.gfl,17)*2*pi;
end
if n_gfl == 0
    gfl = 0;
   set_param(file1+str_gfl,'commented','on');
end

%% Device 5: GFM %%
if n_gfm > 0
    set_param(file1+str_gfm,'commented','off');

    % constants
    gfm.p_ref = p_ref(busnum.gfm);
    gfm.q_ref = q_ref(busnum.gfm);
    gfm.v_ref = v_ref(busnum.gfm);

    % initial conditions
    gfm.v_d = real(v_bus_dq(busnum.gfm));
    gfm.v_q = imag(v_bus_dq(busnum.gfm));
    gfm.i_d = real(i_inj_dq(busnum.gfm));
    gfm.i_q = imag(i_inj_dq(busnum.gfm));
    gfm.theta = theta_bus(busnum.gfm);

    % multipliers
    gfm.c = c_bus(busnum.gfm);
    gfm.i_loop_w = sys_data(busnum.gfm,9)*2*pi;
    gfm.p_meas_w = sys_data(busnum.gfm,10)*2*pi;
    gfm.q_meas_w = sys_data(busnum.gfm,11)*2*pi;

end
if n_gfm == 0
    gfm = 0;
   set_param(file1+str_gfm,'commented','on');
end
%% LINE %%
% constants
% initial conditions
line.i_D = real(i_line_DQ);
line.i_Q = imag(i_line_DQ);
% multipliers
line.r = r_line;
line.l = x_line/w;
end