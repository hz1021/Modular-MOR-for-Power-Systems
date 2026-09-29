function gfm_gains = design_gains_gfm_simple(sys_data, gfm,ft)

gfm_gains = struct();


tau = 1./(2*pi*ft);

Kp = gfm.c./tau;
Ki = 0*ones(length(gfm.p_ref),1);

gfm_gains.vd_kp = Kp;
gfm_gains.vd_ki = Ki;

gfm_gains.vq_kp = Kp;
gfm_gains.vq_ki = Ki;

temp_bus = (1:size(sys_data,1))';
device_data = sys_data(:,[1,2,7]);
busnum.gfm = temp_bus(device_data(:,3)==200,1); % getting the gfm number

gfm_gains.dpw = sys_data(busnum.gfm,12); % it must be nonzero for all gfms
gfm_gains.dqv = sys_data(busnum.gfm,13); % it must be nonzero for gfms operating in pq mode