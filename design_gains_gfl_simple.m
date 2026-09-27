function gfl_pv_gains = design_gains_gfl_simple(sys_data,fpll, zeta, fp, fv, xg)

gfl_pv_gains = struct();



%%%%%%%%%%%%%% this is to know the number of each components %%%%%%%%%%%%

device_data = sys_data(:,[1,2,7]);
temp_bus = (1:size(sys_data,1))';
busnum.gfl = temp_bus(device_data(:,3)==100,1); % getting the gfl pv number

% from Adria's notes
wbw = 2*pi*fpll;
wn = wbw.*sqrt(2./(sqrt((4*zeta.^2+2).^2+4)+4.*zeta.^2+2));
gfl_pv_gains.kppll = -2.*zeta.*wn;
gfl_pv_gains.kipll = -wn.^2;

% from Adria's notes we only need current loop and P loop bw to get gains
fc = sys_data(busnum.gfl,9);
gfl_pv_gains.kpp = fp./fc;
gfl_pv_gains.kip = 2*pi.*fp;

gfl_pv_gains.kpv = -fv./(fc.*xg);
gfl_pv_gains.kiv = -2*pi*fp./xg;