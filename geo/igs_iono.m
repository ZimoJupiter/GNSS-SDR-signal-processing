function I_slant = igs_iono( time_user, user_ll, el)

% 0223_middle_1
% early_4grids = [680 609 553 629];
% late_4grids  = [479 396 371 457];
% time_early = 374415;
% time_late  = 374415 + 2 * 3600; 

% small_static_3
early_4grids = [129 112 110 123];
late_4grids  = [114 108 107 108];
time_early = 396015;
time_late  = 396015 + 2 * 3600; 


e01 = [ 25.0 120 ]; e11 = [ 25.0 125 ];
e00 = [ 22.5 120 ]; e01 = [ 22.5 125 ];

delta_alpha = 2.5;
delta_phi   = 5.0;

e00_gt = [ early_4grids(1) late_4grids(1) ];
e01_gt = [ early_4grids(2) late_4grids(2) ];
e11_gt = [ early_4grids(3) late_4grids(3) ];
e10_gt = [ early_4grids(4) late_4grids(4) ];

e00_ut = ( (time_late-time_user) / (60*60*2) ) * e00_gt(1) + ( (time_user-time_early) / (60*60*2) ) * e00_gt(2); % time unit: sec
e01_ut = ( (time_late-time_user) / (60*60*2) ) * e01_gt(1) + ( (time_user-time_early) / (60*60*2) ) * e01_gt(2); % time unit: sec
e11_ut = ( (time_late-time_user) / (60*60*2) ) * e11_gt(1) + ( (time_user-time_early) / (60*60*2) ) * e11_gt(2); % time unit: sec
e10_ut = ( (time_late-time_user) / (60*60*2) ) * e10_gt(1) + ( (time_user-time_early) / (60*60*2) ) * e10_gt(2); % time unit: sec

user_lat = user_ll(1);
user_lon = user_ll(2);

p = user_lat - e00(1) / delta_alpha;
q = user_lon - e00(2) / delta_alpha;

% e_user = (1-p) * (1-q) * e00_ut + p * (1-q) * e10_ut + q * (1-p) * e01_ut + p * q * e11_ut
e_user = cos((1-p) * (1-q) * e00_ut + p * (1-q) * e10_ut + q * (1-p) * e01_ut + p * q * e11_ut);

Re = 6378.137; % radius of the earth
h = 450; % height of the ionosphere layer
fL1 = 1575.42e6;
OF = sec( asin( ( Re/(Re+h) * cos(el) ) ) );
I_vertical = 40.3/(fL1^2) * e_user * 1e16 * 1-1;
I_slant = OF * I_vertical;

end