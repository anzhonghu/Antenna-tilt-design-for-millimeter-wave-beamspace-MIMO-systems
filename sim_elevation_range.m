clear;
naz = 16;
nel = 76;
n_arr = naz * nel;
K = 100;
Mb = 4*K;
Rmin = 10;
Rmax = 100;
sigman2 = 10^(-7.2)*1e-3;
sigmas_store = [1e0; 1e2;  1e4; 1e6];
height = 10;
f = 80*10^9;
lambda = 3 * 10^8 / f;
miu = 0.5;
Nite = 1e2;
capacity = zeros(length(sigmas_store), 2);
angle = zeros(naz * nel, 2);
betam_min = atan(height/Rmax);
betam_max = atan(height/Rmin);
mechanical_tilt_opt = 0;
betam_st = [0.5*(betam_min+betam_max); mechanical_tilt_opt];
angle_store = zeros(1e4, 2);
diff = 1e4 * ones(n_arr, 1);
N_inex = 100;
index_store = zeros(n_arr, N_inex);
for n = 1 :2
    betam_n = betam_st(n, 1);
    for n1 = 1 : 100
        r = Rmin + (Rmax - Rmin) / 100 * n1;
        for n2 = 1 : 100
            eta = (-60 + 120 / 100 * n2) / 180 * pi;
            nn = (n1-1)*100+n2;
            angle_store(nn, n) = asin((r * cos(eta) * sin(betam_n) - height * cos(betam_n)) / sqrt(r^2 + height^2));%el
        end
        disp([n, n1])
    end
end
angle_store = sort(angle_store);