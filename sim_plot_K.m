clear;
naz = 16;
nel = 76;
n_arr = naz * nel;
K = 100;
Rmin = 10;
Rmax = 100;
sigman2 = 10^(-7.2)*1e-3;
K_store = [10; 40; 70; 100;];%
rho = 1e1;%rho/N in the paper
height = 10;
f = 80*10^9;
lambda = 3 * 10^8 / f;
miu = 0.5;
betam_min = atan(height/Rmax);
betam_max = atan(height/Rmin);
mechanical_tilt_opt = 0;
beta_m = [0.5*(betam_min+betam_max); mechanical_tilt_opt;];%-phi_m in the paper
elevation_range = [-0.4703, -0.7424; 0.3429, -0.0997];
N_inex = 1e2;%for electrical tilt
Nite = 1e2;
capacity = zeros(length(K_store), length(beta_m) * 2);
angle = zeros(naz * nel, 2);
for nx = 1 : naz
    for ny = 1 : nel
        n = (ny - 1) * naz + nx;
        angle(n, 2) = (-1+2*ny/nel);%el
        angle(n, 1) = (-1+2*nx/naz);%az
    end
end
for signalpo_n = 1 : length(K_store)
    K = K_store(signalpo_n, 1);
    Mb = 4*K;
    beam_numberforMS = floor(Mb/K);
    c_km = zeros(K, beam_numberforMS);
    el_in = zeros(Mb, 1);
    az_in = zeros(Mb, 1);
    for ii = 1 : Nite
        pos = zeros(K, 3);
        theta = zeros(K, length(beta_m));
        phi = zeros(K, length(beta_m));
        beta = zeros(K, length(beta_m) * 2);
        beam_in = zeros(Mb, length(beta_m));
        H = zeros(n_arr, K * length(beta_m));
        for k = 1 : K
            pos_temp = zeros(1, 2);
            while norm(pos_temp) < Rmin || norm(pos_temp) > Rmax || abs(atan(pos_temp(1, 2) / pos_temp(1, 1))) > pi / 3
                pos_temp(1, 1) = rand(1, 1) * Rmax;
                pos_temp(1, 2) = (rand(1, 1) * 2 - 1) * Rmax;
            end
            pos(k, 1:2) = pos_temp;
            pos(k, 3) = norm(pos_temp);
            pos(k, 3) = norm([pos(k, 3), height]);%distance
            for capa_n = 1 : length(beta_m)
                phi(k, capa_n) = asin(pos_temp(1, 2) / sqrt(pos_temp(1, 2)^2 + (pos_temp(1, 1) * cos(beta_m(capa_n, 1)) + height * sin(beta_m(capa_n, 1)))^2));%az
                theta(k, capa_n) = asin((pos_temp(1, 1) * sin(beta_m(capa_n, 1)) - height * cos(beta_m(capa_n, 1))) / pos(k, 3));%el
                Aaz = -min(12 * phi(k, capa_n)^2 / (70/180*pi)^2, 25);
                Ael = -min(12 * theta(k, capa_n)^2 / (7/180*pi)^2, 20);
                D0 = -min(-Aaz-Ael, 25);
                D0 = 10^(D0*0.1);
                beta(k, capa_n) = sqrt(D0 * lambda^2 / (16 * pi^2 * pos(k, 3)^2)) * exp(1i * rand(1,1) * 2 * pi);
                h = zeros(n_arr, 1);
                for nx = 0 : naz-1
                    for ny = 0 : nel-1
                        n = ny * naz + 1 + nx;
                        h(n, 1) = exp(-1i * 2 * pi * miu * ((nx-0.5*(naz-1)) * sin(phi(k, capa_n)) * cos(theta(k, capa_n)) + (ny-0.5*(nel-1)) * sin(theta(k, capa_n))));
                    end
                end
                H(:, K * (capa_n - 1) + k) = h * beta(k, capa_n);
                diff = 1e4 * ones(n_arr, 1);
                for n = 1 : n_arr
                    flag = 0;
                    for m = 1 : beam_numberforMS * (k-1)
                        if abs(n-beam_in(m, capa_n)) < 1e-2
                            flag = 1;
                        else
                        end
                    end
                    if abs(flag-1) < 1e-2
                    else
                        diff(n, 1) = norm([sin(phi(k, capa_n)) * cos(theta(k, capa_n))-angle(n, 1), sin(theta(k, capa_n))-angle(n, 2)]);
                    end
                end
                [~, index] = sort(diff);
                beam_in(beam_numberforMS*(k-1)+1:beam_numberforMS*k, capa_n) = index(1:beam_numberforMS, 1);
            end
        end
        for capa_n = 1 : length(beta_m)
            U = zeros(n_arr, Mb);
            for nk = 1 : Mb
                if abs(beam_in(nk, capa_n)/naz - floor(beam_in(nk, capa_n)/naz))<1e-2
                    x_temp = naz;
                    y_temp = floor(beam_in(nk, capa_n)/naz);
                else
                    y_temp = floor(beam_in(nk, capa_n)/naz) + 1;
                    x_temp = beam_in(nk, capa_n) - (y_temp-1) * naz;
                end
                el_in(nk, 1) = -1+2*y_temp/nel;
                az_in(nk, 1) = -1+2*x_temp/naz;
                for nx = 0 : naz-1
                    for ny = 0 : nel-1
                        n = ny * naz + nx + 1;
                        U(n, nk) = exp(-1i * 2 * pi * miu * ((nx-0.5*(naz-1)) * az_in(nk, 1) + (ny-0.5*(nel-1)) * el_in(nk, 1))) / sqrt(n_arr);
                    end
                end
            end
                Hb = U' * H(:, K * (capa_n - 1) + 1 : K * capa_n);
                F = Hb / (Hb' * Hb + K * sigman2 / rho * eye(K));
                alpha = sqrt(K * rho / trace(F * F'));
                A = Hb' * F;
                for k = 1 : K
                    interf = 0;
                    for kk = 1 : K
                        if kk == k
                        else
                            interf = interf + abs(A(k, kk))^2;
                        end
                    end
                    interf = alpha^2 * interf / K + sigman2;
                    sinr = alpha^2 * abs(A(k, k))^2 / K / interf;
                    capacity(signalpo_n, capa_n) = capacity(signalpo_n, capa_n) + log2(1 + sinr);
                end
            for k = 1 : K
                for m = 1 : beam_numberforMS
                    dkmaz = sin(phi(k, capa_n)) * cos(theta(k, capa_n)) - az_in((k-1)*beam_numberforMS + m, 1);
                    dkmel = sin(theta(k, capa_n)) - el_in((k-1)*beam_numberforMS + m, 1);
                    c_km(k, (k-1)*beam_numberforMS + m) = abs(1 - exp(1i * pi * naz * dkmaz)) * abs(1 - exp(1i * pi * nel * dkmel)) / (pi^2 * n_arr * abs(dkmaz * dkmel));
                end
            end
            capacity_temp = zeros(N_inex, 1);
            for n = 1 : N_inex
                betae_n = elevation_range(1, capa_n) + (elevation_range(2, capa_n) - elevation_range(1, capa_n)) / (N_inex-1) * (n - 1);
                ctemp = zeros(K, 1);
                ano_temp = 0;
                for k = 1 : K
                    Aaz = -min(12 * phi(k, capa_n)^2 / (70/180*pi)^2, 25);
                    Ael = -min(12 * (theta(k, capa_n)-betae_n)^2 / (7/180*pi)^2, 20);
                    D0 = -min(-Aaz-Ael, 25);
                    D0 = 10^(D0*0.1);
                    beta(k, length(beta_m) + capa_n) = sqrt(D0 * lambda^2 / (16 * pi^2 * pos(k, 3)^2)) * beta(k, capa_n) / abs(beta(k, capa_n));
                    ctemp(k, 1) = abs(beta(k, length(beta_m) + capa_n))^2 * sum((c_km(k, (k-1)*beam_numberforMS + 1 : k*beam_numberforMS)).^2);
                    ano_temp = ano_temp + ctemp(k, 1) / (ctemp(k, 1) + sigman2 * K / rho / n_arr)^2;
                end
                for k = 1 : K
                    capacity_temp(n, 1) = capacity_temp(n, 1) + log2(1 + rho * n_arr / ano_temp * (ctemp(k, 1) / (ctemp(k, 1) + sigman2 * K / rho / n_arr))^2 / sigman2) ;
                end
            end
            [~, capa_index] = max(capacity_temp(:, 1));
            betae = elevation_range(1, capa_n) + (elevation_range(2, capa_n) - elevation_range(1, capa_n)) / (N_inex-1) * (capa_index - 1);
            for k = 1 : K
                Aaz = -min(12 * phi(k, capa_n)^2 / (70/180*pi)^2, 25);
                Ael = -min(12 * (theta(k, capa_n)-betae)^2 / (7/180*pi)^2, 20);
                D0 = -min(-Aaz-Ael, 25);
                D0 = 10^(D0*0.1);
                beta(k, length(beta_m) + capa_n) = sqrt(D0 * lambda^2 / (16 * pi^2 * pos(k, 3)^2)) * beta(k, capa_n) / abs(beta(k, capa_n));
                Hb(:, k) = U' * H(:, K * (capa_n - 1) + k) / beta(k, capa_n) * beta(k, length(beta_m) + capa_n);
            end
            F = Hb / (Hb' * Hb + K * sigman2 / rho * eye(K));
            alpha = sqrt(K * rho / trace(F * F'));
            A = Hb' * F;
            for k = 1 : K
                interf = 0;
                for kk = 1 : K
                    if kk == k
                    else
                        interf = interf + abs(A(k, kk))^2;
                    end
                end
                interf = alpha^2 * interf / K + sigman2;
                sinr = alpha^2 * abs(A(k, k))^2 / K / interf;
                capacity(signalpo_n, length(beta_m) + capa_n) = capacity(signalpo_n, length(beta_m) + capa_n) + log2(1 + sinr);
            end
        end
        disp([signalpo_n, ii])
    end
end
capacity = capacity / Nite;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure;
set(h,'PaperType','A4');
xx = axes('FontSize',16);
semilogy(K_store, capacity(1:length(K_store), 1), 'k--d','LineWidth',2,'MarkerSize',14)
hold on
semilogy(K_store, capacity(1:length(K_store), 3), 'k--^','LineWidth',2,'MarkerSize',14)
semilogy(K_store, capacity(1:length(K_store), 4), 'k-s','LineWidth',2,'MarkerSize',14)

xlim([min(K_store), max(K_store)])
le = legend('fixed mechanical tilt','dynamic electrical tilt', 'the proposed tilts', 'Location', 'northwest');
set(le,'Fontsize',16,'Fontname','Times')
set(gca,'XTick',K_store)
xlabel('Number of MSs K','Fontsize',20,'Fontname','Times')
ylabel('Spectrum efficiency (bps/Hz)','Fontsize',20,'Fontname','Times')
grid on
print(h,'-dpdf','capacity_K')