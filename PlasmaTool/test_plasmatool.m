% test_plasmatool.m
% Basic sanity checks for plasma.computeAll

clear; clc;

base = struct('n_val', 1e19, 'n_unit', 'm^-3', 'T_val_eV', 1e3, ...
    'Z', 1, 'MW_amu', 2.0141, 'B_T', 2.0, 'lnLambda_mode', 'auto', 'lnLambda', 15);

% 1) Debye length decreases with density
p1 = base; p1.n_val = 1e18;
p2 = base; p2.n_val = 1e20;
r1 = plasma.computeAll(p1);
r2 = plasma.computeAll(p2);
assert(r2.lambda_D_m < r1.lambda_D_m, 'Debye length should decrease with density.');

% 2) Electron gyrofrequency scales linearly with B
p3 = base; p3.B_T = 1;
p4 = base; p4.B_T = 3;
r3 = plasma.computeAll(p3);
r4 = plasma.computeAll(p4);
assert(abs(r4.Omega_ce_rad_s/r3.Omega_ce_rad_s - 3) < 1e-12, 'Omega_ce should scale with B.');

% 3) nu_ei ~ T^(-3/2) (for fixed lnLambda-ish)
p5 = base; p5.T_val_eV = [1e3; 8e3]; p5.lnLambda_mode = 'manual'; p5.lnLambda = 15;
r5 = plasma.computeAll(p5);
ratio = r5.nu_ei_Hz(1) / r5.nu_ei_Hz(2);
assert(ratio > 10 && ratio < 30, 'nu_ei should decrease strongly with T^(3/2).');

% 4) Pressure linear in T
p6 = base; p6.T_val_eV = [1e3; 2e3];
r6 = plasma.computeAll(p6);
assert(abs(r6.p_Pa(2)/r6.p_Pa(1) - 2) < 1e-12, 'Pressure should scale linearly with T.');

% 5) eta = 1/sigma consistency
assert(max(abs(r6.eta_par_Ohm_m .* r6.sigma_par_S_m - 1)) < 1e-12, 'eta*sigma should equal 1.');

fprintf('All PlasmaTool sanity checks passed.\n');
