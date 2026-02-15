function results = computeAll(params)
%PLASMA.COMPUTEALL Compute core plasma parameters using SI internally.
% Inputs (params fields):
%   n_val, n_unit ('m^-3'|'cm^-3'), T_val_eV (scalar/vector eV), Z, MW_amu,
%   B_T, lnLambda_mode ('auto'|'manual'), lnLambda (if manual)
%
% Formulas (Spitzer/Braginskii baseline):
%   p = n*kB*T_K = n*T_eV*e
%   v_te = sqrt(2*kB*T/me) = sqrt(2*e*T_eV/me)
%   v_ti = sqrt(2*kB*T/mi)
%   omega_pe = sqrt(n*e^2/(eps0*me))
%   omega_pi = sqrt(n*Z^2*e^2/(eps0*mi))
%   lambda_D = sqrt(eps0*kB*T/(n*e^2))
%   Omega_ce = e*B/me, Omega_ci = Z*e*B/mi
%   nu_ei ≈ 2.91e-6 * n_e[cm^-3] * Z * lnL / Te[eV]^(3/2)
%   nu_ii ≈ 4.80e-8 * Z^4 * n_i[cm^-3] * lnL / (sqrt(A)*Ti[eV]^(3/2))
%   sigma_par = n_e*e^2/(me*nu_ei), eta = 1/sigma
%   kappa_e_par ~= 3.16 * n_e*kB^2*T_K/(me*nu_ei)
%   kappa_i_par ~= 3.9  * n_i*kB^2*T_K/(mi*nu_ii)

c = plasma.constants();
validateParams(params);

n_m3 = params.n_val;
if strcmpi(params.n_unit, 'cm^-3')
    n_m3 = n_m3 * 1e6;
end

Te_eV = params.T_val_eV(:);
Ti_eV = Te_eV; % Assume thermal equilibrium if not separately supplied.
T_K = Te_eV * c.e / c.kB;

Z = params.Z;
A = params.MW_amu;
mi = A * c.amu;
B = params.B_T;

if isfield(params, 'lnLambda_mode') && strcmpi(params.lnLambda_mode, 'manual')
    lnL = params.lnLambda * ones(size(Te_eV));
else
    % Simple practical approximation for fusion/lab plasma (NRL-form style).
    n_cm3 = n_m3 / 1e6;
    lnL = 23 - log(sqrt(max(n_cm3, eps))*Z ./ max(Te_eV, eps).^(3/2));
    lnL = max(2, min(25, lnL));
end

% Core scalar/vector outputs.
p_Pa = n_m3 .* Te_eV .* c.e;
v_te = sqrt(2 .* c.e .* Te_eV ./ c.me);
v_ti = sqrt(2 .* c.e .* Ti_eV ./ mi);
omega_pe = sqrt(n_m3 * c.e^2 / (c.eps0 * c.me)) .* ones(size(Te_eV));
omega_pi = sqrt(n_m3 * (Z*c.e)^2 / (c.eps0 * mi)) .* ones(size(Te_eV));
lambda_D = sqrt(c.eps0 .* c.e .* Te_eV ./ (n_m3 * c.e^2));

Omega_ce = (c.e * B / c.me) .* ones(size(Te_eV));
Omega_ci = (Z * c.e * B / mi) .* ones(size(Te_eV));

n_cm3 = n_m3 / 1e6;
nu_ei = 2.91e-6 .* n_cm3 .* Z .* lnL ./ max(Te_eV, eps).^(3/2);
nu_ii = 4.80e-8 .* (Z^4) .* n_cm3 .* lnL ./ (sqrt(A) .* max(Ti_eV, eps).^(3/2));
nu_ie = nu_ii; % same-order approximation for reporting convenience

beta_e = Omega_ce ./ nu_ei;
beta_i = Omega_ci ./ nu_ii;

sigma_par = n_m3 * c.e^2 ./ (c.me .* nu_ei);
eta_par = 1 ./ sigma_par;

kappa_e_par = 3.16 .* n_m3 .* c.kB.^2 .* T_K ./ (c.me .* nu_ei);
kappa_i_par = 3.90 .* n_m3 .* c.kB.^2 .* T_K ./ (mi .* nu_ii);

% Optional simple magnetized tensor reductions (Lorentz model style)
sigma_perp = sigma_par ./ (1 + beta_e.^2);
sigma_hall = sigma_par .* beta_e ./ (1 + beta_e.^2);
kappa_e_perp = kappa_e_par ./ (1 + beta_e.^2);
kappa_i_perp = kappa_i_par ./ (1 + beta_i.^2);

results = struct();
results.inputs = params;
results.n_m3 = n_m3;
results.T_eV = Te_eV;
results.T_K = T_K;
results.mi_kg = mi;
results.lnLambda = lnL;

results.p_Pa = p_Pa;
results.v_te_ms = v_te;
results.v_ti_ms = v_ti;
results.omega_pe_rad_s = omega_pe;
results.omega_pi_rad_s = omega_pi;
results.lambda_D_m = lambda_D;
results.Omega_ce_rad_s = Omega_ce;
results.Omega_ci_rad_s = Omega_ci;

results.nu_ei_Hz = nu_ei;
results.nu_ii_Hz = nu_ii;
results.nu_ie_Hz = nu_ie;
results.beta_e = beta_e;
results.beta_i = beta_i;

results.sigma_par_S_m = sigma_par;
results.eta_par_Ohm_m = eta_par;
results.kappa_e_par_W_mK = kappa_e_par;
results.kappa_i_par_W_mK = kappa_i_par;
results.sigma_perp_S_m = sigma_perp;
results.sigma_hall_S_m = sigma_hall;
results.kappa_e_perp_W_mK = kappa_e_perp;
results.kappa_i_perp_W_mK = kappa_i_perp;
end

function validateParams(params)
required = {'n_val','n_unit','T_val_eV','Z','MW_amu','B_T'};
for k = 1:numel(required)
    if ~isfield(params, required{k})
        error('Missing required field: %s', required{k});
    end
end
if any(params.T_val_eV <= 0)
    error('Temperature must be > 0 eV.');
end
if params.n_val <= 0
    error('Density must be > 0.');
end
if params.Z <= 0
    error('Mean charge state Z must be > 0.');
end
if params.MW_amu <= 0
    error('Molecular weight (amu) must be > 0.');
end
if params.B_T < 0
    error('Magnetic field B must be >= 0 T.');
end
if isfield(params,'lnLambda_mode') && strcmpi(params.lnLambda_mode,'manual') && params.lnLambda <= 0
    error('Manual lnLambda must be > 0.');
end
end
