# PlasmaTool (MATLAB R2022b+)

A lightweight App Designer-style GUI and modular physics backend to compute basic plasma parameters versus temperature.

## File tree

```text
PlasmaTool/
  PlasmaTool.mlapp
  test_plasmatool.m
  README.md
  +plasma/
    constants.m
    units.m
    computeAll.m
    plotSweep.m
```

## Data model

### `params` struct

```matlab
params = struct( ...
  'n_val',        1e19,        ... % numeric density value
  'n_unit',       'm^-3',      ... % ''m^-3'' | ''cm^-3''
  'T_val_eV',     1e3,         ... % scalar or vector [eV]
  'Z',            1,           ... % mean charge state
  'MW_amu',       2.0141,      ... % ion molecular weight [amu]
  'B_T',          2.0,         ... % magnetic field [T]
  'lnLambda_mode','auto',      ... % ''auto'' | ''manual''
  'lnLambda',     15);             % used if manual
```

### `results` struct (selected fields)

```matlab
results.T_eV
results.p_Pa
results.v_te_ms
results.v_ti_ms
results.omega_pe_rad_s
results.omega_pi_rad_s
results.lambda_D_m
results.Omega_ce_rad_s
results.Omega_ci_rad_s
results.nu_ei_Hz
results.nu_ii_Hz
results.beta_e
results.beta_i
results.sigma_par_S_m
results.eta_par_Ohm_m
results.kappa_e_par_W_mK
results.kappa_i_par_W_mK
```

## Usage

1. Add `PlasmaTool/` to the MATLAB path.
2. Launch the app:
   ```matlab
   app = PlasmaTool;
   ```
3. Set inputs:
   - Density `n` with units (`m^-3` or `cm^-3`)
   - Temperature mode:
     - **single**: one temperature (`eV` or `keV`)
     - **sweep**: log-spaced 1 eV to 100 keV with configurable point count
   - Mean charge state `Z`
   - Ion species (`H, D, T, He, Ar`) and/or MW override
   - Magnetic field `B` in Tesla
   - Coulomb logarithm `lnΛ` mode (`auto` or `manual`)
4. Click **Compute** for scalar outputs or **Plot Sweep** for tabs/plots.
5. Use **Export CSV** to save currently computed/swept data.

## Physics model and formulas (SI internal)

Let `T_eV` be electron temperature in eV and `n` in `m^-3`.

- `k_B T = e * T_eV`
- Pressure:
  - `p = n * k_B T = n * e * T_eV`
- Thermal speeds:
  - `v_te = sqrt(2 e T_eV / m_e)`
  - `v_ti = sqrt(2 e T_i,eV / m_i)`
- Plasma frequencies:
  - `ω_pe = sqrt(n e^2 / (ε0 m_e))`
  - `ω_pi = sqrt(n Z^2 e^2 / (ε0 m_i))`
- Debye length:
  - `λ_D = sqrt(ε0 k_B T / (n e^2))`
- Gyrofrequencies:
  - `Ω_ce = e B / m_e`
  - `Ω_ci = Z e B / m_i`
- Collision frequencies (Spitzer practical form, with `n` in `cm^-3`, `T` in `eV`):
  - `ν_ei ≈ 2.91×10^-6 n_e Z lnΛ / T_e^(3/2)`
  - `ν_ii ≈ 4.80×10^-8 Z^4 n_i lnΛ / (sqrt(A) T_i^(3/2))`
- Magnetization/Hall parameters:
  - `β_e = Ω_ce / ν_ei`
  - `β_i = Ω_ci / ν_ii`
- Parallel transport (Lorentz/Spitzer-Braginskii style approximations):
  - `σ_|| = n_e e^2 / (m_e ν_ei)`
  - `η_|| = 1 / σ_||`
  - `κ_e,|| ≈ 3.16 n_e k_B^2 T / (m_e ν_ei)`
  - `κ_i,|| ≈ 3.90 n_i k_B^2 T / (m_i ν_ii)`
- Bonus cross-field/Hall reductions:
  - `σ_⊥ = σ_|| / (1 + β_e^2)`
  - `σ_H = σ_|| β_e / (1 + β_e^2)`

### Constants and sources

- CODATA physical constants (exact SI redefinition values where applicable)
- Spitzer-like collision frequency fits (commonly used in NRL Plasma Formulary-style engineering estimates)
- Braginskii transport scaling coefficients (`3.16`, `3.9`) for parallel thermal conductivity approximations

## Notes and assumptions

- Quasi-neutral single-ion plasma with `n_i = n_e / Z` not explicitly separated in this first-pass engineering model.
- `T_i = T_e` assumed unless extended.
- `lnΛ` auto mode uses a bounded practical estimate (`2..25`) suitable for broad scans.
- Outputs are vectorized for sweeps.

## Tests

Run:

```matlab
cd PlasmaTool
test_plasmatool
```

Included sanity checks:
1. `λ_D` decreases with `n`.
2. `Ω_ce` scales linearly with `B`.
3. `ν_ei` decreases with `T^(3/2)`.
4. Pressure scales linearly with `T`.
5. `ησ = 1` consistency.
