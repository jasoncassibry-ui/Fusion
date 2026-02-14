function plotSweep(axStruct, results)
%PLASMA.PLOTSWEEP Plot sweep results to provided axes handles.
% axStruct fields: Scalars, Collisions, Transport, Magnetization

T = results.T_eV;

% Scalars
axes(axStruct.Scalars); cla(axStruct.Scalars);
loglog(T, results.lambda_D_m, 'LineWidth', 1.5); hold on;
loglog(T, results.v_te_ms, 'LineWidth', 1.5);
loglog(T, results.v_ti_ms, 'LineWidth', 1.5); hold off;
grid on; xlabel('T_e [eV]');
ylabel('Value (SI)');
legend('\lambda_D [m]', 'v_{te} [m/s]', 'v_{ti} [m/s]', 'Location', 'best');
title('Scalar quantities');

% Collisions
axes(axStruct.Collisions); cla(axStruct.Collisions);
loglog(T, results.nu_ei_Hz, 'LineWidth', 1.5); hold on;
loglog(T, results.nu_ii_Hz, 'LineWidth', 1.5); hold off;
grid on; xlabel('T_e [eV]'); ylabel('\nu [s^{-1}]');
legend('\nu_{ei}', '\nu_{ii}', 'Location', 'best');
title('Collision frequencies');

% Transport
axes(axStruct.Transport); cla(axStruct.Transport);
loglog(T, results.sigma_par_S_m, 'LineWidth', 1.5); hold on;
loglog(T, results.eta_par_Ohm_m, 'LineWidth', 1.5);
loglog(T, results.kappa_e_par_W_mK, 'LineWidth', 1.5);
loglog(T, results.kappa_i_par_W_mK, 'LineWidth', 1.5); hold off;
grid on; xlabel('T_e [eV]'); ylabel('Transport coeff. (SI)');
legend('\sigma_{||}', '\eta_{||}', '\kappa_{e,||}', '\kappa_{i,||}', 'Location', 'best');
title('Parallel transport');

% Magnetization
axes(axStruct.Magnetization); cla(axStruct.Magnetization);
loglog(T, results.beta_e, 'LineWidth', 1.5); hold on;
loglog(T, results.beta_i, 'LineWidth', 1.5);
loglog(T, abs(results.sigma_hall_S_m), '--', 'LineWidth', 1.2); hold off;
grid on; xlabel('T_e [eV]'); ylabel('Dimensionless / SI');
legend('\beta_e', '\beta_i', '|\sigma_H| [S/m]', 'Location', 'best');
title('Magnetization / Hall response');
end
