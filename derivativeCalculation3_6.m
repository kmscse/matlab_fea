% MATLAB codes for Finite Element Analysis
% Derivatice Calculation
% Author: Kaung Myat San
% Location: England, the United Kingdom
% Time: 07:10 P.M London Time 22/01/2024

x0 = 1.0;   % the point where the derivative is calculated
dx = 10.^[-3:.1:0];     % delta x values varying from 1 to 10^-3
for k=1:length(dx)
    % forward difference scheme
    dudx_f(k) = ((x0+dx(k))^3-(x0+dx(k))^2 - x0^3 + x0^2)/dx(k);
    % backward difference scheme
    dudx_c(k) = ((x0+dx(k))^3 - (x0+dx(k))^2 ...
        -(x0-dx(k))^3+(x0-dx(k))^2)/(2*dx(k));
end

% plot the results
figure(1);
semilogx(dx, dudx_f, 'k--', 'LineWidth',2);
hold on
semilogx(dx, dudx_b,'k:','LineWidth',2);
semilogx(dx, dudx_c,'k:','LineWidth',2);
set(gca,' fontsize ', 16);
xlabel('x', ' fontsize ', 18);
ylabel('du/dx', ' fontsize ', 18);
legend('Forward difference', 'Backward difference', ...
    'Central difference', 'Location','northwest');