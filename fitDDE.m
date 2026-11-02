% fit dde model M636

% data read in
data = readtable('blowfly_data.csv');
t_data = data.time;
N_data = data.total;
N0 = N_data(1);

% initial data plot
figure;
plot(t_data, N_data, 'k-', 'LineWidth', 1.3);
hold on;
plot(t_data, N_data, 'b.', 'MarkerSize', 15, 'MarkerFaceColor', 'b');
hold off;
grid on;
set(gcf, 'Color', 'w');
ax = gca;
set(ax, 'Color', 'w');
set(ax, 'XColor', 'k', 'YColor', 'k');
set(ax, 'Box', 'on');
xlabel('Days');
ylabel('Number of adult blowflies');
title('Nicholsons blowfly population', 'Color', 'k')
grid on;


%*************** Quasi Cycle Hypothesis ***************%
PTd = 23.5;
dTd = 3.0;
Td = 14.4;
Nd = 1500; %not given

params = [PTd, dTd, Td, Nd];

bounds.lower = [23.5-4.5, 3.0-0.7, 14.8-0.4,   1];
bounds.upper = [23.5+4.5, 3.0+0.7, 14.8+0.4, 1e5];

modelFun = @(p, t) ddeRun(p, t, N0);

opts = optimoptions('lsqcurvefit', 'Display', 'iter', ...
    'MaxFunctionEvaluations', 800, 'MaxIterations', 400, 'FiniteDifferenceType', 'central');

[pFit, err, ~, ~, ~, ~, J] = ...
    lsqcurvefit(modelFun, params, t_data, N_data, bounds.lower, bounds.upper, opts);

% SE
Npts = numel(N_data);
Npar = numel(pFit);
dof  = max(Npts - Npar, 1);
sigma = err / dof;
std_err = sqrt(diag(sigma*inv(full(J' * J))));

fprintf('PTd = %.4f (std err %.4f)\n', pFit(1), std_err(1));
fprintf('dTd = %.4f (std err %.4f)\n', pFit(2), std_err(2));
fprintf('Td = %.4f (std err %.4f)\n', pFit(3), std_err(3));
fprintf('Nd = %.4f (std err %.4f)\n', pFit(4), std_err(4));
fprintf('P = %.4f\n', pFit(1)/pFit(3));
fprintf('d = %.4f\n', pFit(2)/pFit(3));

% plot quasi cycle fit
fit_quasi = modelFun(pFit, t_data);

figure;
plot(t_data, N_data, '.', 'MarkerSize', 10, 'DisplayName', 'Data'); hold on;
plot(t_data, fit_quasi, '-', 'LineWidth', 1.8, 'DisplayName', 'Model fit');
set(gcf, 'Color', 'w');
ax = gca;
set(ax, 'Color', 'w');
set(ax, 'XColor', 'k', 'YColor', 'k');
set(ax, 'Box', 'on');
xlabel('Days');
ylabel('Number of Adult Blowflies');
legend('Location', 'best');
grid on;
title('Data Fit for Quasi Cycle Hypothesis', 'Color', 'k');


%*************** Limit Cycle Hypothesis ***************%
PTd = 150;
dTd = 2.9;
Td = 14.8;
Nd = 500; %not given

params = [PTd, dTd, Td, Nd];

bounds.lower = [150-70, 2.9-0.5, 14.8-0.4,   1];
bounds.upper = [150+70, 2.9+0.5, 14.8+0.4, 1e5];

modelFun = @(p, t) ddeRun(p, t, N0);

opts = optimoptions('lsqcurvefit', 'Display', 'iter', ...
    'MaxFunctionEvaluations', 800, 'MaxIterations', 400, 'FiniteDifferenceType', 'central');

[pFit, err, ~, ~, ~, ~, J] = ...
    lsqcurvefit(modelFun, params, t_data, N_data, bounds.lower, bounds.upper, opts);

% SE
Npts = numel(N_data);
Npar = numel(pFit);
dof  = max(Npts - Npar, 1);
sigma = err / dof;
std_err = sqrt(diag(sigma*inv(full(J' * J))));

fprintf('PTd = %.4f (std err %.4f)\n', pFit(1), std_err(1));
fprintf('dTd = %.4f (std err %.4f)\n', pFit(2), std_err(2));
fprintf('Td = %.4f (std err %.4f)\n', pFit(3), std_err(3));
fprintf('Nd = %.4f (std err %.4f)\n', pFit(4), std_err(4));
fprintf('P = %.4f\n', pFit(1)/pFit(3));
fprintf('d = %.4f\n', pFit(2)/pFit(3));

% plot fit
fit_limit = modelFun(pFit, t_data);

figure;
plot(t_data, N_data, '.', 'MarkerSize', 10, 'DisplayName', 'Data'); hold on;
plot(t_data, fit_limit, '-', 'LineWidth', 1.8, 'DisplayName', 'Model fit');
set(gcf, 'Color', 'w');
ax = gca;
set(ax, 'Color', 'w');
set(ax, 'XColor', 'k', 'YColor', 'k');
set(ax, 'Box', 'on');
xlabel('Days');
ylabel('Number of Adult Blowflies');
legend('Location', 'best');
grid on;
title('Data Fit for Limit Cycle Hypothesis', 'Color', 'k');


%dde function
function modRun = ddeRun(p, time, N0)

    PTd = p(1);
    dTd = p(2);
    Td = p(3);
    Nd = p(4);
    P = PTd/Td;
    delta = dTd/Td;

    % DDE equation
    dde = @(t, N, Z) P.*Z(1).*exp(-Z(1)./Nd)-delta.*N;

    history = @(t) N0;

    tAll = [0, max(time)];
    sol = dde23(dde, Td, history, tAll);

    %eval
    modRun = deval(sol, time);
    modRun = modRun(1, :).';
end
