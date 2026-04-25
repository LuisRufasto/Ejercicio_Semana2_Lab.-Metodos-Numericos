% =========================================================================
% Interpolación: Métodos de Lagrange y Newton (Versión Lineal)
% =========================================================================
clear; clc; close all;

% 1. DATOS 
x = [-0.9; -0.2; 0.3; 0.7; 1.5; 2.2; 3.2; 3.9; 4.5; 5.1; 5.7; 6.2];
y = [1.841; 4.872; 4.757; 3.873; 1.625; 1.128; 7.048; 18.689; 35.375; 59.621; 92.723; 128.008];
n = length(x);

% 2. VECTORES PARA LA GRÁFICA
xx = linspace(min(x), max(x), 500);
yy_lagrange = zeros(size(xx));
yy_newton = zeros(size(xx));

% =========================================================================
% 3. CÁLCULO DE INTERPOLACIÓN PARA LAS GRÁFICAS
% =========================================================================

% --- Preparar Coeficientes de Newton (Diferencias Divididas) ---
tabla = zeros(n, n);
tabla(:,1) = y;
for j = 2:n
    for i = 1:(n-j+1)
        tabla(i,j) = (tabla(i+1,j-1) - tabla(i,j-1)) / (x(i+j-1) - x(i));
    end
end
coef_newton = tabla(1,:);

% --- Evaluar todos los puntos para las curvas ---
for k = 1:length(xx)
    x_val = xx(k);
    
    % Evaluar Lagrange
    y_lag = 0;
    for i = 1:n
        L = 1;
        for j = 1:n
            if j ~= i
                L = L * (x_val - x(j)) / (x(i) - x(j));
            end
        end
        y_lag = y_lag + L * y(i);
    end
    yy_lagrange(k) = y_lag;
    
    % Evaluar Newton
    y_newt = coef_newton(1);
    p = 1;
    for i = 2:n
        p = p * (x_val - x(i-1));
        y_newt = y_newt + coef_newton(i) * p;
    end
    yy_newton(k) = y_newt;
end

% =========================================================================
% 4. GRAFICACIÓN DE LOS POLINOMIOS
% =========================================================================
figure('Name', 'Interpolación Polinómica', 'Color', 'w');

% Curvas
plot(xx, yy_lagrange, 'r-', 'LineWidth', 5, 'DisplayName', 'Polinomio Lagrange'); hold on;
plot(xx, yy_newton, 'b--', 'LineWidth', 2, 'DisplayName', 'Polinomio Newton');

% Puntos
plot(x, y, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7, 'DisplayName', 'Datos Originales');

legend('Location', 'northwest', 'FontSize', 11);
title('Interpolación: Método de Lagrange vs Newton', 'FontSize', 14);
xlabel('Eje X', 'FontSize', 12); 
ylabel('Eje Y (f(x))', 'FontSize', 12);
grid on;

% =========================================================================
% 5. EVALUACIÓN EN x = 3.5
% =========================================================================
x_eval = 3.5;

% Evaluar Lagrange en 3.5
val_lagrange = 0;
for i = 1:n
    L = 1;
    for j = 1:n
        if j ~= i
            L = L * (x_eval - x(j)) / (x(i) - x(j));
        end
    end
    val_lagrange = val_lagrange + L * y(i);
end

% Evaluar Newton en 3.5
val_newton = coef_newton(1);
p = 1;
for i = 2:n
    p = p * (x_eval - x(i-1));
    val_newton = val_newton + coef_newton(i) * p;
end

% Imprimir resultados
fprintf('\n==================================================\n');
fprintf('      RESULTADOS DE LA EVALUACIÓN EN x = 3.5      \n');
fprintf('==================================================\n');
fprintf('Valor interpolado con Lagrange : %.6f\n', val_lagrange);
fprintf('Valor interpolado con Newton   : %.6f\n', val_newton);
fprintf('==================================================\n\n');