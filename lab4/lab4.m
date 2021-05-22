function lab4()
    a = 0;
    b = 1;
    epsArr = [0.01 0.0001 0.000001];
    
    fprintf('eps          | x*          | f* = f(x*)  | N           \n');
    fprintf('-------------+-------------+-------------+-------------\n');
    fprintf('Модифицированный метод Ньютона с конечно-разностной    \n');
    fprintf('аппроксимацией производных                             \n');
    fprintf('-------------+-------------+-------------+-------------\n');
    
    for eps = epsArr
        [xRes, yRes, xCalc, xResList] = newton(a, b, @targetFunc, eps);
        xRange = a:eps:b;
        
        fprintf('%13.6f|', eps);
        fprintf('%13.6f|', xRes);
        fprintf('%13.6f|', yRes);
        fprintf('%13d\n', length(xCalc));
        
        % Команда figure создаёт новое (добавочное) графическое окно.
        figure('Name', 'eps = ' + string(eps)); 
        hold on; % Команда удержания текущего графическоро окна. hold on
                 % присваивает свойству NextPlot для текущих объектов
                 % figure и axes значение add.
        title('Метод поразрядного поиска');
        grid on;
        plot(xRange, targetFunc(xRange), '-g', 'LineWidth', 2);
        plot(xCalc, targetFunc(xCalc), 'xk', 'LineWidth', 2, 'MarkerSize', 10);
        plot(xResList, targetFunc(xResList), 'xr', 'LineWidth', 2, 'MarkerSize', 10);
        plot(xRes, yRes, 'dm','LineWidth', 5, 'MarkerSize', 10);
        legend('Целевая функция f(x)', ...
               'Точки расчёта целевой функции', ...
               'Последовательность приближений (xi; f(xi))', ...
               'Найденная точка минимума (x*; f(x*))');
    end
    
    fprintf('-------------+-------------+-------------+-------------\n');
    fprintf('Стандартная функция fminbnd пакета MatLAB              \n');
    fprintf('-------------+-------------+-------------+-------------\n');
    
    for eps = epsArr
        [x, fval, ~, output] = ...
            fminbnd(@targetFunc, a, b, optimset('TolX', eps));
        fprintf('%13.6f|', eps);
        fprintf('%13.6f|', x);
        fprintf('%13.6f|', fval);
        fprintf('%13d\n', output.iterations);
    end
end

%% Целевая функция
function y = targetFunc(x)
    y = sin((x.^4 + x.^3 - 3 .* x + 3 - 30.^(1 ./ 3)) ./ 2) + ...
        tanh((4 .* 3.^0.5 .* x.^3 - 2 .* x - 6 .* 2.^0.5 + 1) ./  ...
             (-2 .* 3.^0.5 .* x.^3 + x + 3 .* 2.^0.5)) + ...
        1.2;
end

%% Аппроксимация первой производной с использованием центральной разности.
function df = dTargetFunc(f_prev, f_next, delta)
    df = (f_next - f_prev) / 2 / delta;
end

%% Аппроксимация второй производной.
function d2f = d2TargetFunc(f_prev, f_middle, f_next, delta)
    d2f = (f_next - 2 * f_middle + f_prev) / delta / delta;
end

%% Модифицированный метод Ньютона (f''(x_k) = f''(x_0) = const).
function [xRes, yRes, xCalc, xResList] = newton(a, b, f, eps)
    % Для начального приближения используется метод золотого сечения.
	[x0, f0, xCalc] = getInitialApproximation(a, b, f);
    
    x_prev = x0 - eps;
    x_next = x0 + eps;
    f_prev = f(x_prev);
    f_next = f(x_next);
    xCalc = [xCalc x_prev x_next];
    d2x0 = d2TargetFunc(f_prev, f0, f_next, eps);
    
    xPrevStep = x0;
    xResList = [xPrevStep];
    
    isFirstIteration = 1;
    
    while 1
        if isFirstIteration
            dxPrevStep = dTargetFunc(f_prev, f_next, eps);
            isFirstIteration = 0;
        else
            x_prev = xPrevStep - eps;
            x_next = xPrevStep + eps;
            xCalc = [xCalc x_prev x_next];
            dxPrevStep = dTargetFunc(f(x_prev), f(x_next), eps);
        end
        
        xRes = xPrevStep - dxPrevStep / d2x0;
        xResList = [xResList xRes];
        
        if abs(xRes - xPrevStep) <= eps
            break
        else
            xPrevStep = xRes;
        end
    end
    
    yRes = f(xRes);
    xCalc = [xCalc xRes];
end

%% Выбор начального приближения с использованием метода золотого сечения.
function [xInit, fInit, xCalc] = getInitialApproximation(a, b, f)
    tau = (sqrt(5) - 1) / 2;
    l = b - a;
    
    x1 = b - tau * l;
    x2 = a + tau * l;
    
    f1 = f(x1);
    f2 = f(x2);
    
    xCalc = [x1 x2];
    
    if f2 > f1   % Переход к отрезку [a; x2].
        [xInit, fInit] = deal(x1, f1);
    else         % Переход к отрезку [x1; b].
        [xInit, fInit] = deal(x2, f2);
    end
end