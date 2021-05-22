function lab2()
    a = 0;
    b = 1;
    epsArr = [0.01 0.0001 0.000001];
    
    fprintf('eps          | x*          | f* = f(x*)  | N           \n');
    fprintf('-------------|-------------|-------------|-------------\n');
    
    for eps = epsArr
        [xRes, yRes, xCalc] = goldenRatio(a, b, @targetFunc, eps);
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
        plot(xRes, yRes, 'dm','LineWidth', 5, 'MarkerSize', 10);
        legend('Целевая функция f(x)', ...
               'Последовательность приближений (xi; f(xi))', ...
               'Найденная точка минимума (x*; f(x*))');
    end
end

function y = targetFunc(x)
    y = sin((x.^4 + x.^3 - 3 .* x + 3 - 30.^(1 ./ 3)) ./ 2) + ...
        tanh((4 .* 3.^0.5 .* x.^3 - 2 .* x - 6 .* 2.^0.5 + 1) ./  ...
             (-2 .* 3.^0.5 .* x.^3 + x + 3 .* 2.^0.5)) + ...
        1.2;
end

% Функция написана по текстовому описанию алгоритма.
function [xRes, yRes, xCalc] = goldenRatio(a, b, f, eps)
    tau = (sqrt(5) - 1) / 2;
    eps_n = (b - a) / 2;
    
    x1 = (a - 2 * tau + 2) * eps_n;
    x2 = (a + 2 * tau) * eps_n;
    
    f1 = f(x1);
    f2 = f(x2);
    
    xCalc = [x1 x2];
    
    while eps_n > eps
        if f1 <= f2
            % Переход к отрезку [a; x2].
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = b - tau * (b - a);
            f1 = f(x1);
            xCalc = [xCalc x1];
        else
            % Переход к отрезку [x1; b].
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + tau * (b - a);
            f2 = f(x2);
            xCalc = [xCalc x2];
        end
        
        eps_n = tau * eps_n;
    end
    
    xRes = (b + a) / 2;
    yRes = f(xRes);
    xCalc = [xCalc xRes];
end

% Функция написана по схеме алгоритма.
function [xRes, yRes, xCalc] = goldenRatio1(a, b, f, eps)
    tau = (sqrt(5) - 1) / 2;
    l = b - a;
    
    x1 = b - tau * l;
    x2 = a + tau * l;
    
    f1 = f(x1);
    f2 = f(x2);
    
    xCalc = [x1 x2];
    
    while l >= 2 * eps
        if f1 < f2
            % Переход к отрезку [a; x2].
            b = x2;
            l = b - a;
            x2 = x1;
            f2 = f1;
            x1 = b - tau * l;
            f1 = f(x1);
            xCalc = [xCalc x1];
        else
            % Переход к отрезку [x1; b].
            a = x1;
            l = b - a;
            x1 = x2;
            f1 = f2;
            x2 = a + tau * l;
            f2 = f(x2);
            xCalc = [xCalc x2];
        end
    end
    
    xRes = (b + a) / 2;
    yRes = f(xRes);
    xCalc = [xCalc xRes];
end