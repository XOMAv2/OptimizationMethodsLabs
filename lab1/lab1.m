function lab1()
    a = 0;
    b = 1;
    epsArr = [0.01 0.0001 0.000001];
    
    fprintf('eps          | x*          | f* = f(x*)  | N           \n');
    fprintf('-------------|-------------|-------------|-------------\n');
    
    for eps = epsArr
        [xRes, yRes, xCalc] = bitwiseSearch(a, b, @targetFunc, eps);
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
  y = sin((x.^4 + x.^3 - 3 * x + 3 - 30.^(1 / 3)) / 2) + ...
      tanh((4 * 3.^0.5 * x.^3 - 2 * x - 6 * 2.^0.5 + 1) /  ...
           (-2 * 3.^0.5 * x.^3 + x + 3 * 2.^0.5)) + ...
      1.2;
end

function [xRes, yRes, xCalc] = bitwiseSearch(a, b, f, eps)
    delta = (b - a) / 4;

    x0 = a;
    f0 = f(x0);
    
    xCalc = [x0];

    x1 = x0;
    f1 = f0;

    while abs(delta) > eps
        x0 = x1;
        f0 = f1;

        x1 = x0 + delta;
        f1 = f(x1);
        
        xCalc = [xCalc x1];
        
        if f1 >= f0 || x1 <= a || x1 >= b
            delta = -delta / 4;
        end
    end

    xRes = x1;
    yRes = f1;
end