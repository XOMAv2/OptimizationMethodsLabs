function lab3()
    a = 0;
    b = 1;
    epsArr = [0.01 0.0001 0.000001];
    
    fprintf('eps          | x*          | f* = f(x*)  | N           \n');
    fprintf('-------------|-------------|-------------|-------------\n');
    
    for eps = epsArr
        [xRes, yRes, xCalc] = parabolasMethod(a, b, @targetFunc, eps);
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

%% Метод парабол.
function [xRes, yRes, xCalc] = parabolasMethod(a, b, f, eps)
	[x1, y1, x2, y2, x3, y3] = getBasePoints(a, b, f);
    xCalc = [x1 x2 x3];
    
    isFirstIteration = 1;
    
    while 1
        xAvg = getXAvg(x1, y1, x2, y2, x3, y3);
        
        if isFirstIteration
            isFirstIteration = 0;
        else
            if abs(xAvg - xPrev) <= eps
                xRes = xAvg;
                yRes = f(xAvg);
                xCalc = [xCalc xAvg];
                break;
            end
        end
        
        xPrev = xAvg;
        
        if x2 == xAvg
            % Вместо xAvg можно взять любую точку интервала (x1; x3)
            % например, (x1 + x2) / 2 или (x2 + x3) / 2.
            xAvg = (x1 + x2) / 2;
        end

        yAvg = f(xAvg);
        xCalc = [xCalc xAvg];

        [x1, y1, x2, y2, x3, y3] = ...
            getNextPoints(x1, y1, x2, y2, x3, y3, xAvg, yAvg);
    end
end

%% Вычисление точки минимума xAvg квадратного трёхчлена.
function xAvg = getXAvg(x1, f1, x2, f2, x3, f3)
	a1 = (f2 - f1) / (x2 - x1);
    a2 = ((f3 - f1) / (x3 - x1) - a1) / (x3 - x2);
	xAvg = (x1 + x2 - a1 / a2) / 2;
end

%% Выбор "базовых" точек с использованием метода золотого сечения.
function [x1, y1, x2, y2, x3, y3] = getBasePoints(a, b, f)
    tau = (sqrt(5) - 1) / 2;
    l = b - a;
    
    x1 = b - tau * l;
    x2 = a + tau * l;
    
    f1 = f(x1);
    f2 = f(x2);
    
    if f2 > f1   % Переход к отрезку [a; x2].
        [x1, y1, x2, y2, x3, y3] = deal(a, f(a), x1, f1, x2, f2);
    else         % Переход к отрезку [x1; b].
        [x1, y1, x2, y2, x3, y3] = deal(x1, f1, x2, f2, b, f(b));
    end
end

%% Определение тройки чисел x1, x2, x3 для второй и последующих итераций.
% Используется метод исключения отрезков.
% Делается предположение, что x2 не равен xAvg.
function [x1, y1, x2, y2, x3, y3] = ...
    getNextPoints(x_1, y_1, x_2, y_2, x_3, y_3, xAvg, yAvg)

    if x_2 > xAvg
        if y_2 > yAvg
            [x1, y1, x2, y2, x3, y3] = deal(x_1, y_1, xAvg, yAvg, x_2, y_2);
        else
            [x1, y1, x2, y2, x3, y3] = deal(xAvg, yAvg, x_2, y_2, x_3, y_3);
        end
    else
        if yAvg > y_2
            [x1, y1, x2, y2, x3, y3] = deal(x_1, y_1, x_2, y_2, xAvg, yAvg);
        else
            [x1, y1, x2, y2, x3, y3] = deal(x_2, y_2, xAvg, yAvg, x_3, y_3);
        end
    end
end