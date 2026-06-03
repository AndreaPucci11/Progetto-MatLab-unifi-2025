%   x = CHEBYSHEV(a, b, n) 
%
%   Input :
%     a,b  = estremi dell'intervallo [a, b]
%     n    = grado del polinomio (numero di suddivisioni - 1)
%
%   Output:
%     x    = nodi di Chebyshev nell'intervallo [a, b]
%                                           2-11-2025
function x = chebyshev(a, b, n)

    if n < 0
        error('Il parametro n deve essere un intero >= 0');
    end
    
    x_cheb_std = zeros(n + 1, 1);   
    x = zeros(n + 1, 1);            
    
    for k = 0:n
        x_cheb_std(k + 1) = cos(((2 * k + 1) * pi) / (2 * (n + 1)));
    end
    
    for k = 0:n
        x(k + 1) = ((a + b) / 2) + ((b - a) / 2) * x_cheb_std(k + 1);
    end

end
