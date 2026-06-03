function [x,valFunz] = newtonStandard(f,f1,x0,tol,itmax)

% [x,valFunz] = newtonStandard(f,f1,x0,tol,itmax)
% 
% Input:
% f = funzione di partenza
% f1 = derivata prima della funzione f
% x0 = punto di partenza
% tol = tolleranza (default 1e-12)
% itmax = numero di iterazioni massime (default 1000)
% 
% Output:
% x = approssimazione dello zero della funzione
% valFunz = valutazioni funzionali svolte durante l'esecuzione 
% del codice
%
%                                         Rel. 2-11-2025

if nargin < 5
    itmax = 1000;
elseif itmax<1 
    error('maxit errato');
end

if nargin < 4
     tol = 1e-12;
elseif tol < 0 
    error('tolleranza negativa');
end

if nargin < 3
    error('numero argomenti di ingresso errato');
end

for i = 1:itmax
    f1x = f1(x0);
    if f1x==0 
        error('il metodo non converge');
    end
    fx = f(x0);
    x = x0-fx/f1x;
    err = abs(x-x0);
    if err <= tol 
        break; 
    end
    x0 = x;
 end
 if err > tol
     warning('tolleranza richiesta non soddisfatta');
 end
 valFunz= i*2;
 return

