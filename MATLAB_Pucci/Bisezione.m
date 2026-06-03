function [x,valFunz]= Bisezione(f,a,b,tol) 

% bisezione: metodo della bisezione per trovare uno zero di f(x) in [a, b]
%
%   Input:
%   f = funzione anonima 
%   a,b = estremi dell'intervallo iniziale [a, b], con f(a)*f(b) < 0
%   tol       - (opzionale) tolleranza (default: 1e-12)
% 
%   Output:
%   x = approssimazione dello zero
%   valutazioniFunz = numero di iterazioni effettuate
%
%                                                  Rel. 2-11-2025


    if nargin == 3
        tol =1e-12;
    elseif nargin<3
            error('dati in input non presenti');
    end

    if b<= a 
        error('intervallo iniziale errato');
    end

    if tol <= 0
        error ('tolleranza errata'); 
    end

    fa = f(a); 
        valFunz=1;
    if fa == 0 
        x=a; 
        return; 
    end
 
    fb = f(b);
       valFunz=valFunz+1;
    if fb == 0
        x=b;
        return;
    end
    if fa*fb > 0
        error('metodo non applicabile');
    end
    maxit =ceil(-log2(tol / (b- a)));

    for it = 1:maxit
        x=(a+b)/2;
        fx = f(x);
        f1 = abs(fb-fa)/(b-a);
        if abs(fx) <= f1*tol 
            break;
        
        elseif fa*fx<0
                b=x;
                fb=fx;
        else
            a=x;
            fa=fx;
        end
        
    end
    valFunz=valFunz + it;  
    return 
    end

