function [rho] = section_doree(U,d,P,C,D,alpha,beta)
tau=(sqrt(5)-1)/2; %ratio de la section doree
x1=a + (1- tau)*(b-a); %1er point d'evaluation
x2=a+ tau * (b-a);%2eme point d'evaluation
%evaluation de la fonction en x1 et x2
f1=calcul_JN_gradJN(U+x1*d, P,C,D,alpha,beta);
f2=calcul_JN_gradJN(U+x2*d, P,C,D,alpha,beta);

%la boucle
while abs(b-a)> 1e-6 %on arrete si l'intervalle est de mesure 10^-6 donc petit
    if f1< f2
        %le minimum est dans [a,x2]
        b=x2;
        x2=x1;
        f2=f1;
        x1 = a + (1-tau)*(b-a);
        f1=calcul_JN_gradJN(U+x1*d, P,C,D,alpha,beta);
    else
        %le minimum est dans [x1,b]
        a=x1;
        x1=x2;
        f1=f2;
        x2= a + tau*(b-a);
        f2= calcul_JN_gradJN(U+x2*d,P,C,D,alpha,beta);
    end
end
%enfin on retourne le milieu de l'intervalle final comme pas optimal
rho =(a+b)/2;
end

