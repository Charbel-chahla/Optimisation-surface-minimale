function U = descente(P,C,g,alpha,beta,epsilon,max)
N=size(P,1);
U=zeros(N,1);
%initialisation
for i=D
    U(i)=g(P(i,1),P(i,2));
end
for k=1:max
    [J,gradJ]=calcul_JN_gradJN(U,P,C,D,alpha,beta);
    if norm(gradJ) < epsilon
        break;
    end
    rho= section_doree(U,-gradJ,P,C,D,alpha,beta); %Methode de la section doree
    U=U - rho*gradJ;
end
end

