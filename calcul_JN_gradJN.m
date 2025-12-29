function [J,gradJ] = calcul_JN_gradJN(U,P,C,D,alpha,beta)
N=size(P,1);
J=0;
L=size(C,1);
gradJ=zeros(N,1);
for l=1:L
     %chaque ligne du tableau C correspond a un triangle
     %chaque ligne du tableau P correspond a un sommet
     %notre variable sommets est un tableau des sommets de Tl
    sommets=P(C(l,:),:);

    %U_l=Ukl du projet
    U_l=U(C(l,:));

    %calcul de gradients des Tau k l
    %calcul des vecteurs normals
    n1=[sommets(2,2)-sommets(1,2); sommets(1,1)-sommets(2,1)];
    n2=[sommets(3,2)-sommets(2,2); sommets(2,1)-sommets(3,1)];
    n3=[sommets(3,2)-sommets(1,2); sommets(1,1)-sommets(3,1)];

    
    %dot=produit scalaire
    %le calcul des gradients a été transformé en simplifiant le calcul avec
    %le facteur M
    M=dot(sommets(1,:)-sommets(2,:),n3');
   %M ne peut pas etre nulle sinon on a le cas d'un triangle plat->probleme
   %avec la bijection
    grad_tau=(1/M).*[n2';-n3';-n1'];

    %calcul de Sl
    grad_uh= sum(U_l.* grad_tau,1);
    Sl=norm(grad_uh)^2;


    %Calcul de ml mesure du triangle Tl;
    ml=polyarea(sommets(:,1),sommets(:,2));
    %autrement : ml=0.5*norm(M); mais polyarea est deja optimisée par
    %Matlab
    
    %caclul de J
    J=J+ ml*(alpha*sqrt(1+Sl)+beta*Sl);

    %Calcul du gradient
    for k=1:3
        i=C(l,k);
        if ~ismember(i,D) %le sommet n'est pas au bord
            %derivee de Sl
            dS_dUk = 2*dot(grad_uh,grad_tau(k,:));
            gradJ(i) = gradJ(i) + area *(alpha * dS_dUk /sqrt(1+Sl)+ 2 * beta * dS_dUk);
            
        end
    end
end