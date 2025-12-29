function [P,C,D] = maillagerectangle(a,b,c,d,Lx,Ly)
x=a:(b-a)/Lx:b;
y=c:(b-a)/Ly:d;
[X,Y]=meshgrid(x,y);
N=size(X,1)*size(X,2);
Xr=reshape(X,N,1);
Yr=reshape(Y,N,1);
DT=delaunayTriangulation(Xr,Yr);
P=DT.Points;
C=DT.ConnectivityList;
Q=P(:,1)==b |P(:,1)== a| P(:,2)==d |P(:,2)==c;
D=find(Q);

end