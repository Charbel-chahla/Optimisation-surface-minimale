clear all


function[P,C]=Triangulation(Lx,Ly)
x=0:1/Lx:Lx;
y=0:1/Ly:Ly;
[X,Y]=meshgrid(x,y);
N=size(X,1)
