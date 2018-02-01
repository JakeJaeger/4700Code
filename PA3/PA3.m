xmax=20;
ymax=20;
V=zeros(ymax,xmax);
steps=200;

%set boundary conditions
V(:,1)=1;
V(:,xmax)=1;
V(1,:)=0;
V(ymax,:)=0;

for i=1:steps
    Vold=V;
   for x=2:(xmax-1)
       %V(1,x)= (Vold(1,x-1)+Vold(1,x+1)+Vold(2,x))/3; %bottom
       %V(ymax,x)=(Vold(ymax,x-1)+Vold(ymax,x+1)+Vold(ymax-1,x))/3;%top
       for y=2:(ymax-1)
        V(x,y)= (Vold(x-1,y)+Vold(x+1,y)+Vold(x,y-1)+Vold(x,y+1))/4;
       end
   end
   figure(1)
   surf(V)
   pause(0.05)
end 
[Ex,Ey]=gradient(V);
figure(2)
quiver(Ex,Ey)
hold on
contour(V)