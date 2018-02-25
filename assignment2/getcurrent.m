function [ current ] = getcurrent( meshSize,neckSize,sigo,sigi,plots)
W=meshSize;
L=round(3/2*W);
Wb=[round((W-neckSize)/2) round((W+neckSize)/2)];
Lb=[round(L/3) round(L-(L/3))];


G=sparse(L*W);
B=zeros(1,L*W);
sigmatrix=sigo*ones(W,L);
for i=Lb(1):Lb(2)
    for j=1:Wb(1)
        sigmatrix(j,i)=sigi;
    end
    for j=Wb(2):W
        sigmatrix(j,i)=sigi;
    end
end
for i =1:L
    for j=1:W
        n=j+(i-1)*W;
        if i==1
            G(n,:)=0;
            G(n,n)=1;
            B(n)=1;
        elseif i==L
            G(n,:)=0;
            G(n,n)=1;
            B(n)=0;
        elseif j==1
            Rup=(sigmatrix(j,i)+sigmatrix(j+1,i))/2;
            Rleft=(sigmatrix(j,i)+sigmatrix(j,i-1))/2;
            Rright=(sigmatrix(j,i)+sigmatrix(j,i+1))/2;
            G(n,n)=-(Rup+Rleft+Rright);
            G(n,n+1)=Rup;
            G(n,n+W)=Rright;
            G(n,n-W)=Rleft;
        elseif j==W
            Rdown=(sigmatrix(j,i)+sigmatrix(j-1,i))/2;
            Rleft=(sigmatrix(j,i)+sigmatrix(j,i-1))/2;
            Rright=(sigmatrix(j,i)+sigmatrix(j,i+1))/2;
            G(n,n)=-(Rdown+Rleft+Rright);
            G(n,n-1)=Rdown;
            G(n,n+W)=Rright;
            G(n,n-W)=Rleft;
        else
            Rup=(sigmatrix(j,i)+sigmatrix(j+1,i))/2;
            Rdown=(sigmatrix(j,i)+sigmatrix(j-1,i))/2;
            Rleft=(sigmatrix(j,i)+sigmatrix(j,i-1))/2;
            Rright=(sigmatrix(j,i)+sigmatrix(j,i+1))/2;
            G(n,n)=-(Rup+Rdown+Rleft+Rright);
            G(n,n+1)=Rup;
            G(n,n-1)=Rdown;
            G(n,n+W)=Rright;
            G(n,n-W)=Rleft;
        end       
    end
end
V=G\B';
Vmatrix=zeros(W,L);
for i =1:L
    for j=1:W
        n=j+(i-1)*W;
        Vmatrix(j,i)=V(n);
    end 
end


[Ex,Ey]=gradient(Vmatrix);
Ex=-Ex;
Ey=-Ey;





Jx=Ex.*sigmatrix;
Jy=Ey.*sigmatrix;


current=sum(Jx(:,1))/W;

if plots
figure
surf(Vmatrix)
title('Figure 4: V(x,y)')
xlabel('Length')
ylabel('Width')

figure 
quiver(Ex,Ey)
title('Figure 5: E(x,y)')
xlabel('Length')
ylabel('Width')

figure
surf(sigmatrix)
title('Figure 6: \sigma (x,y)')
xlabel('Length')
ylabel('Width')

figure
quiver(Jx,Jy)
title('Figure 7: J(x,y)')
xlabel('Length')
ylabel('Width')
end

end