clear
close all
x=0;
v=0;
a=1;
tmax=1000;
t=0;
while t<tmax
    xp=x;
    vp=v;
    if rand()>0.05
        v=v+a;
    else 
        v=0;
    end
    x=x+v;
    driftv=x/t;
    plot([x,xp],[v,vp],'r')
    hold on
    title(['Drift Velocity: ',num2str(driftv)])
    pause(0.1)
    t=t+1;
end


