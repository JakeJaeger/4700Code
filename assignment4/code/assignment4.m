%% Question 1: PA9 Work
% 
% In order to model the circuit with a finite difference model, we need
% the differential equations which represent the circuit
% 
% $V_1 = V_{in}$
% 
% $g_1(V_2-V_1)+C\frac{d(V_2-V_1)}{dt}+G_2V_2-I_L = 0$
% 
% $V_2-V_3- L \frac{dI_L}{dt} = 0$
% 
% $-I_L+g_3V_3 = 0$
% 
% $V_4-\alpha I_3 = 0$
% 
% $G_3V_3- I_3 = 0$
% 
% $G_4(V_O-V_4)+G_OV_O = 0$
% 
% Rewritten in the frequency domain these become
% 
% $V_1 = V_{in}$
% 
% $g_1(V_2-V_1)+j\omega C(V_2-V_1)+G_2V_2-I_L = 0$
% 
% $V_2-V_3- L j\omega I_L = 0$
% 
% $-I_L+g_3V_3 = 0$
% 
% $V_4-\alpha I_3 = 0$
% 
% $G_3V_3- I_3 = 0$
% 
% $G_4(V_O-V_4)+G_OV_O = 0$
% 
% the resultant C and G matrices, along with the F vector can be seen in
% the code block below. the circuit can then be modelled using
% 
% $$(G+j\omega C)V=F(\omega )$$
% 
% 
close all
g1=1;
g2=1/2;
g3=1/10;
g4=1/0.1;
g5=1/1000;
a=100;
c=0.25;
L=0.2;
%V=[V1;V2;Il;V3;I3;V4;V0]
G=[1,0,0,0,0,0,0;
  -g2,g1+g2,-1,0,0,0,0;
  0,1,0,-1,0,0,0;
  0,0,-1,g3,0,0,0;
  0,0,0,0,-a,1,0;
  0,0,0,g3,-1,0,0;
  0,0,0,0,0,-g4,g4+g5]

C=[0,0,0,0,0,0,0;
    -c,c,0,0,0,0,0;
    0,0,-L,0,0,0,0;
    0,0,0,0,0,0,0;
    0,0,0,0,0,0,0;
    0,0,0,0,0,0,0;
    0,0,0,0,0,0,0]


plotV0=[];
plotV3=[];
plotVin=[];
for Vin=-10:10
F=[Vin;0;0;0;0;0;0];    
V=G\F;
plotV0=[plotV0,V(7)];
plotV3=[plotV3,V(4)];
plotVin=[plotVin,Vin];
end
figure(1)
hold on
plot (plotVin,plotV0)
plot (plotVin,plotV3)
title("Figure 1: DC Sweep of input voltage")
xlabel("Input Voltage (V)")
ylabel("Output Voltage (V)")
legend("Vout","V3")

plotV0=[];
plotw=[];
Vin=1;
F=[Vin;0;0;0;0;0;0];
for w=0:10000  
    V=(G+C*1i*w)\F;
    plotV0=[plotV0,V(7)];
    plotw=[plotw,w];
end
figure (2)
semilogx (plotw,abs(plotV0));
title("Figure 2: Output voltage over an AC sweep")
xlabel("frequency (rad/s)")
ylabel("Vout (V)")

figure (3)
semilogx (plotw,20*log10(abs(plotV0)))
title("Figure 3: Gain over an AC sweep")
xlabel("frequency (rad/s)")
ylabel("Gain (dB)")

plotgain=[];
plotc=[];
for i=1:2000
    c=0.25+0.05*randn();
    C=[0,0,0,0,0,0,0;
    -c,c,0,0,0,0,0;
    0,0,-L,0,0,0,0;
    0,0,0,0,0,0,0;
    0,0,0,0,0,0,0;
    0,0,0,0,0,0,0;
    0,0,0,0,0,0,0];
    V=(G+C*1i*pi)\F;
    plotgain=[plotgain,20*log10(abs(V(7)))];
    plotc=[plotc,c];
end

figure (4)
histogram(plotgain)
title("Figure 4: Histogram of Gain for random pertubations of Frequency")
xlabel("Frequency (rad/s)")
ylabel("Count")


%% Question 2 Transient Circuit Simulation
% 
% The circuit will act as a sort of Band Pass filter, with the low side
% rolloff slowed by the R1 in parallel with the capacitor.  thus we expect
% the gain to decrease as the frequency increases above the centre
% frequency.  The current
% controlled voltage source is a source of gain, giving the circuit the
% properties of an amplifier as well. Using a numerical finite difference
% solution, the circuit is simulated in the code block below for a step,
% sinusiodal, and gaussian input.  The input and output voltages for these
% simulations are shown in figures 5, 7, and 9. The signals were also put
% through a fourier transform to get their frequency domain
% representations.  The magnitude plots are shown in figures 6, 8, and 10.
% when the time step is increased by a factor of 10, the simulation becomes
% very innacurate.  This can be seen by comparing figure 11 to figure 9.  


clear
g1=1;
g2=1/2;
g3=1/10;
g4=1/0.1;
g5=1/1000;
a=100;
c=0.25;
L=0.2;

G=[1,0,0,0,0,0,0;
  -g2,g1+g2,-1,0,0,0,0;
  0,1,0,-1,0,0,0;
  0,0,-1,g3,0,0,0;
  0,0,0,0,-a,1,0;
  0,0,0,g3,-1,0,0;
  0,0,0,0,0,-g4,g4+g5];

C=[0,0,0,0,0,0,0;
    -c,c,0,0,0,0,0;
    0,0,-L,0,0,0,0;
    0,0,0,0,0,0,0;
    0,0,0,0,0,0,0;
    0,0,0,0,0,0,0;
    0,0,0,0,0,0,0];

%V=[V1;V2;Il;V3;I3;V4;V0]
t=linspace(0,1,1000);
dt=t(2)-t(1);
V=zeros(7,1000);
plotV0=[];
plotVin=[];
Vin=0;
for i=1:1000
    if t(i)>0.03
    Vin=1;
    end
    F=[Vin;0;0;0;0;0;0];
    if i==1
        V(:,i)=((C/dt + G)^-1)*(F);
    else
        V(:,i)=((C/dt + G)^-1)*(C*(V(:,i-1)/dt)+F);
    end
    plotV0=[plotV0,V(7,i)];
    plotVin=[plotVin,Vin];
end
figure(5)
plot(t,plotV0)
hold on
plot(t,plotVin)
title("Figure 5: Voltage over Time for Step Input")
xlabel("Time (s)")
ylabel("Voltage (V)")
legend("Vout","Vin")

figure(6)
semilogy(linspace(-500,500,1000),fftshift(abs(fft(plotV0))))
hold on
semilogy(linspace(-500,500,1000),fftshift(abs(fft(plotVin))))
title("Figure 6: Magnitude in Frequency Domain for Step Input")
xlabel("Frequency (Hz)")
ylabel("Magnitude (V)")
legend("output","input")

V=zeros(7,1000);
plotV0=[];
plotVin=[];
Vin=0;
f=1/(0.03); %frequency of sin
for i=1:1000
    Vin=sin(2*pi*f*t(i));
    F=[Vin;0;0;0;0;0;0];
    if i==1
        V(:,i)=((C/dt + G)^-1)*(F);
    else
        V(:,i)=((C/dt + G)^-1)*(C*(V(:,i-1)/dt)+F);
    end
    plotV0=[plotV0,V(7,i)];
    plotVin=[plotVin,Vin];
end
figure(7)
plot(t,plotV0)
hold on
plot(t,plotVin)
title("Figure 7: Voltage over Time for sin Input")
xlabel("Time (s)")
ylabel("Voltage (V)")
legend("Vout","Vin")

figure(8)
semilogy(linspace(-500,500,1000),fftshift(abs(fft(plotV0))))
hold on
semilogy(linspace(-500,500,1000),fftshift(abs(fft(plotVin))))
title("Figure 8: Magnitude in Frequency Domain for sin Input")
xlabel("Frequency (Hz)")
ylabel("Magnitude (V)")
legend("output","input")

V=zeros(7,1000);
plotV0=[];
plotVin=[];
for i=1:1000
    Vin=exp(-((t(i)-0.06)^2)/(2*0.03^2));
    F=[Vin;0;0;0;0;0;0];
    if i==1
        V(:,i)=((C/dt + G)^-1)*(F);
    else
        V(:,i)=((C/dt + G)^-1)*(C*(V(:,i-1)/dt)+F);
    end
    plotV0=[plotV0,V(7,i)];
    plotVin=[plotVin,Vin];
end
figure(9)
plot(t,plotV0)
hold on
plot(t,plotVin)
title("Figure 9: Voltage over Time for Gaussian Input")
xlabel("Time (s)")
ylabel("Voltage (V)")
legend("Vout","Vin")

figure(10)
semilogy(linspace(-500,500,1000),fftshift(abs(fft(plotV0))))
hold on
semilogy(linspace(-500,500,1000),fftshift(abs(fft(plotVin))))
title("Figure 10: Magnitude in Frequency Domain for Gaussian Input")
xlabel("Frequency (Hz)")
ylabel("Magnitude (V)")
legend("output","input")
V=zeros(7,500);
plotV0=[];
plotVin=[];
Vin=0;
f=1/(0.03); %frequency of sin
for i=1:100
    Vin=exp(-((t(i)-0.06)^2)/(2*0.03^2));
    F=[Vin;0;0;0;0;0;0];
    if i==1
        V(:,i)=((C/dt + G)^-1)*(F);
    else
        V(:,i)=((C/(10*dt) + G)^-1)*(C*(V(:,i-1)/(10*dt))+F);
    end
    plotV0=[plotV0,V(7,i)];
    plotVin=[plotVin,Vin];
end
figure(11)
newt=linspace(0,1,100);
plot(newt,plotV0)
hold on
plot(newt,plotVin)
title("Figure 11: Voltage over Time for Gaussian Input with Reduced Timestep")
xlabel("Time (s)")
ylabel("Voltage (V)")
legend("Vout","Vin")

%% Question 3: Circuit with Noise
% 
% To simulate thermal noise in R3, we add a current source and capacitor in
% parralel to it.  The current source simulates the noise, while the
% capacator serves to bandwidth limit that noise.  This requires an extra
% equation to be added, adding an extra line to the C and G matricies.  the
% new matricies are included in the code block below, and output into the
% report.  Figure 12 shows the input and output voltage for a gaussian
% input.  Figure 13 shows the fourier transforms of those signals.
% increasing Cn lowers the noise level, this can be most clearly seen in
% the fourier transform plots in figures 13, 15, and 17. The voltage plots 
% are in figures 12, 14, and 16. Increasing the time step once again causes
% the signals to be very innacurate.  It also has the additional effect of
% lessening the frequency of the noise, as can be seen by comparing figure
% 18 to figure 12.  

clear
g1=1;
g2=1/2;
g3=1/10;
g4=1/0.1;
g5=1/1000;
a=100;
c=0.25;
L=0.2;
cn=0.00001;

G=[1,0,0,0,0,0,0,0;
  -g2,g1+g2,-1,0,0,0,0,0;
  0,1,0,-1,0,0,0,0;
  0,0,-1,g3,0,0,0,-1;
  0,0,0,0,-a,1,0,0;
  0,0,0,g3,-1,0,0,0;
  0,0,0,0,0,-g4,g4+g5,0;
  0,0,0,0,0,0,0,1]

C=[0,0,0,0,0,0,0,0;
    -c,c,0,0,0,0,0,0;
    0,0,-L,0,0,0,0,0;
    0,0,0,cn,0,0,0,0;
    0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0]

%V=[V1;V2;Il;V3;I3;V4;V0;In]
t=linspace(0,1,1000);
dt=t(2)-t(1);
V=zeros(8,1000);
plotV0=[];
plotVin=[];
for i=1:1000
    Vin=exp(-((t(i)-0.06)^2)/(2*0.03^2));
    In=0.001*randn();
    F=[Vin;0;0;0;0;0;0;In];
    if i==1
        V(:,i)=((C/dt + G)^-1)*(F);
    else
        V(:,i)=((C/dt + G)^-1)*(C*(V(:,i-1)/dt)+F);
    end
    plotV0=[plotV0,V(7,i)];
    plotVin=[plotVin,Vin];
end
figure(12)
plot(t,plotV0)
hold on
plot(t,plotVin)
title({"Figure 12: Voltage over Time for Gaussian Input with noise", "Cn=0.00001"})
xlabel("Time (s)")
ylabel("Voltage (V)")
legend("Vout","Vin")
figure(13)
semilogy(linspace(-500,500,1000),fftshift(abs(fft(plotV0))))
hold on
semilogy(linspace(-500,500,1000),fftshift(abs(fft(plotVin))))
title({"Figure 13: Magnitude in Frequency Domain for Gaussian Input with noise","Cn=0.00001"})
xlabel("Frequency (Hz)")
ylabel("Magnitude (V)")
legend("output","input")

cn=0.00005;
C=[0,0,0,0,0,0,0,0;
    -c,c,0,0,0,0,0,0;
    0,0,-L,0,0,0,0,0;
    0,0,0,cn,0,0,0,0;
    0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0];
V=zeros(8,1000);
plotV0=[];
plotVin=[];
for i=1:1000
    Vin=exp(-((t(i)-0.06)^2)/(2*0.03^2));
    In=0.001*randn();
    F=[Vin;0;0;0;0;0;0;In];
    if i==1
        V(:,i)=((C/dt + G)^-1)*(F);
    else
        V(:,i)=((C/dt + G)^-1)*(C*(V(:,i-1)/dt)+F);
    end
    plotV0=[plotV0,V(7,i)];
    plotVin=[plotVin,Vin];
end
figure(14)
plot(t,plotV0)
hold on
plot(t,plotVin)
title({"Figure 14: Voltage over Time for Gaussian Input with noise","Cn=0.00005"})
xlabel("Time (s)")
ylabel("Voltage (V)")
legend("Vout","Vin")
figure(15)
semilogy(linspace(-500,500,1000),fftshift(abs(fft(plotV0))))
hold on
semilogy(linspace(-500,500,1000),fftshift(abs(fft(plotVin))))
title({"Figure 15: Magnitude in Frequency Domain for Gaussian Input with noise","Cn=0.00005"})
xlabel("Frequency (Hz)")
ylabel("Magnitude (V)")
legend("output","input")

cn=0.0005;
C=[0,0,0,0,0,0,0,0;
    -c,c,0,0,0,0,0,0;
    0,0,-L,0,0,0,0,0;
    0,0,0,cn,0,0,0,0;
    0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0];
V=zeros(8,1000);
plotV0=[];
plotVin=[];
for i=1:1000
    Vin=exp(-((t(i)-0.06)^2)/(2*0.03^2));
    In=0.001*randn();
    F=[Vin;0;0;0;0;0;0;In];
    if i==1
        V(:,i)=((C/dt + G)^-1)*(F);
    else
        V(:,i)=((C/dt + G)^-1)*(C*(V(:,i-1)/dt)+F);
    end
    plotV0=[plotV0,V(7,i)];
    plotVin=[plotVin,Vin];
end
figure(16)
plot(t,plotV0)
hold on
plot(t,plotVin)
title({"Figure 14: Voltage over Time for Gaussian Input with noise","Cn=0.0005"})
xlabel("Time (s)")
ylabel("Voltage (V)")
legend("Vout","Vin")
figure(17)
semilogy(linspace(-500,500,1000),fftshift(abs(fft(plotV0))))
hold on
semilogy(linspace(-500,500,1000),fftshift(abs(fft(plotVin))))
title({"Figure 17: Magnitude in Frequency Domain for Gaussian Input with noise","Cn=0.0005"})
xlabel("Frequency (Hz)")
ylabel("Magnitude (V)")
legend("output","input")

cn=0.00001;
C=[0,0,0,0,0,0,0,0;
    -c,c,0,0,0,0,0,0;
    0,0,-L,0,0,0,0,0;
    0,0,0,cn,0,0,0,0;
    0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0];
V=zeros(8,200);
t=linspace(0,1,200);
dt=t(2)-t(1);
plotV0=[];
plotVin=[];
for i=1:200
    Vin=exp(-((t(i)-0.06)^2)/(2*0.03^2));
    In=0.001*randn();
    F=[Vin;0;0;0;0;0;0;In];
    if i==1
        V(:,i)=((C/dt + G)^-1)*(F);
    else
        V(:,i)=((C/dt + G)^-1)*(C*(V(:,i-1)/dt)+F);
    end
    plotV0=[plotV0,V(7,i)];
    plotVin=[plotVin,Vin];
end
figure(18)
plot(t,plotV0)
hold on
plot(t,plotVin)
title({"Figure 18: Voltage over Time for Gaussian Input with noise","Reduced Timestep"})
xlabel("Time (s)")
ylabel("Voltage (V)")
legend("Vout","Vin")
figure(19)
semilogy(linspace(-500,500,200),fftshift(abs(fft(plotV0))))
hold on
semilogy(linspace(-500,500,200),fftshift(abs(fft(plotVin))))
title({"Figure 19: Magnitude in Frequency Domain for Gaussian Input with noise","Reduced Timestep"})
xlabel("Frequency (Hz)")
ylabel("Magnitude (V)")
legend("output","input")

%% Question 4: Non-linearity
% 
% If the current controlled voltage source on the output stage were
% non-linear, e.g modeled by $V=\alpha I_3+\beta I^2_3 +\gamma I^3_3$ then
% an extra non-linear component would need to be added to the equation to
% account for this.  in the matrix solution, a non-linear B-vector would be
% added such that the time domain equation becomes
% 
% $$C\frac{dV}{dt}+GV+B=F$$
% 
% with the V4 component of B being $\beta {(G_3V_3)}^2+\gamma {(G_3V_3)}^3$. B
% is a function of V3 at the current time, and thus in addition to iterating over time, we will
% need to iterate within each timestep to converge on a value for B.  
