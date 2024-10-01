clc
clear all
close all
% bus line from  to   R        X            B
   lidata= [1    2   0.02     0.06          0.03
            1    3   0.08     0.24          0.025
            2    3   0.06     0.25          0.02
            2    4   0.06     0.18          0.02
            2    5   0.04     0.12          0.015
            3    4   0.01     0.03          0.01
            4    5   0.08     0.24          0.025] ;
      

nf=(lidata(:,1));            %line form bus 
nt=(lidata(:,2));            %line to bus
nbus=max(max(nf),max(nt));   %taking max no of bus 
lines=length(nf) ;           % no of total lines
R=lidata(:,3) ;              %Resistance from line data
X=lidata(:,4);               %reactance from line data
B=complex(0,lidata(:,5)) ;   %line charging from line data

Z=complex(R,X) ;             %calculate impedance of lines 
Y=ones(lines,1)  ;           %initiaize admittance of the lines
Y=Y./Z ;                     %calculatee admittance of the line 1/z
Ybus=zeros(nbus,nbus) ;      %initialize n bus (5x5)

 %calculate non diagonal element
for k=1:lines                  
    Ybus(nf(k),nt(k))=-Y(k);
    Ybus(nt(k),nf(k))=-Y(k);
end

% calulate for diagonal element
for n= 1:lines
    for k=1:lines
        if nf(k)==n || nt(k)==n
            Ybus(n,n)=Ybus(n,n)+Y(k)+B(k);
        end
    end
end

Ybus; %display Ybus

%inatialize bus data

%       Bus Bus  Voltage Angle   ---Load---- -------Generator---
%        No  code Mag.    Degree  MW    Mvar  MW   Mvar 
busdata=[1   1    1.06    0.0    00.0  00.0  00.0  00.0   
         2   0    1.0     0.0    20.0  10.0  40.0  30.0  
         3   0    1.0     0.0    45.0  15.0  00.0  00.0  
         4   0    1.0     0.0    40.0  05.0  00.0  00.0   
         5   0    1.0     0.0    60.0  10.0  00.0  00.0];

% Declare Matrix
Sb=100;
S=zeros(nbus,1); %initialize complex power
% for calculating schedule data (generation - load)
pl=busdata(:,5)./Sb ;%load real power
ql=busdata(:,6)./Sb ;%load reactive power
pg=busdata(:,7)./Sb; %generation real power
qg=busdata(:,8)./Sb ;%generation reactive power
S=complex((pg-pl),(qg-ql)); % complex schedule power of buses 
V=complex(busdata(:,3),0); % bus voltages

V_old = zeros(nbus);% store the privious itteration

tolerance= 0.000000000000000000001;
max_err= 5000;
itt=0;

while (max_err>tolerance)           % gauss seidal 

    for i=1:nbus
        V_old(i)=V(i);
        sum=0;
        for k=1:nbus
            
            if(k~=i)
                sum=sum+Ybus(i,k)*V(k);
            end     
        end
        if i~=1
        V(i) = ((conj(S(i))/conj(V(i)))-sum)/Ybus(i,i);
        err(i) = abs(V(i)-V_old(i));
        end
    end
    max_err= max(err);
    itt=itt+1;
end
V;
theta=angle(V);
disp("Voltage Magnitude");
disp(abs(V))
disp("Angle");
disp(rad2deg(angle(V)))
% Calculat the real power and reactive power


%Power to line
for k = 1:nbus
    from_bus = nf(k);
    to_bus = nt(k);
    Yline = Y(k)+0.5*B(k);
    I_line = Yline * (V(from_bus) - V(to_bus));
    S_line = V(from_bus) * conj(I_line);  % Power flow in the line
    fprintf('From Bus %d to Bus %d: P = %.4f MW, Q = %.4f MVAR\n', from_bus, to_bus, real(S_line)*100, imag(S_line)*100);
end

% for schedule power calculation
P=zeros(nbus,1);
Q=zeros(nbus,1);

for i=1:nbus
    T_P=0;
    T_Q=0;
    for j=1:nbus
        T_P=T_P+abs(V(i))*abs(V(j))*abs(Ybus(i,j))*cos(angle(Ybus(i,j))+theta(j)-theta(i));
        T_Q=T_Q+abs(V(i))*abs(V(j))*abs(Ybus(i,j))*sin(angle(Ybus(i,j))+theta(j)-theta(i));
    end
    P(i)= T_P*100;
    Q(i)=-T_Q*100;
end

disp("");

disp("Real Power (Schedule)")
P
disp("Apparent Power (Schedule)")
Q
itt


    