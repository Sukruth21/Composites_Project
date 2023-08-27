clc
clear
close all

syms tht
fprintf('Program to Calculate the First Ply Failure Load of a given Symmetric Laminate using Classical Laminate Theory\n\n');

%Material Properties
E1=38.6e9;
E2=8.27e9;
v12=0.28;
G12=4.14e9;
t=0.125e-3;
n=8; %Number of Plies
v21=E2*v12/E1;

fprintf('Material and Geometric Properties:\n\n');
fprintf('E1 = %0.1f N/m^2\n',E1);
fprintf('E2 = %0.1f N/m^2\n',E2);
fprintf('mu = %0.3f\n',v12);
fprintf('G12 = %0.1f N/m^2\n',G12);
fprintf('Thickness of each ply = %f m\n',t);
fprintf('Number of plies = %d \n',n);

%Ultimate Strengths
sig1tu=1062e6;
sig1cu=610e6;
sig2tu=31e6;
sig2cu=118e6;
tau12u=72e6;

fprintf('\n\nUltimate Strengths of the Ply:\n\n');
fprintf('sigma1tu = %d N/m^2\n',sig1tu);
fprintf('sigma1cu = %d N/m^2\n',sig1cu);
fprintf('sigma2tu = %d N/m^2\n',sig2tu);
fprintf('sigma2cu = %d N/m^2\n',sig2cu);
fprintf('tau12u = %d N/m^2\n',tau12u);

%Angular Orientations of Plies in order
theta=[0 pi/4 -pi/4 pi/2];
thetaf=[0 pi/4 -pi/4 pi/2 pi/2 -pi/4 pi/4 0];
fprintf('\nOrientations (in degrees) of Plies in order from top to bottom:\n\n');
for i=1:n
    fprintf('\t%d',thetaf(i)*180/pi);
end
fprintf('\n');

%Load Vectors
N=[100e3;0;0];
M=[0;0;0];

%Positions of Different Plies
z=[-0.5e-3 -0.375e-3 -0.25e-3 -0.125e-3 0e-3 0.125e-3 0.25e-3 0.375e-3 0.5e-3];
Q11=E1/(1-v12*v21);
Q12=(v12*E2)/(1-v12*v21);
Q21=Q12;
Q22=E2/(1-v12*v21);
Q66=G12;

%Q matrix for the Laminate
Q=[Q11 Q12 0;
   Q21 Q22 0;
   0 0 Q66];
fprintf('\n\n Q matrix for the laminate:\n\n');
disp(Q);

%Transformation Matrix Definition
T=[cos(tht)^2 sin(tht)^2 sin(2*tht);
    sin(tht)^2 cos(tht)^2 -1*sin(2*tht);
    -0.5*sin(2*tht) 0.5*sin(2*tht) cos(2*tht)];

R=[1 0 0;
    0 1 0;
    0 0 2];

%Generating Qbar matrices for every ply
for i=1:n/2
    subs(T,tht,theta(i));
    Qb{i}=inv(subs(T,tht,theta(i)))*Q*R*subs(T,tht,theta(i))*inv(R);
end

for i=1:n
    Qb{n+1-i}=Qb{i};
end
fprintf('\nQbar matrix for differently oriented plies:\n\n');
for i=1:n/2
    fprintf('For theta=%d\n',theta(i)*180/pi);
    disp(Qb{i});
    fprintf('\n\n');
end

npf=zeros(n/2,1); %Number of Ply Failures

fprintf('__________________________________________________________________________\n')

for iter=1:n/2 
    fprintf('Iteration = %d\n\n',iter)

    %Generating A Matrix
    A=zeros(3);
    for i=1:n
        A=A+(Qb{i}.*(z(i+1)-z(i)));
    end
    fprintf('The A Matrix:\n\n');
    disp(A);
    
    %Generating B Matrix
    B=zeros(3);
    for i=1:n
        B=B+0.5.*(Qb{i}.*((z(i+1))^2-(z(i))^2));
    end
    fprintf('The B Matrix:\n\n');
    disp(B);
    
    
    %Generating D Matrix
    D=zeros(3);
    for i=1:n
        D=D+(Qb{i}.*((z(i+1))^3-(z(i))^3))./3;
    end
    fprintf('The D Matrix:\n\n');
    disp(D);
    
    %Forming the ABBD Matrix
    ABBD=zeros(6);
    ABBD(1:3,1:3)=A;
    ABBD(1:3,4:6)=B;
    ABBD(4:6,1:3)=B;
    ABBD(4:6,4:6)=D;
    fprintf('\nThe Final ABBD Matrix is:\n\n');
    disp(ABBD);
    
    %Load Vector
    NM=zeros(6,1);
    NM(1:3,1)=N;
    NM(4:6,1)=M;
    fprintf('\nThe Load Vector is:\n\n');
    disp(NM);
    
    %Calculating the Mid Surface Strains
    epsxyo=inv(ABBD)*NM;
    fprintf('\nThe Mid Surface Strains are:\n\n');
    fprintf('Epsxo=%f\n',epsxyo(1));
    fprintf('Epsyo=%f\n',epsxyo(2));
    fprintf('Gammaxyo=%f\n',epsxyo(3));
    fprintf('\nThe Curvatures are:\n\n');
    fprintf('Kx=%f\n',epsxyo(4));
    fprintf('Ky=%f\n',epsxyo(5));
    fprintf('Kxy=%f\n',epsxyo(6));
    
    %Caluculating epsxy for every ply
    for k=1:n
        epsxy{k}=epsxyo(1:3,1)+((z(k+1)+z(k))/2).*epsxyo(4:6,1);
    end
    
    %Calculating sigxy for every ply
    for k=1:n/2
        sigxy{k}=Qb{k}*epsxy{k};
    end
    
    %Calculating sig12 for every ply
     X=[1 0 0;
        0 1 0;
        0 0 0.5];
    for k=1:n/2
        subs(T,tht,thetaf(k));
        sig12{k}=inv(X)*subs(T,tht,thetaf(k))*X*sigxy{k};
    end  

    fprintf('\n\n');
    for k=1:n/2
        fprintf('For Ply %d & %d Stress in Material Axis are:\n',k,n+1-k);
        disp(sig12{k});
    end
    
    %Calculating Stress Ratios
    for k=1:n/2
        if sig12{k}(1,1)>=0
            SR{k}(1,1)=sig12{k}(1,1)/sig1tu;
        else
            SR{k}(1,1)=sig12{k}(1,1)/sig1cu;
        end
    
        if sig12{k}(2,1)>=0
            SR{k}(2,1)=sig12{k}(2,1)/sig2tu;
        else
            SR{k}(2,1)=sig12{k}(2,1)/sig2cu;
        end
    
        if sig12{k}(3,1)>=0
            SR{k}(3,1)=sig12{k}(3,1)/tau12u;
        else
            SR{k}(3,1)=-1*sig12{k}(3,1)/sig2cu;
        end
    end
    
    for k=1:n/2
        fprintf('\nThe Stress Ratios for %d & %d Ply are:\n',k,n+1-k);
        disp(SR{k});
    end
    
    maxvpos=zeros(n/2,3);
    for k=1:n/2
        [maxk,idxk]=max(SR{k});
        maxvpos(k,1)=k;
        maxvpos(k,2)=maxk;
        maxvpos(k,3)=idxk;
    end
    maxvpos;
    [m,i]=max(maxvpos(1:n/2,2));
    
    PFL=N(1,1)/m;
    N(1,1)=PFL;
    npf(iter,1)=PFL;
    fprintf('\nThe Ply Failure Load for Ply %d = %0.2f N/m^2\n\n',i,vpa(npf(iter,1)));
    Qb{maxvpos(i,1)}=zeros(3);
    Qb{n+1-maxvpos(i,1)}=zeros(3);
    fprintf('__________________________________________________________________________\n')
end

fprintf('\nThe First Ply Failure Load FPF = %0.2f N/m^2\n',vpa(npf(1,1)));
fprintf('\nThe Last Ply Failure Load LPF = %0.2f N/m^2\n',vpa(npf(n/2,1)));