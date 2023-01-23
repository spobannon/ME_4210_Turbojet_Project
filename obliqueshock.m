function [ M1,beta,ro2ro1,p2p1,T2T1,M2,p02p01,Theta,Delta_S] = obliqueshock( Flag1,Flag2,var1,var2 )
%This function solves the oblique shock wave relation based on the
%relations given in "Modern Compressible Flow With Historical Prespective-John D. Anderson Jr" 
% You can solve for different cases depending on the given data for your problem
% You can also decided whether you want the strong or weak solution
% the values assigned for the flags determine the case you want to solve for.
% variables names ending with 1 refer to the stream before deflection.
% variables names ending with 2 refer to the stream after deflection.
% Outputs: M1: Mach number before deflection.
%          M2: Mach number after deflection
%          beta: wave angle in radian
%          Theta: deflection angle in radian
%          ro2ro1: density ratio
%          p2p1: pressure ratio
%          T2T1: Temperature ratio.
%          p02p01: total pressure ratio
%          Delta_S: Change in entropy
%   Flag1=1>>var1=M1 and var2=beta
%   Flag1=2>>var1=theta and var2=beta
%   Flag1=3>>var1=M1 and var2=Theta
%   Flag2=1>>weak solution
%   Flag2=2>>strong solution
gma=1.4;
if Flag1==1
    M1=var1;
    beta=var2;
    Mn1=M1.*sin(beta);
    ro2ro1=((gma+1).*(Mn1.^2))/((gma-1).*(Mn1.^2)+2);
    p2p1=1+(2*gma)/(gma+1)*((Mn1.^2)-1);
    Mn2=sqrt(((Mn1.^2)+5)/(7*(Mn1.^2)-1));
    T2T1=p2p1/ro2ro1;
    Theta=atan(2*cot(beta).*(((M1.^2).*(sin(beta).^2)-1)/((M1.^2).*(gma+cos(2*beta))+2)));
    M2=Mn2/sin(beta-Theta);
    p02p01=(((6*(Mn1.^2))/((Mn1.^2)+5)).^3.5).*((6/(7*(Mn1.^2)-1)).^2.5);
    Delta_S=287*log(T2T1.^(3.5)/p2p1);
    
    
elseif Flag1==2
    Theta=var1;
    beta=var2;
    syms x
    s(1)=vpasolve(tan(Theta)-2*cot(beta).*(((x.^2).*(sin(beta).^2)-1)/((x.^2).*(gma+cos(2*beta))+2)), x,[0,inf]);
    M1=abs(s(1));
    Mn1=M1.*sin(beta);
    ro2ro1=((gma+1).*(Mn1.^2))/((gma-1).*(Mn1.^2)+2);
    p2p1=1+(2*gma)/(gma+1)*((Mn1.^2)-1);
    Mn2=sqrt(((Mn1.^2)+5)/(7*(Mn1.^2)-1));
    T2T1=p2p1./ro2ro1;
    M2=Mn2/sin(beta-Theta);
    p02p01=(((6*(Mn1.^2))/((Mn1.^2)+5)).^3.5).*((6/(7*(Mn1.^2)-1)).^2.5);
    Delta_S=287*log(T2T1.^(3.5)/p2p1);   
    
elseif Flag1==3
    M1=var1;
    Theta=var2;
    syms x
    s(1)=vpasolve(tan(Theta)-2*cot(x).*(((M1.^2).*(sin(x).^2)-1)/((M1.^2).*(gma+cos(2*x))+2)), x,[1.0821,pi/2]);
    s(2)=vpasolve(tan(Theta)-2*cot(x).*(((M1.^2).*(sin(x).^2)-1)/((M1.^2).*(gma+cos(2*x))+2)), x,[0,s(1)]);
    if abs(s(1))<abs(s(2))
        weaksol=abs(s(1));
        strongsol=abs(s(2));
    else
        weaksol=abs(s(2));
        strongsol=abs(s(1));
    end
    if Flag2==1
        beta=weaksol;
        Mn1=M1.*sin(beta);
        ro2ro1=((gma+1).*(Mn1.^2))/((gma-1).*(Mn1.^2)+2);
        p2p1=1+(2*gma)/(gma+1)*((Mn1.^2)-1);
        Mn2=sqrt(((Mn1.^2)+5)/(7*(Mn1.^2)-1));
        T2T1=p2p1/ro2ro1;
        M2=Mn2/sin(beta-Theta);
        p02p01=(((6*(Mn1.^2))/((Mn1.^2)+5)).^3.5).*((6/(7*(Mn1.^2)-1)).^2.5);
        Delta_S=287*log(T2T1.^(3.5)/p2p1);

    elseif Flag2==2
     beta=strongsol; 
     Mn1=M1.*sin(beta);
        ro2ro1=((gma+1).*(Mn1.^2))/((gma-1).*(Mn1.^2)+2);
        p2p1=1+(2*gma)/(gma+1)*((Mn1.^2)-1);
        Mn2=sqrt(((Mn1.^2)+5)/(7*(Mn1.^2)-1));
        T2T1=p2p1/ro2ro1;
        M2=Mn2/sin(beta-Theta);
        p02p01=(((6*(Mn1.^2))/((Mn1.^2)+5)).^3.5).*((6/(7*(Mn1.^2)-1)).^2.5);
        Delta_S=287*log(T2T1.^(3.5)/p2p1);
    else 
        warning('wrong Flag2 value,please enter 1 for weak solution or 2 for strong solution')
        
    end

    

else
    warning('please enter a suitable Flag number')
end
M2=vpa(M2);
end

