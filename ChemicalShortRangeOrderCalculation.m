clear all;
clc;

Li=1; Be=2; B=3; C=4; N=5;  Na=6; Mg=7; Al=8; Si=9; P=10; Ca=11; Sc=12; Ti=13; V=14; 
Cr=15; Mn=16; Fe=17; Co=18; Ni=19; Cu=20; Zn=21; Ga=22; Ge=23; Sr=24; Y=25; Zr=26;
Nb=27; Mo=28; Tc=29; Rh=30; Pd=31; Ag=32; Cd=33; La=34; Ce=35; Nd=36; Eu=37; Gd=38;
Tb=39; Dy=40; Er=41; Yb=42; Hf=43; Ta=44; W=45; Pt=46; Au=47; Pb=48; Sn=49; Ho=50; Tm=51;
Lu=52; Re=53; Ru=54; Bi=55; In=56; Sb=57; Os=58;

%read the atomic size in the database
enthalpy = xlsread('mixingenthalpy.xlsx','B2:BG59'); 

%constituent elements 
ele=[Co Cr Fe Ni];

%chemical composition in atomic ratio
com=[1 1 1 1];

%convert chemical composition in atomic percentage
for i=1:length(ele)
    c(i)=com(i)/sum(com);
end

N1=length(ele);

% gas constant R J/mol/K
R=8.3144598;
%Temperature K
T1=5000:-50:50;
T1=T1';

%P_ij is the number of i-j bond p_ij=P_ij/zN is the probability of i-j bond
p0=zeros(N1);
for i=1:N1-1
    for j=i+1:N1
        p0(i,j)=c(i)*c(j);
        p0(j,i)=p0(i,j);
    end
end

for i=1:N1
    p0(i,i)=c(i)/2-(sum(p0(i,:))-p0(i,i))/2;
end

%Phi_th is the calculated ordering parameter
%Pij_th is the calculated atomic bond number
%Smix is the calculated mixing entropy

%initialization column vector x0
for i=1:N1-1
    for j=i+1:N1
        x0(convert(N1,i,j),1)=p0(i,j);
    end
end

%Enlthpy and Omega values KJ/mol
Omg=zeros(length(ele));
Hmix=zeros(length(ele));
for i=1:length(ele)
    for j=i+1:length(ele)
        Hmix(i,j)=enthalpy(ele(i),ele(j));
        if Hmix(i,j)==0
            Hmix(i,j)=0;
        end
        Omg(i,j)=4*Hmix(i,j);
        Hmix(j,i)=Hmix(i,j);
        Omg(j,i)=Omg(i,j);
    end
end

%loop for temperature
for ii=1:length(T1)
    T=T1(ii);
    % calculate A(i,j)
    for i=1:length(ele)
        for j=i+1:length(ele)
            A(i,j)=exp(Omg(i,j)*1000/R/T);
            A(j,i)=A(i,j);
        end
    end
    if ii==1
        % initiation Broyden matrix 
        B0=eye(length(x0));
        % initiation function value
        F0=Function(x0,N1,c,A);

        % start iteration
        xk=x0; Fk=F0; Bk=B0;
    elseif ii>1
        % initiation Broyden matrix 
        B0=eye(length(x0));
        % initiation function value
        Fk=Function(xk,N1,c,A);
        Bk=B0;
    end

    
    tol=1e-4;
    N=1e4;
    for m=1:N
        xk1=xk-Bk*Fk;
        Fk1=Function(xk1,N1,c,A);
        sk=xk1-xk;
        yk=Fk1-Fk;
        Bk1=Bk+(sk-Bk*yk)*sk'*Bk/(sk'*Bk*yk);
        % sulotion judgement
        for i=1:length(x0)
            ere(i,1)=abs(xk1(i)-xk(i))/xk(i);
        end

        if max(max(ere))<tol
            disp('Congratulation. Convergence!!!');
            disp(['Number of Iterations:',num2str(m)]);
            if min(xk1)<0
                disp('Minus value occurs!!!');
            end
            % calculation ordering parameter phi_ij
            for i=1:length(x0)
                temp=convert2(N1,i);
                Phi_th(ii,i)=(xk1(i)-c(temp(1))*c(temp(2)))/(min(c(temp(1)),c(temp(2)))-c(temp(1))*c(temp(2)));
                Pij_th(ii,:)=xk1';
            end
            % calculation mixing entropy
            p1=zeros(N1);
            for i=1:length(x0)
                temp=convert2(N1,i);
                p1(temp(1),temp(2))=xk1(i);
                p1(temp(2),temp(1))=p1(temp(1),temp(2));
            end
            for i=1:N1
                p1(i,i)=c(i)/2-(sum(p1(i,:))-p1(i,i))/2;
            end
       
            Sum1=0;
            Sum2=0;
            for i=1:N1
                Sum1=Sum1-2*p1(i,i)*log(2*p1(i,i)/c(i));
            end
            for i=1:N1
                for j=1:N1
                    if i~=j
                        Sum2=Sum2-p1(j,i)*log(p1(j,i)/c(i));
                    end
                end
            end
            Smix(ii,1)=Sum1+Sum2;  
            break;
        end
        if m==N
            disp('Non Convergence!');
            return;
        end 
        xk=xk1;
        Bk=Bk1;
        Fk=Fk1;
    end% end of iteration
    
end
           
            
