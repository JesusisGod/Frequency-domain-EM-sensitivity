%% WRITTEN BY JIDE NOSAKARE OGUNBO in: 
%% Ogunbo, J. N., 2018, Computation of analytical sensitivity matrix for the 
%% frequency-domain EM data: MATLAB code: Computers and Geosciences


%% Conductivity vector removes air and implicitly takes sig_air=0
function dHz_dsig=ANALYTIC_SENSITIVITY_FREQ_DOMAIN_EM(rho,h,TxA,I,f)

lf=length(f);
TRS=sqrt(TxA);             % Transmitter side/length
a=TRS/sqrt(pi);            % TR=TxR, transmitter radius
    

%% Filter coefficient from Key, K., 2012. Is the fast Hankel transform faster than quadrature? 
%% Geophysics 77(3), F21-F30. doi: 10.1190/GEO2011-0237.1.
%% For usage please: Visit http://software.seg.org/disclaimer.txt

Fr=load('kk201Hankel.txt');
Filter.base=Fr(:,1);         % Filter abscissae
Filter.J0=Fr(:,2);           % Bessel function zero order
Filter.J1=Fr(:,3);           % Bessel function first order
lambda=Filter.base/a;        % % Lambda (integration variable): Ensuring discrete integration step
H1=Filter.J1;
J1=H1.';
%% 
sigma=1./rho;                  % conductivity
NL=length(sigma);
mu0=4.*pi.*1e-7;


%%
    

lla=length(lambda);

Y_j=zeros(lla,NL); 

u_j=zeros(lla,NL);
Ej=zeros(lla,NL-1);

dA0_dsigj=zeros(lla,NL-1);
dHz_dsig=zeros(lf,NL-1);
Rj=zeros(lla,NL-1);           % Length is 1 less because air's not used
    
for ifreq=1:lf
    w=2*pi*f(ifreq);

    %% Variables for each layer have length of number of layers
    iwu=1i*w*mu0;
    %% 
    for k=1:NL                  
        u_j(:,k)=sqrt(lambda.^2+iwu.*sigma(k)); % Eqn 6b
    end
        
        %%
    if NL==1 % Homogeneous
        u1=u_j(:,1);
        du1_dsig1=0.5*iwu./u1;  % Eqn 11a
        dA0_dsigj=-2*lambda./(lambda+u1).^2.*du1_dsig1;
    else
        for j=1:NL-1
            n=NL-j;
            nume=u_j(:,n)-u_j(:,n+1);                            % Numerator of eqn (7)
            deno=u_j(:,n)+u_j(:,n+1);                            % Denominator of eqn (7)
            if(sum(nume)==0), nume=1e-6*ones(size(nume)); end    % Avoid NaN from consecutive layers with equal layer resistivity in dFj_dsigj from S2nd when Rj=Fjp1=0 from nume=0
            Rj(:,n)=nume./deno;                                  % Eqn (7)

            nume2=Rj(:,n)+Y_j(:,n+1);                            % Part of numerator in eqn (9)
            deno2=1+Rj(:,n).*Y_j(:,n+1);                         % Denominator of eqn (9)
            Ej(:,n)=exp(-2*u_j(:,n)*h(n));                       % E_j in eqn (9)          
            Y_j(:,n)=nume2./deno2.*Ej(:,n);                      % Eqn (9)
        end

     %% Now for the surface h=0
     
         u1=u_j(:,1);
         R0=(lambda-u1)./(lambda+u1);                            % Eqn 8         
         Y1=Y_j(:,1);
            
%%

            uj  =u_j(:,1:end-1); % uj
            sjp1=u_j(:,2:end  ); % u(j+1)
            Yj  =Y_j(:,1:end-1); % Yj
            Yjp1=Y_j(:,2:end  ); % Y(j+1)
            hj  =repmat(h,lla,1);
            
            dA0_dR0=(1-Y1.^2)./(1+R0.*Y1).^2;                      % Eqn 13
            dA0_dY1=(1.0-R0.^2)./(1.0+R0.*Y1).^2;                  % Eqn 12
            
            du1_dsig1=0.5*iwu./u1;                                 % Eqn 14
            duj1_dsigj1=0.5*iwu./sjp1;                             % Eqn 14
            
            coef=-Yj.*(iwu)./uj;                                   % Coefficient in eqn 21
            F1st=hj;                                               % First term in eqn 21
            s2nd=sjp1./((uj+sjp1).^2.*(Rj+Yjp1));                  % Second term in eqn 21
            T3rd=(sjp1.*Yjp1)./( (uj+sjp1).^2.*(1.0+Rj.*Yjp1) );   % Third term in eqn 21            
            dYj_dsigj=coef.*(F1st-s2nd+T3rd);                      % Eqn 21
            
            dYj_dYjp1=Ej.*(1.0-Rj.^2)./(1.0+Rj.*Yjp1).^2;          % Eqn 22
            
            dYj_dRj=Ej.*(1-Yjp1.^2)./(1+Rj.*Yjp1).^2;              % Eqn 23 
            dYj_dhj=-2*uj.*Yj;                                     % Eqn 24           
            
            dR0_dsig1=-2*lambda./(lambda+u1).^2.*du1_dsig1;        % A compact form of Eqn 25 for first layer            
            dRj_dsigj1=-2*uj./(uj+sjp1).^2.*duj1_dsigj1;           % A compact form of Eqn 25 for jth layer
            
                    
            
            % NOTE YOU HAVE COLLECTED ALL VARIABLES FOR LAYER 1 TO NL-1 (OR 1:j IN THE LITERATURE)
            % REMAINING THE LAST LAYER NL, OR j+1 IN THE LITERATURE
            
            Ejj=Ej(:,end);
            ujj=uj(:,end);
            ujjp1=sjp1(:,end);
            dunp1_dsignp1=duj1_dsigj1(:,end);
            dYn_dsignp1=-2*Ejj.*ujj./(ujj+ujjp1).^2.*dunp1_dsignp1; % A compact form of Eqn 27
            
            %% Perform the product            
            
            % FIRST LAYER
            dA0_dsigj(:,1)=dA0_dY1.*dYj_dsigj(:,1)+dA0_dR0.*dR0_dsig1; % Eqn 29 (First layer Resistivity sensitivity)
            dA0_dsigj(:,NL+1)=dA0_dY1.*dYj_dhj(:,1); % Eqn 34 (First layer thickness sensitivity) 
            
            if NL>2                
                % SECOND LAYER
                dF1_dsig2b=dYj_dYjp1(:,1).*dYj_dsigj(:,2)+dYj_dRj(:,1).*dRj_dsigj1(:,1);% Bracket term in Eqn 30
                dA0_dsigj(:,2)=dA0_dY1.*dF1_dsig2b;            % Eqn 30 (2nd layer Resistivity sensitivity)

                dA0_dsigj(:,NL+2)=dA0_dY1.*dYj_dYjp1(:,1).*dYj_dhj(:,2); % Eqn 35 (Second layer thickness sensitivity)              
                
                cum_prod=1;
                for ii=3:NL-1 % THIRD LAYER TO NL-1(=jth layer in literature) LAYER
                    cum_prod=cum_prod.*dYj_dYjp1(:,ii-2);                    
                    dFj_dsigjb=dYj_dYjp1(:,ii-1).*dYj_dsigj(:,ii)+dYj_dRj(:,ii-1).*dRj_dsigj1(:,ii-1);
                    dA0_dsigj(:,ii)   =dA0_dY1.*cum_prod.*dFj_dsigjb; % Eqn 31 for 3rd layer resistivity but Eqn 32 for layer >3 resistivity sensitivity
                    dA0_dsigj(:,NL+ii)=dA0_dY1.*cum_prod.*dYj_dYjp1(:,ii-2+1).*dYj_dhj(:,ii); % Eqn 36 for 3rd to NL-1 layer thickness sensitivity
                end                
                
            end            

            %% LAST LAYER RESISTIVITY SENSITIVITY (NL=j+1 in literature)
            if NL>2
                dA0_dsigj(:,NL)=dA0_dY1.*cum_prod.*dYj_dYjp1(:,end-1).*dYn_dsignp1; % Eqn 33
            elseif(NL==2)
                dA0_dsigj(:,NL)=dA0_dY1.*dYn_dsignp1;
            end
    end
        
            %% SOLVE HANKEL TRANSFORM FOR THE SENSITIVITY MATRIX
           
            for i=1:2*NL-1                
                dHz_dsig(ifreq,i)=0.5*I*a*sum(lambda.*(dA0_dsigj(:,i)).*J1.')/a; % Eqn 28 Digital filtering of Hankel transform               
            end

end





