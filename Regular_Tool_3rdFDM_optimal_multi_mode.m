%《Milling stability prediction of variable pitch and variable helix tools using third-order updated full-discretization method》
% Updated 3rd Full-discretization method for milling SLD prediction.

% REF: “Third-order updated full-discretization method for milling stability prediction”
% REF: On the accurate calculation of milling stability limits using third-order full-discretization method
% The model parameters are from “A full-discretization method for prediction of milling stability,2010,MTM”.

% If only considering 1 mode in x direction and 1 mode in y direction,
% matrix "G" equals to unit matrix I, and matrix "G" can be deleted.
% In this way, the calculation times can be significantly reduced.
% see "Regular_Tool_3rdFDM_single_mode.m"

%%
close all; clear all; clc;
tic
% the model parameters
N = 4;                   % number of teeth
Kt= 6.97e8;                % tangential cutting force coefficient (N/m^2)
Kn= 2.558e8;                % normal cutting force coefficient (N/m^2)
aD = 0.5;              % radial depth of cut ratio, for aD = 1, 0.1, and 0.05


% 2-DOF
wx1=441.64*2*pi;  wx2=563.6*2*pi;  wx3=778.56*2*pi;       wy1=516.2*2*pi;
mx1=11.125;       mx2=1.4986;      mx3=13.063;            my1=1.199;
zetax1=0.028722;  zetax2=0.055801; zetax3=0.058996;       zetay1=0.025004;

num_modal_x=3;    num_modal_y=1;    % number of modals 
num_modal=num_modal_x+num_modal_y;

G=[ones(1,num_modal_x),zeros(1,num_modal_y);
   zeros(1,num_modal_x),ones(1,num_modal_y)]; % for multi-mode using

% 2-DOF
M=diag([mx1,mx2,mx3,my1]);% mass matrix
invM=inv(M);
C=diag([2*mx1*zetax1*wx1,2*mx2*zetax2*wx2,2*mx3*zetax3*wx3,2*my1*zetay1*wy1]);% damping matrix
K=diag([mx1*wx1^2,mx2*wx2^2,mx3*wx3^2,my1*wy1^2]);% stiffness matrix


up_or_down = -1;            % 1: up-milling, -1: down-milling
if up_or_down == 1          % up-milling
    fist= 0;                  % start angle (rad)
    fiex = acos(1-2*aD) ;  % exit angle (rad)
elseif up_or_down == -1    % down-milling 
    fist= acos(2*aD-1) ;    % start angle (rad)
    fiex = pi;                % exit angle (rad)
end

% the simulation parameters
stx = 200;                  % steps of spindle speed
sty = 100;                  % steps of depth of cut
w_st = 0e-3;                % starting depth of cut (m)
w_fi = 10e-3;               % final depth of cut (m)
o_st = 2e3;                 % starting spindle speed (rpm)
o_fi = 10e3;                % final spindle speed (rpm)
% the computational parameters
m = 72;                      % number of discretization interval over one tooth pass period, i.e. one time delay T

    
D = zeros(num_modal*(m+2),num_modal*(m+2)); % initialize matrix D, the discrete map matrix
D(2*num_modal+1:3*num_modal,1:num_modal)=eye(num_modal);
D(3*num_modal+1:end,(2*num_modal+1):(end-num_modal))=eye(num_modal*(m-1));


hxx=zeros(m+1,1);       % initialize h
hxy=zeros(m+1,1);
hyx=zeros(m+1,1);
hyy=zeros(m+1,1);

% Discretization of the specific cutting force coefficient h(t) in Eq. (29)
for i = 1 : m+1
    dtr = 2*pi/N/m;                     % Δt
    for j = 1 : N                        % loop for tooth j
        fi = i*dtr +(j-1)*2*pi/N;
        if (fi >= fist)&&(fi <= fiex)
            g = 1;                        % tooth is in the cut
        else
            g = 0;                        % tooth is out of cut
        end
        hxx(i) = hxx(i) + g*(Kt*cos(fi)+Kn*sin(fi))*sin(fi);
        hxy(i) = hxy(i) + g*(Kt*cos(fi)+Kn*sin(fi))*cos(fi);
        hyx(i) = hyx(i) + g*(-Kt*sin(fi)+Kn*cos(fi))*sin(fi);
        hyy(i) = hyy(i) + g*(-Kt*sin(fi)+Kn*cos(fi))*cos(fi);
    end
end
%-------- Begin of the proposed method --------
A0=[-invM*C/2, invM;
    C*invM*C/4-K,-C*invM/2];
I=eye (size (A0)) ;
invA0 = inv(A0);                               % invA0 = (A0)^(-1)


% start of computation
for x = 1 : stx+1                               % sweeping spindle speeds
    o = o_st +(x-1)*(o_fi-o_st)/stx;         % spindle speed
    tau= 60/o/N;                                 % time delay
    dt = tau/m;                                  % time step
    %------- Calculation of Φ0,Φ1,Φ2,Φ3,Φ4,Φ5 -------
    Fi0 = expm(A0*dt);                          % Φ0
    Fi1=  invA0 * (Fi0 - I);                     % Φ1
    Fi2 = invA0 * (Fi0 * dt - Fi1);             % Φ2
    Fi3 = invA0 * (Fi0 * dt^2 - 2 * Fi2) ;      % Φ3
    Fi4 = invA0 * (Fi0 * dt^3 - 3 * Fi3) ;      % Φ4
    Fi5 = invA0 * (Fi0 * dt^4 - 4 * Fi4) ;      % Φ5
    %----- The end of calculation of Φ0,Φ1,Φ2,Φ3,Φ4,Φ5  -----

    %------- Calculation of coefficients -------
    Q_minus_2_1k=(Fi2/(3*dt)-5*Fi3/(6*dt^2)+2*Fi4/(3*dt^3)-Fi5/(6*dt^4));
    Q_minus_2_0k=(Fi3/(3*dt^2)-Fi4/(2*dt^3)+Fi5/(6*dt^4));
    Q_minus_1_1k=(-3*Fi2/(2*dt)+7*Fi3/(2*dt^2)-5*Fi4/(2*dt^3)+Fi5/(2*dt^4));
    Q_minus_1_0k=(-3*Fi3/(2*dt^2)+2*Fi4/(dt^3)-Fi5/(2*dt^4));
    Q_0_1k=(3*Fi2/(dt)-11*Fi3/(2*dt^2)+3*Fi4/(dt^3)-Fi5/(2*dt^4));
    Q_0_0k=(3*Fi3/(dt^2)-5*Fi4/(2*dt^3)+Fi5/(2*dt^4));
    Q_1_1k=(Fi1-17*Fi2/(6*dt)+17*Fi3/(6*dt^2)-7*Fi4/(6*dt^3)+Fi5/(6*dt^4));
    Q_1_0k=(Fi2/dt-11*Fi3/(6*dt^2)+Fi4/(dt^3)-Fi5/(6*dt^4));

    H_0_1j=(Fi2/(3*dt)+Fi3/(6*dt^2)-Fi4/(3*dt^3)-Fi5/(6*dt^4));
    H_0_0j=(Fi3/(3*dt^2)+Fi4/(2*dt^3)+Fi5/(6*dt^4));
    H_minus_1_1j=(Fi1-Fi2/(2*dt)-3*Fi3/(2*dt^2)+Fi4/(2*dt^3)+Fi5/(2*dt^4));
    H_minus_1_0j=(Fi2/dt+Fi3/(2*dt^2)-Fi4/(dt^3)-Fi5/(2*dt^4));
    H_minus_2_1j=(-Fi2/(dt)+3*Fi3/(2*dt^2)-Fi5/(2*dt^4));
    H_minus_2_0j=(-Fi3/(dt^2)+Fi4/(2*dt^3)+Fi5/(2*dt^4));
    H_minus_3_1j=(Fi2/(6*dt)-Fi3/(6*dt^2)-Fi4/(6*dt^3)+Fi5/(6*dt^4));
    H_minus_3_0j=(Fi3/(6*dt^2)-Fi5/(6*dt^4));
    %------- The end of Calculation of coefficients -------

    for y = 1 : sty+1                               % sweeping depth of cuts
        w = w_st +(y-1)*(w_fi-w_st)/sty;         % depth of cut
        Fi= eye(num_modal*(m+2),num_modal*(m+2));  % construct transition matrix Fi

        for i = 1 : m
              
            hi0=w*[hxx(i) hxy(i) 
                hyx(i) hyy(i)];
            hi0 = G'* hi0 * G;
            
            hi1=w*[hxx(i+1) hxy(i+1) 
                 hyx(i+1) hyy(i+1)];
            hi1 = G'* hi1 * G;
            
            A_zeros=zeros(num_modal,num_modal);
            
            Ak0=[A_zeros,A_zeros
                 -hi0,A_zeros];
            Ak1=[A_zeros,A_zeros
                 -hi1,A_zeros];
            Bk0=-Ak0;     
            Bk1=-Ak1;     

            Q_minus_2=Q_minus_2_1k*Ak1+ Q_minus_2_0k*Ak0;
            Q_minus_1=Q_minus_1_1k*Ak1+ Q_minus_1_0k*Ak0;
            Q_0=Q_0_1k*Ak1+ Q_0_0k*Ak0;
            Q_1=Q_1_1k*Ak1+Q_1_0k*Ak0; 
            H_0=H_0_1j*(-Bk1)+ H_0_0j*(-Bk0);     
            H_minus_1=H_minus_1_1j*(-Bk1)+ H_minus_1_0j*(-Bk0);    
            H_minus_2=H_minus_2_1j*(-Bk1)+H_minus_2_0j*(-Bk0);
            H_minus_3=H_minus_3_1j*(-Bk1)+H_minus_3_0j*(-Bk0);


            inv0fimk1 = inv(Q_1-I);  

            D(1:2*num_modal,1:1:2*num_modal)=-inv0fimk1*(Q_0+Fi0);
            D(1:2*num_modal,2*num_modal+1:3*num_modal)=-inv0fimk1*Q_minus_1(:, 1:num_modal);
            D(1:2*num_modal,3*num_modal+1:4*num_modal)=-inv0fimk1*Q_minus_2(:, 1:num_modal);

            D(1:2*num_modal,(num_modal*(m-2)+1):(num_modal*(m-1)))=inv0fimk1*H_minus_3(:, 1:num_modal);
            D(1:2*num_modal,(num_modal*(m-1)+1):(num_modal*m))=inv0fimk1*H_minus_2(:, 1:num_modal);
            D(1:2*num_modal,(num_modal*(m)+1):(num_modal*(m+1)))=inv0fimk1*H_minus_1(:, 1:num_modal);
            D(1:2*num_modal,(num_modal*(m+1)+1):(num_modal*(m+2)))=inv0fimk1*H_0(:, 1:num_modal);

            
            Fi = D * Fi;   

        end
    ss(x,y) = o;                         % matrix of spindle speeds
    dc(x,y) = w;                         % matrix of depth of cuts
    ei(x,y) = max(abs(eig(Fi)));      % matrix of eigenvalues        

    end                                   % End of sweeping depth of cuts
   stx+1-x
end                                       % End of sweeping spindle speeds
%-------- The end of the proposed method --------
toc
figure;
contour (ss, dc, ei, [1, 1], 'k' );   % moduluses of all the eigenvalues of the transition matrix Φ
xlabel ('Spindle speed (rpm)' );
ylabel ('Axial depth of cut (m)' );