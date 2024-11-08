%《Milling stability prediction of variable pitch and variable helix tools using third-order updated full-discretization method》
% Updated 3rd Full-discretization method for milling SLD prediction.

% REF: “Third-order updated full-discretization method for milling stability prediction”
% REF: On the accurate calculation of milling stability limits using third-order full-discretization method
% The model parameters are from “A full-discretization method for prediction of milling stability,2010,MTM”.


%%
clear,clc;
tic
% the model parameters
N = 2;                   % number of teeth
Kt= 6e8;                % tangential cutting force coefficient (N/m^2)
Kn= 2e8;                % normal cutting force coefficient (N/m^2)
w0 = 922*2*pi;         % angular natural frequency (rad/s)
zeta = 0.011;          % relative damping (1)
m_t = 0.03993;         % modal mass (kg)
aD = 1;              % radial depth of cut ratio, for aD = 1, 0.1, and 0.05

up_or_down = -1;            % 1: up-milling, -1: down-milling
if up_or_down == 1          % up-milling
    fist= 0;                  % start angle (rad)
    fiex = acos(l-2*aD) ;  % exit angle (rad)
elseif up_or_down == -1    % down-milling 
    fist= acos(2*aD-1) ;    % start angle (rad)
    fiex = pi;                % exit angle (rad)
end
stx = 200;                  % steps of spindle speed
sty = 100;                  % steps of depth of cut
w_st = 0e-3;                % starting depth of cut (m)
w_fi = 5e-3;               % final depth of cut (m)
o_st = 2e3;                 % starting spindle speed (rpm)
o_fi = 10e3;                % final spindle speed (rpm)
% the computational parameters
m = 60;                      % number of discretization interval over one tooth pass period, i.e. one time delay T

% Discretization of the specific cutting force coefficient h(t) in Eq. (29)
for i = 1 : m+1
    dtr = 2*pi/N/m;                     % Δt
    h (i) = 0;
    for j = 1 : N                        % loop for tooth j
        fi = i*dtr +(j-1)*2*pi/N;
        if (fi >= fist)&&(fi <= fiex)
            g = 1;                        % tooth is in the cut
        else
            g = 0;                        % tooth is out of cut
        end
        h(i) = h(i) + g*(Kt*cos(fi)+Kn*sin(fi))*sin(fi);
    end
end
%-------- Begin of the proposed method --------
A0 = [-zeta * w0, 1/m_t;                      % A0 in Eq. (33)
    m_t * ((zeta* w0)^2 - w0^2), -zeta* w0];
I=eye (size (A0)) ;
invA0 = inv(A0);                               % invA0 = (A0)^(-1)

D = zeros(m+2,m+2);        % initialize matrix D, the discrete map matrix
d = ones(m+1, 1);
d(1 : 2) = 0;
D = D+diag(d,-1);
D(3, 1) = 1;

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
    for y = 1 : sty+1                               % sweeping depth of cuts
        w = w_st +(y-1)*(w_fi-w_st)/sty;         % depth of cut
        Fi= eye(m+2,m+2);  

        for i = 1 : m
            Ak0 = [0 0; -w*h(i) 0];                   
            Ak1 = [0 0; -w*h(i+1) 0];
            Bk0 = [0 0; w*h(i) 0];                   
            Bk1 = [0 0; w*h(i+1) 0];
  
            Q_minus_2=(Fi2/(3*dt)-5*Fi3/(6*dt^2)+2*Fi4/(3*dt^3)-Fi5/(6*dt^4))*Ak1+ ...
                       (Fi3/(3*dt^2)-Fi4/(2*dt^3)+Fi5/(6*dt^4))*Ak0;
            Q_minus_1=(-3*Fi2/(2*dt)+7*Fi3/(2*dt^2)-5*Fi4/(2*dt^3)+Fi5/(2*dt^4))*Ak1+ ...
                       (-3*Fi3/(2*dt^2)+2*Fi4/(dt^3)-Fi5/(2*dt^4))*Ak0;
            Q_0=(3*Fi2/(dt)-11*Fi3/(2*dt^2)+3*Fi4/(dt^3)-Fi5/(2*dt^4))*Ak1+ ...
                       (3*Fi3/(dt^2)-5*Fi4/(2*dt^3)+Fi5/(2*dt^4))*Ak0;
            Q_1=(Fi1-17*Fi2/(6*dt)+17*Fi3/(6*dt^2)-7*Fi4/(6*dt^3)+Fi5/(6*dt^4))*Ak1+ ...
                        (Fi2/dt-11*Fi3/(6*dt^2)+Fi4/(dt^3)-Fi5/(6*dt^4))*Ak0;                   
            H_0=(Fi2/(3*dt)+Fi3/(6*dt^2)-Fi4/(3*dt^3)-Fi5/(6*dt^4))*(-Bk1)+ ...
                        (Fi3/(3*dt^2)+Fi4/(2*dt^3)+Fi5/(6*dt^4))*(-Bk0);     
            H_minus_1=(Fi1-Fi2/(2*dt)-3*Fi3/(2*dt^2)+Fi4/(2*dt^3)+Fi5/(2*dt^4))*(-Bk1)+ ...
                        (Fi2/dt+Fi3/(2*dt^2)-Fi4/(dt^3)-Fi5/(2*dt^4))*(-Bk0);    
            H_minus_2=(-Fi2/(dt)+3*Fi3/(2*dt^2)-Fi5/(2*dt^4))*(-Bk1)+ ...
                        (-Fi3/(dt^2)+Fi4/(2*dt^3)+Fi5/(2*dt^4))*(-Bk0);
            H_minus_3=(Fi2/(6*dt)-Fi3/(6*dt^2)-Fi4/(6*dt^3)+Fi5/(6*dt^4))*(-Bk1)+ ...
                        (Fi3/(6*dt^2)-Fi5/(6*dt^4))*(-Bk0);

            inv0fimk1 = inv(Q_1-I);  

            D(1:2,1:2)=-inv0fimk1*(Q_0+Fi0);
            D(1:2,3)=-inv0fimk1*Q_minus_1(1:2, 1:1);
            D(1:2,4)=-inv0fimk1*Q_minus_2(1:2, 1:1);

            D(1:2,m-1)=inv0fimk1*H_minus_3(1:2, 1:1);
            D(1:2,m)=inv0fimk1*H_minus_2(1:2, 1:1);
            D(1:2,m+1)=inv0fimk1*H_minus_1(1:2, 1:1);
            D(1:2,m+2)=inv0fimk1*H_0(1:2, 1:1);

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



