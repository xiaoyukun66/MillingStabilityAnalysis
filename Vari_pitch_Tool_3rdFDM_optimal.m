%此文件为Vari_pitch_Tool_3rdFDM.m的优化文件，以提升效率。
%(1)提前计算各表达式系数  将inv改为左除；ss,dc,ei初始化;    
%(1) 383.72s 降低至 344.57s
%(2)系数变化项与系数不变项的计算，只改变矩阵对应位置的值 286s
%(3)系数不变项的计算,以"先累加再矩阵乘法"替换"先矩阵乘法再累加"。 265s

%考虑螺旋角效应的变齿距立铣颤振稳定性识别,多模态

%%
clear , clc;
tic
%模型参数
beta=30/180*pi;                         %刀具螺旋角
R=19.05/2*10^(-3);                      %刀具半径m
N = 4;                                  %刀齿数
angle_chiju=[70,110,70,110]*pi/180;     %变齿距,rad

%以第一个齿为参考，计算各齿滞后角(用于转动角phi的计算)
angle_zhihou=zeros(1,N);
for j=2:N
    angle_zhihou(j)=sum(angle_chiju(1:j-1));
end

Kt= 6.97e8;                             % 切向切削系数(N/m^2)
Kn= 2.558e8;                            % 法向切削系数(N/m^2)
aD = 0.5;              % radial depth of cut ratio, for aD = 1, 0.1, and 0.05

% 2-DOF
wx1=441.64*2*pi;  wx2=563.6*2*pi;  wx3=778.56*2*pi;       wy1=516.2*2*pi;
mx1=11.125;       mx2=1.4986;      mx3=13.063;            my1=1.199;
zetax1=0.028722;  zetax2=0.055801; zetax3=0.058996;       zetay1=0.025004;

num_modal_x=3;    num_modal_y=1;    % number of modals 
num_modal=num_modal_x+num_modal_y;

G=[ones(1,num_modal_x),zeros(1,num_modal_y);
   zeros(1,num_modal_x),ones(1,num_modal_y)]; % for multi-mode using

%先考虑2自由度多阶模态
M=diag([mx1,mx2,mx3,my1]);%质量矩阵
invM=inv(M);
C=diag([2*mx1*zetax1*wx1,2*mx2*zetax2*wx2,2*mx3*zetax3*wx3,2*my1*zetay1*wy1]);%阻尼矩阵
K=diag([mx1*wx1^2,mx2*wx2^2,mx3*wx3^2,my1*wy1^2]);%刚度矩阵


%顺铣，计算切入角和切出角
up_or_down = -1;            % 1: up-milling, -1: down-milling
if up_or_down == 1          % up-milling
    fist= 0;                  % start angle (rad)
    fiex = acos(1-2*aD) ;  % exit angle (rad)
elseif up_or_down == -1    % down-milling 
    fist= acos(2*aD-1) ;    % start angle (rad)
    fiex = pi;                % exit angle (rad)
end

stx = 200;                  % steps of spindle speed
sty = 100;                  % steps of depth of cut
w_st = 0e-3;                % starting depth of cut (m)
w_fi = 8e-3;               % final depth of cut (m)
o_st = 2e3;                 % starting spindle speed (rpm)
o_fi = 10e3;                % final spindle speed (rpm)
% the computational parameters
m = 72;                      % number of discretization interval over one tooth pass period, i.e. one time delay T


%考虑螺旋效应下，变齿距立铣刀的切削力，与时间、刀齿、切厚均相关
%初始化瞬时切厚
kxx=zeros(stx+1,sty+1,m+1,N);                       
kxy=zeros(stx+1,sty+1,m+1,N);   
kyx=zeros(stx+1,sty+1,m+1,N);   
kyy=zeros(stx+1,sty+1,m+1,N);   
%以d_w离散，计算此时的d_phi
d_w=(w_fi-w_st)/sty;
d_phi=d_w*tan(beta)/R;
%初始化径向瞬时接触角
max_NA=floor(w_fi/d_w);

% 给定转速、切深、时间，计算周期系数项
for x = 1 : stx+1                            % sweeping spindle speeds
    o = o_st +(x-1)*(o_fi-o_st)/stx;         % spindle speed
    dtr = 60/(o*m);                         %时间微元
    for y = 1 : sty+1
        w = w_st +(y-1)*(w_fi-w_st)/sty;        %此时的切削深度
        NA=floor(w/d_w);
        res_w=w-NA*d_w;
        for i=1:m+1                             %时间循环
            fi=zeros(N,NA+1);                  %初始化切触角
            for j=1:N                           %刀齿循环(第一个刀齿无滞后角；后续刀齿滞后角为各齿间角之和)
                %首先获取切触范围
                for k=1:NA+1                       %切深微元循环
                    if k==NA+1
                        fi(j,k) = 2*pi*o/60*i*dtr+angle_zhihou(j) -((k-1)*d_phi+res_w*tan(beta)/R);%计算切触角
                    else
                        fi(j,k) = 2*pi*o/60*i*dtr+angle_zhihou(j) -(k-1)*d_phi;%计算切触角
                    end
                end
                fi=mod(fi,2*pi);%超过2pi的要取余

                %然后进行切触判断,并累加切厚
                for k=1:NA+1
                    if (fi(j,k) >= fist)&&(fi(j,k) <= fiex)
                        g = 1;                        % tooth is in the cut
                    else
                        g = 0;                        % tooth is out of cut
                    end
                    if k==NA+1
                        kxx(x,y,i,j) = kxx(x,y,i,j) + g*res_w*(Kt*cos(fi(j,k))+Kn*sin(fi(j,k)))*sin(fi(j,k));
                        kxy(x,y,i,j) = kxy(x,y,i,j) + g*res_w*(Kt*cos(fi(j,k))+Kn*sin(fi(j,k)))*cos(fi(j,k));
                        kyx(x,y,i,j) = kyx(x,y,i,j) + g*res_w*(-Kt*sin(fi(j,k))+Kn*cos(fi(j,k)))*sin(fi(j,k));
                        kyy(x,y,i,j) = kyy(x,y,i,j) + g*res_w*(-Kt*sin(fi(j,k))+Kn*cos(fi(j,k)))*cos(fi(j,k)); 
                    else
                        kxx(x,y,i,j) = kxx(x,y,i,j) + g*d_w*(Kt*cos(fi(j,k))+Kn*sin(fi(j,k)))*sin(fi(j,k));
                        kxy(x,y,i,j) = kxy(x,y,i,j) + g*d_w*(Kt*cos(fi(j,k))+Kn*sin(fi(j,k)))*cos(fi(j,k));
                        kyx(x,y,i,j) = kyx(x,y,i,j) + g*d_w*(-Kt*sin(fi(j,k))+Kn*cos(fi(j,k)))*sin(fi(j,k));
                        kyy(x,y,i,j) = kyy(x,y,i,j) + g*d_w*(-Kt*sin(fi(j,k))+Kn*cos(fi(j,k)))*cos(fi(j,k));
                    end

                end
            end
        end
    end
    fprintf('(1/2)瞬时切厚计算中：%d\n',stx+1-x)
end

%---------------瞬时切厚已通过等齿距方法验证，结果一致------------%

%-------- Begin of the proposed method --------
%-------- 时滞项、状态项均采用线性插值
A0=[-invM*C/2, invM;
    C*invM*C/4-K,-C*invM/2];
I=eye (size (A0)) ;
invA0 = inv(A0);                                % invA0 = (A0)^(-1)

%初始化kjp
Kjp=zeros(2,2);
Kjp1=zeros(2,2);

% %初始化A_zeros,B_zeros
% A_zeros=zeros(num_modal,num_modal);
% B_zeros=zeros(num_modal,num_modal);

%初始化ss,dc,ei
ss=zeros(stx+1,sty+1);                         
dc=zeros(stx+1,sty+1);
ei=zeros(stx+1,sty+1);

% start of computation
for x = 1 : stx+1                               % sweeping spindle speeds
    o = o_st +(x-1)*(o_fi-o_st)/stx;            % spindle speed
    tau= 60/o;                                  % 周期T
    dt = tau/m;                                 % time step

    tau_zhihou=angle_chiju/(2*pi)*tau;        %每个齿的时滞
    mj_zhihou=floor(tau_zhihou/dt+0.5);       %时滞对应的mj，取整

    mj_zhihou=circshift(mj_zhihou,[1,1]);       %时间滞后为j~j+1齿间的间距
    mj_zhihou_max=max(mj_zhihou);               %时间滞后mj的最大值

    D = zeros(num_modal*(mj_zhihou_max+2),num_modal*(mj_zhihou_max+2)); % initialize matrix D, the discrete map matrix
    D(2*num_modal+1:3*num_modal,1:num_modal)=eye(num_modal);
    D(3*num_modal+1:end,(2*num_modal+1):(end-num_modal))=eye(num_modal*(mj_zhihou_max-1));

    %------- Calculation of Φ0,Φ1,Φ2,Φ3,Φ4,Φ5 -------
    Fi0 = expm(A0*dt);                          % Φ0
    Fi1=  A0 \ (Fi0 - I);                     % Φ1
    Fi2 = A0 \ (Fi0 * dt - Fi1);             % Φ2
    Fi3 = A0 \ (Fi0 * dt^2 - 2 * Fi2) ;      % Φ3
    Fi4 = A0 \ (Fi0 * dt^3 - 3 * Fi3) ;      % Φ4
    Fi5 = A0 \ (Fi0 * dt^4 - 4 * Fi4) ;      % Φ5
    %----- The end of calculation of Φ0,Φ1,Φ2,Φ3,Φ4,Φ5  -----

    %------- 提前计算各表达式系数 -------
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
    %------- 各表达式系数计算结束 -------

    for y = 1 : sty+1                               % sweeping depth of cuts
        w = w_st +(y-1)*(w_fi-w_st)/sty;         % depth of cut
        Fi= eye(num_modal*(mj_zhihou_max+2),num_modal*(mj_zhihou_max+2));  % construct transition matrix Fi

        for i = 1 : m
            A0k=zeros(2*num_modal,2*num_modal);
            A1k=zeros(2*num_modal,2*num_modal);
            %获取系数不变项
            temp_kxx_i0=sum(kxx(x,y,i,:)); 
            temp_kxy_i0=sum(kxy(x,y,i,:)); 
            temp_kyx_i0=sum(kyx(x,y,i,:)); 
            temp_kyy_i0=sum(kyy(x,y,i,:)); 
            Kjp=[temp_kxx_i0 temp_kxy_i0
                 temp_kyx_i0 temp_kyy_i0];
            Kjp_mul0=G'*Kjp*G;
            
            temp_kxx_i1=sum(kxx(x,y,i+1,:)); 
            temp_kxy_i1=sum(kxy(x,y,i+1,:)); 
            temp_kyx_i1=sum(kyx(x,y,i+1,:)); 
            temp_kyy_i1=sum(kyy(x,y,i+1,:)); 
            Kjp1=[temp_kxx_i1 temp_kxy_i1
                 temp_kyx_i1 temp_kyy_i1];
            Kjp_mul1=G'*Kjp1*G;

            A0k(num_modal+1:end,1:num_modal)=-Kjp_mul0;
            A1k(num_modal+1:end,1:num_modal)=-Kjp_mul1;



            Q_minus_2=Q_minus_2_1k*A1k+ Q_minus_2_0k*A0k;
            Q_minus_1=Q_minus_1_1k*A1k+ Q_minus_1_0k*A0k;
            Q_0=Q_0_1k*A1k+ Q_0_0k*A0k;
            Q_1=Q_1_1k*A1k+Q_1_0k*A0k;  

            
            Di_stable=zeros(2*num_modal,num_modal*(mj_zhihou_max+2));
            Di_stable(:,1:4*num_modal) = ...
            [
                -(Q_1-I)\(Q_0+Fi0), ...
                -(Q_1-I)\Q_minus_1(:, 1:num_modal), ...
                -(Q_1-I)\Q_minus_2(:, 1:num_modal)
            ];


            %获取系数随刀齿变化项
            Dj=zeros(2*num_modal,num_modal*(mj_zhihou_max+2));
            Dj_temp=zeros(2*num_modal,num_modal*(mj_zhihou_max+2));
            B0j=zeros(2*num_modal,2*num_modal);
            B1j=zeros(2*num_modal,2*num_modal);
            for j=1:N
                Kjp=[kxx(x,y,i,j) kxy(x,y,i,j) 
                    kyx(x,y,i,j) kyy(x,y,i,j)];
                Kjp_mul0=G'*Kjp*G;

                Kjp1=[kxx(x,y,i+1,j) kxy(x,y,i+1,j) 
                      kyx(x,y,i+1,j) kyy(x,y,i+1,j)];
                Kjp_mul1=G'*Kjp1*G;

                B0j(num_modal+1:end,1:num_modal)=Kjp_mul0;
                B1j(num_modal+1:end,1:num_modal)=Kjp_mul1;


                H_0=H_0_1j*(-B1j)+ H_0_0j*(-B0j);     
                H_minus_1=H_minus_1_1j*(-B1j)+ H_minus_1_0j*(-B0j);    
                H_minus_2=H_minus_2_1j*(-B1j)+H_minus_2_0j*(-B0j);
                H_minus_3=H_minus_3_1j*(-B1j)+H_minus_3_0j*(-B0j);

                
                Dj_temp(:,(num_modal*(mj_zhihou(j)-2)+1):(num_modal*(mj_zhihou(j)+2))) = ...
                [
                 (Q_1-I)\H_minus_3(:, 1:num_modal), ...
                 (Q_1-I)\H_minus_2(:, 1:num_modal), ...
                 (Q_1-I)\H_minus_1(:, 1:num_modal), ...
                 (Q_1-I)\H_0(:, 1:num_modal)
                ];

 
                %Dj累加，只累加首行以缩减时间。
                Dj(:,   (num_modal*(mj_zhihou(j)-2)+1)  :  (num_modal*(mj_zhihou(j)+2)))= ...
                    Dj(:,   (num_modal*(mj_zhihou(j)-2)+1)  :  (num_modal*(mj_zhihou(j)+2)))+ ...
                    Dj_temp(:, (num_modal*(mj_zhihou(j)-2)+1)  :  (num_modal*(mj_zhihou(j)+2)));

%                 Dj=Dj_temp+Dj;
            end
            Di_add=D;
            Di_add(1:2*num_modal,:)=Di_stable+Dj;
            Fi = Di_add * Fi;                                                % Eq. (27)
        end
    ss(x,y) = o;                         % matrix of spindle speeds
    dc(x,y) = w;                         % matrix of depth of cuts
    ei(x,y) = max(abs(eig(Fi)));      % matrix of eigenvalues
    end                                   % End of sweeping depth of cuts
    fprintf('(2/2)求解计算中：%d\n',stx+1-x)
end                                       % End of sweeping spindle speeds
%-------- The end of the proposed method --------
toc
figure;
contour (ss, dc, ei, [1, 1], 'k' );   % moduluses of all the eigenvalues of the transition matrix Φ
xlabel ('Spindle speed (rpm)' );
ylabel ('Axial depth of cut (m)' );


