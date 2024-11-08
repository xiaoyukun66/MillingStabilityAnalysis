%% 变螺旋-变齿距立铣颤振稳定性识别
%%
clear, clc;
tic
%模型参数                     
beta=[30,35,30,35]./180*pi;         %刀具螺旋角
N = 4;                              %刀齿数
R=19.05/2*10^(-3);                  %刀具半径m

%刀尖点处的齿距角；变齿距下，齿距角会随刀具高度变化
angle_chiju_0=[80,100,80,100]*pi/180;  

%以第一个齿为参考，计算各齿刀尖滞后角(用于转动角phi的计算)
angle_zhihou_0=zeros(1,N);
for j=2:N
    angle_zhihou_0(j)=sum(angle_chiju_0(1:j-1));
end


Kt= 6.97e8;                             % 切向切削系数(N/m^2)
Kn= 2.558e8;                            % 法向切削系数(N/m^2)
aD = 0.5;              % radial depth of cut ratio, for aD = 1, 0.1, and 0.05


%先考虑2自由度多阶模态
wx1=441.64*2*pi;  wx2=563.6*2*pi;  wx3=778.56*2*pi;       wy1=516.21*2*pi;
mx1=11.125;       mx2=1.4986;      mx3=13.063;            my1=1.199;
zetax1=0.028722;  zetax2=0.055801; zetax3=0.058996;       zetay1=0.025004;

num_modal_x=3;    num_modal_y=1;    % number of modals 
num_modal=num_modal_x+num_modal_y;

G=[ones(1,num_modal_x),zeros(1,num_modal_y);
   zeros(1,num_modal_x),ones(1,num_modal_y)]; % for multi-mode using

M=diag([mx1,mx2,mx3,my1]);%质量矩阵
invM=inv(M);
C=diag([2*mx1*zetax1*wx1,2*mx2*zetax2*wx2,2*mx3*zetax3*wx3,2*my1*zetay1*wy1]);%阻尼矩阵
K=diag([mx1*wx1^2,mx2*wx2^2,mx3*wx3^2,my1*wy1^2]);%刚度矩阵


%逆铣，计算切入角和切出角
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
m = 60;                      % number of discretization interval over one tooth pass period, i.e. one time delay T

d_w=(w_fi-w_st)/sty;        %刀具轴向离散微元
d_phi=d_w*tan(beta)/R;      %计算刀具轴向微元dw对应的d_phi

%变螺旋变齿距立铣刀的切削力，与时间、刀齿、切厚均相关
%初始化瞬时切厚
max_NA=floor(w_fi/d_w);
kxx=zeros(stx+1,sty+1,m+1,N,max_NA+1);                       
kxy=zeros(stx+1,sty+1,m+1,N,max_NA+1);   
kyx=zeros(stx+1,sty+1,m+1,N,max_NA+1);   
kyy=zeros(stx+1,sty+1,m+1,N,max_NA+1);

%初始化时滞对应的mjl
mjl_zhihou=zeros(stx+1,sty+1,m+1,N,max_NA+1);

% 给定转速、切深、时间，高度，计算径向瞬时接触角、齿距角
for x = 1 : stx+1                            % sweeping spindle speeds
    o = o_st +(x-1)*(o_fi-o_st)/stx;          %主轴转速
    T=60/o;                                   %主轴旋转周期
    dtr = T/m;                                %时间微元
    for y = 1 : sty+1
        w = w_st +(y-1)*(w_fi-w_st)/sty;        %此时的切削深度
        %以d_w离散刀具轴向高度
        NA=floor(w/d_w);
        res_w=w-NA*d_w;
        %计算瞬时切厚
        for i=1:m+1                             %时间循环
            fi=zeros(N,NA+1);                  %初始化切触角
            angle_chiju=zeros(N,NA+1);         %初始化齿距角
            %首先获取切触范围
            for j=1:N                          
                for k=1:NA+1                       %切深微元循环
                    if k==NA+1
                        fi(j,k) = 2*pi*o/60*i*dtr + angle_zhihou_0(j) - ((k-1)*d_phi(j)+res_w*tan(beta(j))/R);%计算切触角
                    else
                        fi(j,k) = 2*pi*o/60*i*dtr + angle_zhihou_0(j) -(k-1)*d_phi(j);%计算切触角
                    end
                end
            end

            %此时拥有刀具径向切触角fi_j,k，那么可以计算不同高度下，j与j+1刃之间的齿距角
            for j=1:N
                if j==1
                    angle_chiju(j,:)=2*pi+fi(j,:)-fi(N,:);
                else
                    angle_chiju(j,:)=fi(j,:)-fi(j-1,:);
                end
            end
            angle_chiju=mod(angle_chiju,2*pi);%齿距角取余，齿1与齿N的齿距角可能为2pi

            fi=mod(fi,2*pi);%为后续判断切触范围，fi超过2pi的要取余

            tau_zhihou=angle_chiju/(2*pi)*T;%齿j在高度k的时滞
            mjl_zhihou(x,y,i,:,1:(NA+1))=floor(tau_zhihou/dtr+0.5);%时滞对应的mj，取整
            mjl_zhihou(x,y,i,:,1:(NA+1))=circshift(mjl_zhihou(x,y,i,:,1:(NA+1)),[1,1]); %刀齿j的时间滞后为j~j+1齿间的间距

            %计算不同转速、切深、时间、刀刃、高度下的切厚
            for j=1:N
                for k=1:NA+1
                    %切触判断
                    if (fi(j,k) >= fist)&&(fi(j,k) <= fiex)
                        g = 1;                        % tooth is in the cut
                    else
                        g = 0;                        % tooth is out of cut
                    end
                    %切厚计算
                    if k==NA+1
                        kxx(x,y,i,j,k) =  g*res_w*(Kt*cos(fi(j,k))+Kn*sin(fi(j,k)))*sin(fi(j,k));
                        kxy(x,y,i,j,k) =  g*res_w*(Kt*cos(fi(j,k))+Kn*sin(fi(j,k)))*cos(fi(j,k));
                        kyx(x,y,i,j,k) =  g*res_w*(-Kt*sin(fi(j,k))+Kn*cos(fi(j,k)))*sin(fi(j,k));
                        kyy(x,y,i,j,k) =  g*res_w*(-Kt*sin(fi(j,k))+Kn*cos(fi(j,k)))*cos(fi(j,k)); 
                    else
                        kxx(x,y,i,j,k) =  g*d_w*(Kt*cos(fi(j,k))+Kn*sin(fi(j,k)))*sin(fi(j,k));
                        kxy(x,y,i,j,k) =  g*d_w*(Kt*cos(fi(j,k))+Kn*sin(fi(j,k)))*cos(fi(j,k));
                        kyx(x,y,i,j,k) =  g*d_w*(-Kt*sin(fi(j,k))+Kn*cos(fi(j,k)))*sin(fi(j,k));
                        kyy(x,y,i,j,k) =  g*d_w*(-Kt*sin(fi(j,k))+Kn*cos(fi(j,k)))*cos(fi(j,k));
                    end

                end
            end

        end%时间循环结束
    end%切深循环结束
    fprintf('(1/2)瞬时切厚计算中：%d\n',stx+1-x)
end%转速循环结束
             

%-------- Begin of the proposed method --------
A0=[-invM*C/2, invM;
    C*invM*C/4-K,-C*invM/2];
I=eye (size (A0)) ;
invA0 = inv(A0);                               % invA0 = (A0)^(-1)

%初始化kjp
Kjp=zeros(2,2);
Kjp1=zeros(2,2);
%初始化B0j,B1j
B0j=zeros(2*num_modal,2*num_modal);
B1j=zeros(2*num_modal,2*num_modal);

%初始化ss,dc,ei
ss=zeros(stx+1,sty+1);                         
dc=zeros(stx+1,sty+1);
ei=zeros(stx+1,sty+1);

% start of computation
for x = 1 : stx+1                               % sweeping spindle speeds
    o = o_st +(x-1)*(o_fi-o_st)/stx;         % spindle speed
    T=60/o;                                  % 周期T
    dt=T/m;                                 % time step

    %------- Calculation of Φ0,Φ1,Φ2,Φ3 -------
    Fi0 = expm(A0*dt);                          % Φ0
    Fi1=  A0 \ (Fi0 - I);                     % Φ1
    Fi2 = A0 \ (Fi0 * dt - Fi1);             % Φ2
    Fi3 = A0 \ (Fi0 * dt^2 - 2 * Fi2) ;      % Φ3
    Fi4 = A0 \ (Fi0 * dt^3 - 3 * Fi3) ;      % Φ4
    Fi5 = A0 \ (Fi0 * dt^4 - 4 * Fi4) ;      % Φ5
    %----- The end of calculation of Φ0,Φ1,Φ2,Φ3 -----

    %提前计算各表达式系数
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


    for y = 1 : sty+1                               % sweeping depth of cuts
        w = w_st +(y-1)*(w_fi-w_st)/sty;         % 此时的切削深度
        NA=floor(w/d_w);                         %以d_w离散刀具轴向高度

        max_mjl_zhihou=max(max(max(mjl_zhihou(x,y,:,:,:))));%计算max(m_jl);
        D = zeros(num_modal*(max_mjl_zhihou+2),num_modal*(max_mjl_zhihou+2)); % initialize matrix D, the discrete map matrix
        D(2*num_modal+1:3*num_modal,1:num_modal)=eye(num_modal);
        D(3*num_modal+1:end,(2*num_modal+1):(end-num_modal))=eye(num_modal*(max_mjl_zhihou-1));
        
        Fi= eye(num_modal*(max_mjl_zhihou+2),num_modal*(max_mjl_zhihou+2));                          % construct transition matrix Fi

        for i = 1 : m
            %初始化系数不变项
            A0k=zeros(2*num_modal,2*num_modal);
            A1k=zeros(2*num_modal,2*num_modal);

            %计算系数不变项
            Kjp(1,1)=sum(sum(kxx(x,y,i,:,1:NA+1)));
            Kjp(1,2)=sum(sum(kxy(x,y,i,:,1:NA+1)));
            Kjp(2,1)=sum(sum(kyx(x,y,i,:,1:NA+1)));
            Kjp(2,2)=sum(sum(kyy(x,y,i,:,1:NA+1)));
            Kjp_mul0=G'*Kjp*G;

            Kjp1(1,1)=sum(sum(kxx(x,y,i+1,:,1:NA+1)));
            Kjp1(1,2)=sum(sum(kxy(x,y,i+1,:,1:NA+1)));
            Kjp1(2,1)=sum(sum(kyx(x,y,i+1,:,1:NA+1)));
            Kjp1(2,2)=sum(sum(kyy(x,y,i+1,:,1:NA+1)));
            Kjp_mul1=G'*Kjp1*G;

            A0k(num_modal+1:end,1:num_modal)=-Kjp_mul0;
            A1k(num_modal+1:end,1:num_modal)=-Kjp_mul1;


            Q_minus_2=Q_minus_2_1k*A1k+ Q_minus_2_0k*A0k;
            Q_minus_1=Q_minus_1_1k*A1k+ Q_minus_1_0k*A0k;
            Q_0=Q_0_1k*A1k+ Q_0_0k*A0k;
            Q_1=Q_1_1k*A1k+Q_1_0k*A0k;    

            
            Di_stable=zeros(2*num_modal,num_modal*(max_mjl_zhihou+2));
            Di_stable(:,1:4*num_modal) = ...
            [
                -(Q_1-I)\(Q_0+Fi0), ...
                -(Q_1-I)\Q_minus_1(:, 1:num_modal), ...
                -(Q_1-I)\Q_minus_2(:, 1:num_modal)
            ];


            %初始化系数变化项
            Djl=zeros(2*num_modal,num_modal*(max_mjl_zhihou+2));
            Djl_temp=zeros(2*num_modal,num_modal*(max_mjl_zhihou+2));


            %计算系数变化项  
            for j=1:N
                %为缩减计算规模，将同一刀齿上，按照时滞量划分区间，合并计算
                mjl_zhihou_find=mjl_zhihou(x,y,i,j,1:NA+1);
                [mjl_zhihou_value,mjl_zhihou_index,~]=unique(mjl_zhihou_find,'stable');
                mjl_zhihou_index_length=length(mjl_zhihou_index);

                %相同时滞量的合并计算
                for equal_mjl_index=1:mjl_zhihou_index_length
                    if equal_mjl_index==mjl_zhihou_index_length

                        Kjp(1,1)=sum(kxx(x,y,i,j,mjl_zhihou_index(equal_mjl_index):NA+1));
                        Kjp(1,2)=sum(kxy(x,y,i,j,mjl_zhihou_index(equal_mjl_index):NA+1));
                        Kjp(2,1)=sum(kyx(x,y,i,j,mjl_zhihou_index(equal_mjl_index):NA+1));
                        Kjp(2,2)=sum(kyy(x,y,i,j,mjl_zhihou_index(equal_mjl_index):NA+1));
            
                        Kjp1(1,1)=sum(kxx(x,y,i+1,j,mjl_zhihou_index(equal_mjl_index):NA+1));
                        Kjp1(1,2)=sum(kxy(x,y,i+1,j,mjl_zhihou_index(equal_mjl_index):NA+1));
                        Kjp1(2,1)=sum(kyx(x,y,i+1,j,mjl_zhihou_index(equal_mjl_index):NA+1));
                        Kjp1(2,2)=sum(kyy(x,y,i+1,j,mjl_zhihou_index(equal_mjl_index):NA+1));

                    else
                        Kjp(1,1)=sum(kxx(x,y,i,j,mjl_zhihou_index(equal_mjl_index):mjl_zhihou_index(equal_mjl_index+1)-1));
                        Kjp(1,2)=sum(kxy(x,y,i,j,mjl_zhihou_index(equal_mjl_index):mjl_zhihou_index(equal_mjl_index+1)-1));
                        Kjp(2,1)=sum(kyx(x,y,i,j,mjl_zhihou_index(equal_mjl_index):mjl_zhihou_index(equal_mjl_index+1)-1));
                        Kjp(2,2)=sum(kyy(x,y,i,j,mjl_zhihou_index(equal_mjl_index):mjl_zhihou_index(equal_mjl_index+1)-1));
            
                        Kjp1(1,1)=sum(kxx(x,y,i+1,j,mjl_zhihou_index(equal_mjl_index):mjl_zhihou_index(equal_mjl_index+1)-1));
                        Kjp1(1,2)=sum(kxy(x,y,i+1,j,mjl_zhihou_index(equal_mjl_index):mjl_zhihou_index(equal_mjl_index+1)-1));
                        Kjp1(2,1)=sum(kyx(x,y,i+1,j,mjl_zhihou_index(equal_mjl_index):mjl_zhihou_index(equal_mjl_index+1)-1));
                        Kjp1(2,2)=sum(kyy(x,y,i+1,j,mjl_zhihou_index(equal_mjl_index):mjl_zhihou_index(equal_mjl_index+1)-1));

                    end
                    Kjp_mul=G'*Kjp*G;
                    Kjp1_mul=G'*Kjp1*G;

                    B0j(num_modal+1:end,1:num_modal)=Kjp_mul;
                    B1j(num_modal+1:end,1:num_modal)=Kjp1_mul;
    
                    H_0=H_0_1j*(-B1j)+ H_0_0j*(-B0j);     
                    H_minus_1=H_minus_1_1j*(-B1j)+ H_minus_1_0j*(-B0j);    
                    H_minus_2=H_minus_2_1j*(-B1j)+H_minus_2_0j*(-B0j);
                    H_minus_3=H_minus_3_1j*(-B1j)+H_minus_3_0j*(-B0j);


                    Djl_temp(:,(num_modal*(mjl_zhihou_value(equal_mjl_index)-2)+1):(num_modal*(mjl_zhihou_value(equal_mjl_index)+2))) = ...
                    [
                     (Q_1-I)\H_minus_3(:, 1:num_modal), ...
                     (Q_1-I)\H_minus_2(:, 1:num_modal), ...
                     (Q_1-I)\H_minus_1(:, 1:num_modal), ...
                     (Q_1-I)\H_0(:, 1:num_modal)
                    ];

                    %Djl累加，只累加首行以缩减时间。
                    Djl(:,   (num_modal*(mjl_zhihou_value(equal_mjl_index)-2)+1)  :  (num_modal*(mjl_zhihou_value(equal_mjl_index)+2)))= ...
                        Djl(:,   (num_modal*(mjl_zhihou_value(equal_mjl_index)-2)+1)  :  (num_modal*(mjl_zhihou_value(equal_mjl_index)+2)))+ ...
                        Djl_temp(:, (num_modal*(mjl_zhihou_value(equal_mjl_index)-2)+1)  :  (num_modal*(mjl_zhihou_value(equal_mjl_index)+2)));
                end


            end
            Di_add=D;
            Di_add(1:2*num_modal,:)=Di_stable+Djl;%总的转移矩阵，系数不变项+系数变化项
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


