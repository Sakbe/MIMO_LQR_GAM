%%
close all
clear all

%%% Desgin Kalman filter
load('shot_47353.mat');
load('ISTTOK_model_Send_pos.mat');
load('ISTTOK_model_Send_neg.mat');
%%
index1=1134;
index2=1394;
%%
I_vert=data.SendToVertical(index1:index2);
I_hor=data.SendToHorizontal(index1:index2);
%  I_vert=data.vert(index1:index2);
%  I_hor=data.hor(index1:index2);
R=data.R0(index1:index2);
Z=data.z0(index1:index2);
R_kal=data.R0_Kalman(index1:index2);
Z_kal=data.z0_Kalman(index1:index2);
time=1e-6*data.time(index1:index2);



%%
Ts=100e-6;
Input=[I_vert,I_hor,R,Z];
Input1=[I_vert,I_hor];
Outputs1=[R,Z];
exp=iddata(outputs,Input1,Ts);
[dummy,dummy1,x0] =compare(ss_pos,exp);

SYS_pos = ss(ss_pos.A,ss_pos.B , ss_pos.C, ss_pos.D, Ts);
B_errorPos = [ss_pos.B,ss_pos.B];
SYS_pos = ss(ss_pos.A,B_errorPos, ss_pos.C, zeros(2,4), Ts);
% [Q,Rr,N]=GetKalmanMatrix(Outputs1,Input1,Ts,time,ss_pos.D)
load('KalmanMAtrixes')
%[kest_pos, L_pos] = kalman(SYS_pos, w_var*eye(2), v_var*eye(2));
[kest_pos, L_pos] = kalman(SYS_pos, Q, Rr,N);
sys_posKAL=ss(kest_pos.A,kest_pos.B,kest_pos.C,kest_pos.D,Ts);	
%Y_pos=lsim(sys_posKAL,Input,time);
[Y_pos,lsimtime,X_pos_k]=lsim(sys_posKAL,Input,time,ss_pos.X0);


%%
SYS_neg = ss(ss_neg.A,ss_neg.B , ss_neg.C, ss_neg.D, Ts);
B_errorNeg = [ss_neg.B,ss_neg.B];
SYS_neg = ss(ss_neg.A,B_errorNeg, ss_neg.C, zeros(2,4), Ts);


%[Qn,Rn,Nn]=GetKalmanMatrix(Outputs1,Input1,Ts,time,ss_neg.D)
load('KalmanMAtrixes_neg');
[kest_neg, L_neg] = kalman(SYS_neg, Qn, Rn,Nn);

sys_negKAL=ss(kest_neg.A,kest_neg.B,kest_neg.C,kest_neg.D,Ts);	
[Y_neg,lsimtime,X_neg_k]=lsim(sys_negKAL,Input,time,ss_neg.X0);
return
%%
figure(1)

subplot(2,1,1)
title('Radial Centroid Position')
plot(Y_pos(:,1),'LineWidth',2)
hold on
plot(R)
plot(R_kal,'.-')
grid on
legend('Kalman','Real','Kalman_GAM')
subplot(2,1,2)
plot(Y_pos(:,2),'LineWidth',2)
hold on
plot(Z)
plot(Z_kal,'.-')
legend('Kalman','Real','KAlman_GAM')
 grid on

%%
close all
%%%%% Desenho do controlador%%%%%

time=0:100e-6:0.03;
R=[1e0,0.4e1];
%R=[1e-5 5e-5];
R=diag(R);
Q=[5e6,1e1,1e-1,1e1,1e1,1e1,1e1,1e1,1e1,1e1];
%Q=1e3*[1e4,1,1,1,1,1,1,1,1,1];
Q=diag(Q);
%Q = ss_pos.C'*ss_pos.C;

[K_pos,S,e] = dlqr(ss_pos.A,ss_pos.B,Q,R) ;
Ac = [(ss_pos.A-ss_pos.B*K_pos)];

N_pos=rscale(ss_pos.A,ss_pos.B,ss_pos.C,ss_pos.D,K_pos);
sys_cl_pos = ss(Ac,ss_pos.B,ss_pos.C,ss_pos.D,Ts);
ref=[0.025,0.0]'*ones(size(time));
[y_cntrl,timel,x_pos]=lsim(sys_cl_pos,N_pos*ref,time,ss_pos.X0);
inputs_pos=N_pos*ref-K_pos*x_pos';

figure(1)
subplot(2,2,1)
plot(time,y_cntrl(:,1))
hold on
plot(time,ref(1,:))
grid on
ylabel('R [m]')
subplot(2,2,2)
plot(time,y_cntrl(:,2))
hold on
plot(time,ref(2,:))
grid on
ylabel('Z [m]')
subplot(2,2,3)
plot(time,inputs_pos(1,:))
grid on
ylabel('Vertical curr [A]')
subplot(2,2,4)
plot(time,inputs_pos(2,:))
grid on
ylabel('Horizontal curr [A]')

%% negativo
close all
time=0:100e-6:0.03;
R=[1e-6,1e-7];
R=diag(R);
Q=[2e4,1e2,1e2,1,1,1,1,1,1,1];
Q=1e1*diag(Q);
%Q = ss_pos.C'*ss_pos.C;

[K_neg,S,e] = dlqr(ss_neg.A,ss_neg.B,Q,R) ;
Ac = [(ss_neg.A-ss_neg.B*K_neg)];

N_neg=rscale(ss_neg.A,ss_neg.B,ss_neg.C,ss_neg.D,K_neg);
sys_cl_neg = ss(Ac,ss_neg.B,ss_neg.C,ss_neg.D,Ts);
ref=([0.025,0.001])'*ones(size(time));
[y_cntrl,timel,x_neg]=lsim(sys_cl_neg,N_neg*ref,time,ss_neg.X0);
inputs_neg=N_neg*ref-K_neg*x_neg';


figure(1)
subplot(2,2,1)
plot(time,y_cntrl(:,1))
hold on
plot(time,ref(1,:))
grid on
ylabel('R [m]')
subplot(2,2,2)
plot(time,y_cntrl(:,2))
hold on
plot(time,ref(2,:))
grid on
ylabel('Z [m]')
subplot(2,2,3)
plot(time,inputs_neg(1,:))
grid on
ylabel('Vertical curr [A]')
subplot(2,2,4)
plot(time,inputs_neg(2,:))
grid on
ylabel('Horizontal curr [A]')

%% Check the states

%%%%%positive
[yl,timel,x]=lsim(ss_pos,Input1,time,ss_pos.X0);

close all
figure(1)

subplot(2,1,1)
plot(X_pos_k(:,5),'LineWidth',2)
hold on
plot(R_kal,'.-')
plot(Y_pos(:,7))
%plot(x(:,3));
grid on
legend('State1','State1_GAM','Output7')
subplot(2,1,2)
plot(X_pos_k(:,6),'LineWidth',2)
hold on
plot(Z_kal,'.-')
plot(Y_pos(:,8))
%plot(x(:,4));
legend('State2','State2_GAM','Output8')
 grid on
 
%%%%% negative