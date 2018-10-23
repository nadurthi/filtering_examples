
Q=0.2*diag([2,2,2,2,2]);
A(1,:)=[1.0688    1.3049    1.1741    1.1739    1.1779    1.1783    1.1784  1.1780]
B(1,:)=[4.9272    5.6395    5.0028    5.0072    4.9460    4.9453    4.9454  4.9450]

Q=0.3*diag([2,2,2,2,2]);
A(2,:)=[ 1.2044    1.5358    1.4718    1.4722    1.4609    1.4652    1.4630 1.4606]
B(2,:)=[ 4.3898    4.3898    4.3898    4.3898    4.3898    4.3898    4.3898 4.3898]

Q=0.4*diag([2,2,2,2,2]);
A(3,:)=[1.3139    1.7724    2.1417    2.2496    1.8011    1.8446    1.7755  1.7682]
B(3,:)=[3.7480    6.0897   11.6603   12.9487    6.9839    8.5146    5.3499  5.3721]

Q=0.5*diag([2,2,2,2,2]);
A(4,:)=[1.7829    2.3009    4.4403    4.2587    3.2334    3.2895    2.8369   2.8580]
B(4,:)=[6.4045    7.3494   24.6523   22.0223   17.8672   19.1010   13.1671   12.9917]

Q=0.6*diag([2,2,2,2,2]);
A(5,:)=[1.6573    2.3534    5.8923    6.7781    3.8212    3.3771    2.9894  3.0719]
B(5,:)=[4.7714    5.7853   28.7206   32.9847   16.7081   11.7563    9.9857   10.7187]


Q=0.7*diag([2,2,2,2,2]);
A(6,:)=[1.8163    2.5227   11.2616   11.5541    3.7580    3.4324    3.4898    3.3718]
B(6,:)=[4.7020    6.7107   39.3925   41.5242   12.0367   10.8983   11.3163    10.5505] 

A(:,[2,8])=[];
B(:,[2,8])=[];

dt=[0.2,0.3,0.4,0.5,0.6,0.7];
figure(1)
plot(dt,A(:,1),'o--',dt,A(:,2),'s--',dt,A(:,3),'*--',dt,A(:,4),'^--',dt,A(:,5),'+--',dt,A(:,6),'d--','linewidth',2.5,'MarkerSize',15)
legend('pf','ckf','ukf','cut4','cut6','cut8')
xlabel('Measurement time interval')
ylabel('||RMSE||_2')
plot_prop_paper

figure(2)
plot(dt,B(:,1),'o--',dt,B(:,2),'s--',dt,B(:,3),'*--',dt,B(:,4),'^--',dt,B(:,5),'+--',dt,B(:,6),'d--','linewidth',2.5,'MarkerSize',15)
legend('pf','ckf','ukf','cut4','cut6','cut8')
xlabel('Measurement time interval')
ylabel('||RMSE||_\infty')
plot_prop_paper

figure(3)
plot3(avg_traj_mupf(:,1),avg_traj_mupf(:,2),avg_traj_mupf(:,3),'-.',avg_traj_ckf(:,1),avg_traj_ckf(:,2),avg_traj_ckf(:,3),avg_traj_ukf(:,1),avg_traj_ukf(:,2),avg_traj_ukf(:,3),avg_traj_cut4(:,1),avg_traj_cut4(:,2),avg_traj_cut4(:,3),avg_traj_cut6(:,1),avg_traj_cut6(:,2),avg_traj_cut6(:,3),avg_traj_cut8(:,1),avg_traj_cut8(:,2),avg_traj_cut8(:,3))
legend('PF-mu','ckf','ukf','cut4','cut6','cut8')

 t=time.tspan;
plot(t,est_fin_mupf(:,p),'-.',t,est_fin_ckf(:,p),t,est_fin_ukf(:,p),t,est_fin_cut4(:,p),t,est_fin_cut6(:,p),t,est_fin_cut8(:,p))
legend('PF-mu','ckf','ukf','cut4','cut6','cut8')


[t,x]=ode45(model.fx,[0:0.01:30],[1.50887,-1.531271,25.46091,10,28]',model.Qt_sq);
plot3(x(:,1),x(:,2),x(:,3))