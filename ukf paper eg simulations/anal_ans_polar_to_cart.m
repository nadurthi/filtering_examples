mu_anal=[amux;amuy];
P_anal=[asigx,asigxy;asigxy,asigy];

mu_ut=[mu_ut(1);mu_ut(2)];
P_ut=P_ut;

mu_ckf=[mu_ckf(1);mu_ckf(2)];
P_ckf=P_ckf;

mu_nm4=[mu_nm4(1);mu_nm4(2)];
P_nm4=P_nm4;

mu_nm6=[mu_nm6(1);mu_nm6(2)];
P_nm6=P_nm6;

mu_gh3=[mu_gh3(1);mu_gh3(2)];
P_gh3=P_gh3;

mu_gh4=[mu_gh4(1);mu_gh4(2)];
P_gh4=P_gh4;

MU=[mu_anal,mu_ut,mu_ckf,mu_nm4,mu_nm6,mu_gh3,mu_gh4]

PU=sqrt([diag(P_anal),diag(P_ut),diag(P_ckf),diag(P_nm4),diag(P_nm6),diag(P_gh3),diag(P_gh4)])

100*abs(MU-repmat(MU(:,1),1,7))./repmat(MU(:,1),1,7)
100*abs(PU-repmat(PU(:,1),1,7))./repmat(PU(:,1),1,7)