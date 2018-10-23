pavg=[65.4286759	72.02136632	79.87460253	79.00713629	79.50538682	79.51241438	79.51459892	79.48129013;
   90.24999243	97.06782081	108.4107063	106.4579647	105.6739032	105.5404136	105.6715454	105.6293123;
   103.0853871	110.5017313	765.7183111	964.3927261	211.9030357	121.0629679	119.6584865	119.7756651;
   115.4791723	123.327338	989.154149	685.9077093	245.3032388	138.8293352	135.8998032	142.8502538];

vavg=[20.78258358	31.43948604	25.59421753	25.17234472	25.24091346	25.23063001	25.23121224	25.21838292;
      23.03521815	35.12102663	211.1159849	29.49112304	28.8729117	28.69644174	28.69031089	28.68892957;
      24.24005636	36.31060805	5843.93763	7873.433814	2773.34731	32.64917142	30.9218027	34.19981549;
      24.16230844	36.19659285	17330.23325	12849.90453	6127.861146	2153.53072	34.73475825	2870.563981];
  
  Omgavg=[0.038512678	0.049787933	0.048943926	0.048517059	0.047431176	0.047406904	0.047406351	0.047404026;
          0.03949075	0.054230534	0.07872852	0.061967224	0.054821878	0.054878811	0.05447234	0.054529476;
          0.039661026	0.057412341	0.775179026	1.171955879	0.511506479	0.065734814	0.054596582	0.085456759;
          0.03909923	0.06120062	2.873242402	2.396186044	1.329138111	0.636459943	0.090532177	0.806625346];
      
pmax=[1281.164786	1281.164786	1281.164786	1281.164786	1281.164786	1281.164786	1281.164786	1281.164786;
      1391.537223	1391.537223	1391.537223	1391.537223	1391.537223	1391.537223	1391.537223	1391.537223;
       1308.673482	1308.673482	4343.332951	4695.087518	1308.673482	1308.673482	1308.673482	1308.673482;
      1525.392581	1525.392581	6052.154647	3470.636425	1525.392581	1525.392581	1525.392581	1525.392581];    

vmax=[66.22953428	106.096552	78.43725227	77.0518492	76.97765209	76.995928	76.99862558	76.97309537;
      70.44512802	107.0865302	2113.412079	83.81964326	81.48516372	80.72987417	80.51431169	80.50901679;
      74.76798755	108.7172451	26109.89729	38396.30175	17135.1672	91.56449803	84.18565946	100.4312904;
      81.71726331	120.557747	96847.82403	60870.55373	30490.58039	11530.74508	101.1470216	13703.46706]; 
  
Omgmax=[0.115692144	0.185323695	0.14760214	0.14555594	0.141445687	0.141694846	0.141700726	0.141668492;
        0.120665954	0.167055982	0.413845613	0.186413196	0.161165687	0.159334821	0.157560552	0.157265456;
        0.12647605	0.171902758	3.082096163	4.748360134	2.339197693	0.271598102	0.154558375	0.277406356;
        0.123294841	0.183960109	11.75269197	9.821763034	5.959563519	3.245037176	0.320932969	3.39044536];
    
    
%     pavg(:,[2,8])=[];
figure
    dt=[1,2.5,3.9,5];
    plot(dt,pavg(:,1),'o--',dt,pavg(:,3),'s--',dt,pavg(:,4),'*--',dt,pavg(:,5),'+--',dt,pavg(:,6),'d--',dt,pavg(:,7),'^--','MarkerSize',16,'Linewidth',2.5)
    xlabel('time interval T(s)')
    ylabel('||RMSE_{pos}||_2')
    legend('pf','ckf','ukf','cut4','cut6','cut8','Location','NorthWest')
    plot_prop_paper
%   saveas(gca,'AirTraffic_pos_avg','eps')   
        figure
    plot(dt,vavg(:,1),'o--',dt,vavg(:,3),'s--',dt,vavg(:,4),'*--',dt,vavg(:,5),'+--',dt,vavg(:,6),'d--',dt,vavg(:,7),'^--','MarkerSize',16,'Linewidth',2.5)
    xlabel('time interval T(s)')
    ylabel('||RMSE_{vel}||_2')
    legend('pf','ckf','ukf','cut4','cut6','cut8','Location','NorthWest')
    plot_prop_paper
%      saveas(gca,'AirTraffic_vel_avg','eps')
            figure
    plot(dt,Omgavg(:,1),'o--',dt,Omgavg(:,3),'s--',dt,Omgavg(:,4),'*--',dt,Omgavg(:,5),'+--',dt,Omgavg(:,6),'d--',dt,Omgavg(:,7),'^--','MarkerSize',16,'Linewidth',2.5)
    xlabel('time interval T(s)')
    ylabel('||RMSE_{Omg}||_2','Interpreter','tex')
    legend('pf','ckf','ukf','cut4','cut6','cut8','Location','NorthWest')
    plot_prop_paper
%      saveas(gca,'AirTraffic_omg_avg','eps')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    dt=[1,2.5,3.9,5];
    plot(dt,pmax(:,1),'o--',dt,pmax(:,3),'s--',dt,pmax(:,4),'*--',dt,pmax(:,5),'+--',dt,pmax(:,6),'d--',dt,pmax(:,7),'^--','MarkerSize',16,'Linewidth',2.5)
    xlabel('time interval T(s)')
    ylabel('||RMSE_{pos}||_{\infty}')
    legend('pf','ckf','ukf','cut4','cut6','cut8','Location','NorthWest')
    plot_prop_paper
%      saveas(gca,'AirTraffic_pos_max','eps')
        figure
    plot(dt,vmax(:,1),'o--',dt,vmax(:,3),'s--',dt,vmax(:,4),'*--',dt,vmax(:,5),'+--',dt,vmax(:,6),'d--',dt,vmax(:,7),'^--','MarkerSize',16,'Linewidth',2.5)
    xlabel('time interval T(s)')
    ylabel('||RMSE_{vel}||_{\infty}')
    legend('pf','ckf','ukf','cut4','cut6','cut8','Location','NorthWest')
    plot_prop_paper
%      saveas(gca,'AirTraffic_vel_max','eps')
            figure
    plot(dt,Omgmax(:,1),'o--',dt,Omgmax(:,3),'s--',dt,Omgmax(:,4),'*--',dt,Omgmax(:,5),'+--',dt,Omgmax(:,6),'d--',dt,Omgmax(:,7),'^--','MarkerSize',16,'Linewidth',2.5)
    xlabel('time interval T(s)')
    ylabel('||RMSE_{Omg}||_{\infty}','Interpreter','tex')
    legend('pf','ckf','ukf','cut4','cut6','cut8','Location','NorthWest')
    plot_prop_paper
%      saveas(gca,'AirTraffic_omg_max','eps')
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%5
    %%
    figure
    dt=[1,2.5,3.9,5];
    plot(dt,pavg(:,1),'o--',dt,pavg(:,5),'+--',dt,pavg(:,6),'d--',dt,pavg(:,7),'^--','MarkerSize',16,'Linewidth',2.5)
    xlabel('time interval T(s)')
    ylabel('||RMSE_{pos}||_2')
    legend('pf','cut4','cut6','cut8','Location','NorthWest')
    plot_prop_paper
%   saveas(gca,'AirTraffic_pos_avg_pfcut','eps')   
        figure
    plot(dt,vavg(:,1),'o--',dt,vavg(:,5),'+--',dt,vavg(:,6),'d--',dt,vavg(:,7),'^--','MarkerSize',16,'Linewidth',2.5)
    xlabel('time interval T(s)')
    ylabel('||RMSE_{vel}||_2')
    legend('pf','cut4','cut6','cut8','Location','NorthWest')
    plot_prop_paper
%      saveas(gca,'AirTraffic_vel_avg_pfcut','eps')
            figure
    plot(dt,Omgavg(:,1),'o--',dt,Omgavg(:,5),'+--',dt,Omgavg(:,6),'d--',dt,Omgavg(:,7),'^--','MarkerSize',16,'Linewidth',2.5)
    xlabel('time interval T(s)')
    ylabel('||RMSE_{Omg}||_2','Interpreter','tex')
    legend('pf','cut4','cut6','cut8','Location','NorthWest')
    plot_prop_paper
%      saveas(gca,'AirTraffic_omg_avg_pfcut','eps')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    dt=[1,2.5,3.9,5];
    plot(dt,pmax(:,1),'o--',dt,pmax(:,5),'+--',dt,pmax(:,6),'d--',dt,pmax(:,7),'^--','MarkerSize',16,'Linewidth',2.5)
    xlabel('time interval T(s)')
    ylabel('||RMSE_{pos}||_{\infty}')
    legend('pf','cut4','cut6','cut8','Location','NorthWest')
    plot_prop_paper
%     saveas(gca,'AirTraffic_pos_max_pfcut','eps')
        figure
    plot(dt,vmax(:,1),'o--',dt,vmax(:,5),'+--',dt,vmax(:,6),'d--',dt,vmax(:,7),'^--','MarkerSize',16,'Linewidth',2.5)
    xlabel('time interval T(s)')
    ylabel('||RMSE_{vel}||_{\infty}')
    legend('pf','cut4','cut6','cut8','Location','NorthWest')
    plot_prop_paper
%      saveas(gca,'AirTraffic_vel_max_pfcut','eps')
            figure
    plot(dt,Omgmax(:,1),'o--',dt,Omgmax(:,5),'+--',dt,Omgmax(:,6),'d--',dt,Omgmax(:,7),'^--','MarkerSize',16,'Linewidth',2.5)
    xlabel('time interval T(s)')
    ylabel('||RMSE_{Omg}||_{\infty}','Interpreter','tex')
    legend('pf','cut4','cut6','cut8','Location','NorthWest')
    plot_prop_paper
%      saveas(gca,'AirTraffic_omg_max_pfcut','eps')