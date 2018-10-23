% A class that containts all targets and all sensors
%         SensMsgs={'OutOfFOV','PseudoUpdate','SimpleTask'};
%         SensStates={'HasToTakeMeas','AlreadyTakenMeas'}
%         SensTypes={'Move','Stationary'};
%         TargTypes={'Move','Stationary','Virtual'};
%         TargStates={'CurrAt_tk','aPrioriAt_tk','aPostAt_tk'};
%         UQType={'UT','CUT4','CUT6','CUT8'};
%         InfoType={'FIM','MI'};
%         AvailSensors={'Range+Bearing','Range','Bearing'};
%         AvailSensors_dim={2,1,1};
%         AvailTargDyn={'UM','CT','NoDyn'};
%         AvailTargDyn_dim={4,5,2};
%         FOVpenaltytypes={'Simple','QuadraticOUT','NoPenalty'};
%         SensConstraintTypes={'1Sensor->1Target','1Sensor->AllTarget','1Sensor->1Target/AllTarget'};
%
%        If a sensor has more than 1 mode , add it as another sensor. Then
%        in the tasking phase put a constraint that only one of these
%        sensors can operate at a given time. The sensor variables are
%        arranged as
%              T        T       T      C        C    ||
%         --   s1      s2      s3      s4      s5    ||
%         t1   s11     s21     s31     s41     s51   ||  1    5   9    13 .
%         t2   s12     s22     s32     s42     s52   ||  2    6   10   14 .
%         t3   s13     s23     s33     s43     s53   ||  3    7   11   15 .
%         t4   s14     s24     s34     s44     s54   ||  4    8   12   16 .

% T means a traking sensors that can see only one object at once
% C means a coverage sensor ... i.e. can see all within the FOV
% Constraints: s11+s12+s13+s14<=1  for T sensors
% S41,s42,s43,s44   for eg: s43=0 if sensor 4 cannot see the 3rd target
% The variables are stacked as [s11,s12,s13,s14  ,s21,s22,s23,s24
% ,s31,s32,s33,s34  ,....] Nt x Ns

classdef TargSens
    properties
        BackUpSimProps
        BackUpTargProps
        BackUpSensProps
        
        SimProps
        TargProps
        SensProps
        
        BackUp_targMeans
        BackUp_targCovss
        
        A
        b
        C
        d
        FIM;
        MUvec
        MU
        
        
        %         xlim
        %         ylim
        %         Ns
        %         Nt
        %         Ntvec_deleted  %operate only on this vector of valid targets
        %         x0
        %         P0
        %         xk  %current state estimates of the targets    xk={{k}{targid}}
        %         Pk  % Pk={{k}{targid}}
        %         xktruth
        %         NEWxktruth
        %         hn
        %         fn
        %         sk  % current position of the sensor
        %         alphak % half anfle of FOV
        %         rmaxk  % max range of FOV
        %         dirk  % pointing angle of sensor
        %         TargType   % (obj.Nt,[1,2])=(Move/Station,Dynamic Model) indexes
        %         TargState  %{time}(time,:)current state : is is propogated, is it meas updated, or nothing done,
        %         % this isvec of structs
        %         SensState   %current state : has to take measurement, already taken meas...
        %         SensType  % {Move/Station,Sensor Model}
        %         SensConstraints  % description of what kind of sensor it is
        %         SensModeConstraintsVec  %if sensor is multimodal  [id of sens1_mode1,id of sens1_mode2; ... ] so mode1 and mode2 cannot work together
        %         Nmodeconstraints % just the size of SensModeConstraintsVec
        %         MU  %sensor-target decision  matrix  targs along rows, sesnsor along columns
        %         MUvec % the vec version of MU decision matrix
        %         MUindex  % matrix that contains the indices of the MUvec
        %         TempSensConfig
        %         R
        %         Q
        %         G  %penalty for FOV out sens
        %         t0
        %         tk  % the current tiems step
        %         tf
        %         nt
        %         Nt_t0  %the starting time of the target. esp for targets added in the middle
        %         Tvec
        %         dt
        %         SensTargtask  %[sensors,targets]
        %         yk % current measurement
        %         quadpts
        %         A  %As<=b constraints
        %         b
        %         C %Cs==d
        %         d
        %         FIM % the latest information matrix
        %         SimMode  % different modes of simulations
    end
    properties(Constant)
        
        
        SensMsgs={'OutOfFOV','PseudoUpdate','SimpleTask'};
        
        SensConfigStates={'HasToTakeMeas','AlreadyTakenMeas'}
        SensTypes={'Move','Stationary'};
        TargTypes={'Move','Stationary','VirtualMove','VirtualMove'};
        TargConfigStates={'CurrAt_tk','aPrioriAt_tk','aPostAt_tk'};
        
        
        
        UQType={'UT','CUT4','CUT6','CUT8'};
        InfoType={'FIM','MI'};
        AvailSensors={'Range+Bearing','Range','Bearing'};
        AvailSensors_dim={2,1,1};
        AvailTargDyn={'UM','CT','NoDyn','diffdynCf1','diffdynCf0','diffdynCf2'};
        AvailTargDyn_dim={4,5,2};
        
        FOVpenaltytypes={'Simple','QuadraticOUT','NoPenalty'};
        
        SensorTaskingMethods={'MIUB','FIM'}
        
        SensConstraintTypes={'1Sensor->1Target','1Sensor->AllTarget','1Sensor->AllTarget/1Target'};
        % 1 sensor can see only 1 target within FOV
        % 1 sensor can see all target within FOV
        % 1 sensor can choose see only 1 target or all target  within FOV
        
        
        
    end
    methods
        
        %============================================================================
        function obj=TargSens(t0,dt,tf)
            
            obj.SimProps.SimMode=[];
            obj.SimProps.SensorMsgsAvail={'OutOfFOV','PseudoUpdate','SimpleTask'};
            obj.SimProps.SensorConfigStatesAvail={'HasToTakeMeas','AlreadyTakenMeas'};
            obj.SimProps.TargetConfigStatesAvail={'CurrAt_tk','aPrioriAt_tk','aPostAt_tk'};
            obj.SimProps.SensorDynTypesAvail={'Move','Stationary'};
            obj.SimProps.TargtDynTypesAvail={'Move','Stationary','VirtualMove','VirtualStationary'};
            obj.SimProps.SigmaPntsTypeAvail={'UT','CUT4','CUT6','CUT8'};
            obj.SimProps.SensorTasksAvail={'MIUB','FIM'};
            obj.SimProps.SensorModelTypesAvail={'Range+Bearing','Range','Bearing','X+Y'};
            obj.SimProps.TargetModelTypesAvail={'UM','CT','NoDyn'};
            obj.SimProps.TargetVisibilityStateAvail={'Lost','Tracked','Invisible','NotTracked'};
            obj.SimProps.SensorDynMotionModelTypes={'UAVgridOnly','UAVgrid+smoothTraj','None'};
            obj.SimProps.SensorModeConstraintTypesAvail={'1ModeOperation'};
            obj.SimProps.SensorFOVtypesAvail={'1Sensor->1Target','1Sensor->AllTarget','1Sensor->AllTarget/1Target'};
            
            obj.SimProps.MUindex=[];
            
            obj.SimProps.Xboundary=[];
            obj.SimProps.Yboundary=[];
            
            obj.SimProps.SigmaPntFunction=[];
            
            obj.SimProps.Time.t0=t0;
            obj.SimProps.Time.tf=tf;
            obj.SimProps.Time_dt=dt;
            obj.SimProps.Time.tvec=t0:dt:tf;
            obj.SimProps.Time.nsteps=length(obj.SimProps.Time.tvec);
            obj.SimProps.Time.tk=1; %set the curr time to 0 or Tvec(1)
            
            obj.SimProps.XY=[];
            obj.SimProps.XYindex=[];
            obj.SimProps.XYADJindex=[];
            obj.SimProps.delx=0;
            
            %Planning/Tasking time steps 
            obj.SimProps.NoPlanningSteps=5;
            
            obj.TargProps.NumbTargets=0;
            obj.TargProps.TargetID=[];
            obj.TargProps.MeanState_allk=cell(1,1); %one cell for each each target
            obj.TargProps.CovState_allk=cell(1,1); %one cell for each each target
            obj.TargProps.VisibilityState_allk=cell(1,1); %one cell for each each target   {'Lost','Tracked','Invisible','NotTracked'};
            obj.TargProps.TruePosState_allk=cell(1,1); %one cell for each each target
            obj.TargProps.TargetDynTypes_allk=cell(1,1); %one cell for each each target {'Move','Stationary','VirtualMove','VirtualStationary'}
            obj.TargProps.Q_allk=cell(1,1); %one cell for each each target
            obj.TargProps.TargetModel_allk=cell(1,1); %one cell for each each target
            obj.TargProps.TargetModelDim_allk=cell(1,1); %one cell for each each target
            obj.TargProps.CovLimit=cell(1,1); % The maximum eigen value of the covariance beyond which it is capped 
            obj.TargProps.TargetModelParas_allk=cell(1,1);
            
            
            obj.SensProps.NumbSensors=0;
            obj.SensProps.SensorsID=[];
            obj.SensProps.SensorDynTypes_allk=cell(1,1); %one cell for each each sensor at each time step
            obj.SensProps.R_allk=cell(1,1); %one cell for each each sensor at each time step
            obj.SensProps.y_allk=cell(1); %each cell has measurements for a time step. at each time step cell(targid,sensid)
            obj.SensProps.SensorNmodes_allk=cell(1,1);
            
            %obj.SensProps.SensorPlannedConfigTasks.SensorPos=cell(1,1); %each cell for 1 sensor.
            
            obj.SensProps.SensorModel_allk=cell(1,1); %one cell for each each target
            obj.SensProps.SensorModelDim_allk=cell(1,1);
            obj.SensProps.SensorFOVtype_allk=cell(1,1); %{'1Sensor->1Target','1Sensor->AllTarget','1Sensor->AllTarget/1Target'};
            obj.SensProps.SensorFOVpenalty_allk=cell(1,1);
            
            obj.SensProps.SensorModeConstraints_allk=cell(1,1);
            obj.SensProps.SensorOnOff_allk=cell(1,1);
            
            
            obj.SensProps.Amat=[];  %Ax<=b
            obj.SensProps.bvec=[];
            obj.SensProps.Cmat=[];  %Cx==d
            obj.SensProps.dvec=[];
            
            % how does the UAV move {'UAVgridOnly', position for all time steps,}
            % {'UAVgrid+smoothTraj', position for all time steps,PlottingTrajectory}
            % obj.SensProps('SensorDynTypes_allk') has to be Move inorder
            % for this to work
            obj.SensProps.SensorDynMotionModel_allk=cell(1,1);   %{'UAVgridOnly','UAVgridOnlyNoDir'}
            
            % the planned trajectories/values are loladed into these
            % vectors
            obj.SensProps.SensorPos_allk=cell(1,1);    %{sensid}(tk,:)
            obj.SensProps.SensorRmax_allk=cell(1,1);   %{sensid}(tk,:)
            obj.SensProps.SensorDir_allk=cell(1,1);    %{sensid}(tk,:)
            obj.SensProps.SensorAlpha_allk=cell(1,1);  %{sensid}(tk,:)
            
            
            obj.SensProps.SensorTaskTargets_allk=cell(1,1); % each cell is one sensor, for each sensor at each time step the list targets it is tasked to
            
            
            
            obj.BackUpSimProps = [];
            obj.BackUpTargProps = [];
            obj.BackUpSensProps = [];
            
            obj.BackUp_targMeans = [];
            obj.BackUp_targCovss = [];
            
            
            %             obj.Ns=0;
            %             obj.Nt=0;
            %             obj.t0=t0;
            %             obj.tf=tf;
            %             obj.dt=dt;
            %             obj.Tvec=t0:dt:tf;
            %             obj.nt=length(obj.Tvec);
            %             obj.tk=1;  %set the curr time to 0 or Tvec(1)
            %             obj.x0={};
            %             obj.P0={};
            %             obj.xk=cell(1,1); %xk{obj.Nt}=[xk0;xk1;xk2].....
            %             obj.Pk=cell(1,1);
            %             obj.xktruth=cell(1,1);
            %             obj.TargType=[0,0];  %(Move/Station,Dynamic Model)
            %             obj.SensType=[0,0];  %(Move/Station,Sensor Model)
            
            %             obj.TargState=cell(1,1);
            %             obj.SensState=cell(1,1);   %current state : has to take measurement, already taken meas...
            
            %             obj.R=cell(1,1); %one for each sensor
            %             obj.Q=cell(1,1);  % one for each target
            %             obj.yk=cell(1,1);  % each time the measurements are loaded into this
            %             obj.SensTargtask=cell(1,1); %each cell  [sensors,targets] is a vector
            %             obj.hn=0; % this is a vector of sens dims
            %             obj.fn=0;
            %             obj.Ntvec_deleted=[];
            %             obj.Nt_t0=0;
            %             obj.sk=cell(1,1);
            %             obj.alphak=cell(1,1);
            %             obj.rmaxk=cell(1,1);
            %             obj.dirk=cell(1,1);
            
            %             obj.SensConstraints=cell(1,1);
            %             obj.C=[];
            %             obj.A=[];
            %             obj.b=[];
            %             obj.d=[];
            %
            %             obj.SensModeConstraintsVec=cell(1,1);
            %             obj.Nmodeconstraints=0;
            %             obj.MUindex=0;
            %
            %             obj.TempSensConfig.Configtk=0; % backed up time step to be careful
            %             obj.TempSensConfig.Config=cell(1,4); %Temporarily store the current config of all the sensors  {alpha,Rmax,pos,dirn}
            %             obj.FIM=0;
            %             obj.SimMode=0;
        end
        %=========================================================================================================
        function obj=SetupGrid(obj,dx)
            
            Xlim=obj.SimProps.Xboundary;
            [X,Y]=meshgrid(Xlim(1):dx:Xlim(2));
            XY=[Y(:),X(:)];
            XYindex=[1:1:size(XY,1)]';
            % [XY,XYindex]
            ng=length(XYindex);
            XYADJindex=cell(ng,1);
            
            for i=1:1:length(XYindex)
                xy=XY(i,:);
                p=[];
                
                r=xy+[dx,0];
                if r(1)>Xlim(2) || r(2) >Xlim(2) || r(1)<Xlim(1) || r(2)<Xlim(1)
                    
                else
                    [~,ind]=min(sqrt(sum((repmat(r,ng,1)-XY).^2,2)));
                    p=horzcat(p,ind(1));
                end
                
                r=xy+[0,dx];
                if r(1)>Xlim(2) || r(2) >Xlim(2) || r(1)<Xlim(1) || r(2)<Xlim(1)
                    
                else
                    [~,ind]=min(sqrt(sum((repmat(r,ng,1)-XY).^2,2)));
                    p=horzcat(p,ind(1));
                end
                
                r=xy+[-dx,0];
                if r(1)>Xlim(2) || r(2) >Xlim(2) || r(1)<Xlim(1) || r(2)<Xlim(1)
                    
                else
                    [~,ind]=min(sqrt(sum((repmat(r,ng,1)-XY).^2,2)));
                    p=horzcat(p,ind(1));
                end
                
                r=xy+[0,-dx];
                if r(1)>Xlim(2) || r(2) >Xlim(2) || r(1)<Xlim(1) || r(2)<Xlim(1)
                    
                else
                    [~,ind]=min(sqrt(sum((repmat(r,ng,1)-XY).^2,2)));
                    p=horzcat(p,ind(1));
                end
                
                XYADJindex{i}=sort(p);
            end

            obj.SimProps.XY=XY;
            obj.SimProps.XYindex=XYindex;
            obj.SimProps.XYADJindex=XYADJindex;
            obj.SimProps.delx=dx;
        end
        %============================================================================
        function obj=Set_Sigmapts(obj,type)
            
            switch(type)
                case 'UT'   %UT
                    obj.SimProps.SigmaPntFunction=@(mu,P)UT_sigmapoints(mu(:),P,2);
                case 'CUT4'   %CUT4
                    disp('lol')
                case 'GH3'   %CUT4
                    obj.SimProps.SigmaPntFunction=@(mu,P)GH_points(mu(:),P,3);
                case 'GH4'   %CUT4
                    obj.SimProps.SigmaPntFunction=@(mu,P)GH_points(mu(:),P,4);    
                case 'GH5'   %CUT4
                    obj.SimProps.SigmaPntFunction=@(mu,P)GH_points(mu(:),P,5);    
            end
        end
        %============================================================================
        function obj=Add_Target(obj,targ)  % targtype=Move/Stationary , truth is for all times
            %this way the target can be initiated at any time
            
            %             if max(strcmpi(TargType,obj.TargTypes))==0
            %                 error('Wrong Targtype');
            %             end
            %             if max(strcmpi(DynType,obj.AvailTargDyn))==0
            %                 error('Wrong DynType');
            %             end
            
            obj.TargProps.NumbTargets=obj.TargProps.NumbTargets+1;
            
            n=obj.TargProps.NumbTargets;
            
            
            obj.TargProps.TargetID(n)=targ.TargetID;
            obj.TargProps.MeanState_allk{n}=targ.MeanState_allk;
            obj.TargProps.CovState_allk{n}=targ.CovState_allk;
            obj.TargProps.VisibilityState_allk{n}=targ.VisibilityState_allk;
            obj.TargProps.TruePosState_allk{n}=targ.TruePosState_allk;
            obj.TargProps.TargetDynTypes_allk{n}=targ.TargetDynTypes_allk;
            obj.TargProps.Q_allk{n}=targ.Q_allk;
            obj.TargProps.TargetModel_allk{n}=targ.TargetModel_allk;
            obj.TargProps.TargetModelDim_allk{n}=targ.TargetModelDim_allk;
            obj.TargProps.CovLimit{n}=targ.CovLimit;
            
            obj.TargProps.TargetModelParas_allk{n}=targ.TargetModelParas_allk;
            
        end
        %============================================================================
        function cont=Load2ContainerCellAt(cont,keyy,NewCellValue,n)
            
            S=cont(keyy);
            S{n}=NewCellValue;
            cont(keyy)=S;
            
        end
        %============================================================================
        function X=GetContCell_Nt(cont,keyy,n)
            S=cont(keyy);
            X=S{n};
        end
        %============================================================================
        function obj=Add_Sensor(obj,sens)
            
            
            obj.SensProps.NumbSensors=obj.SensProps.NumbSensors+1;
            n=obj.SensProps.NumbSensors;
            
            obj.SensProps.SensorsID(n)=sens.SensorsID;
            obj.SensProps.SensorDynTypes_allk{n}=sens.SensorDynTypes_allk;
            obj.SensProps.SensorNmodes_allk{n}=sens.SensorNmodes_allk;
            obj.SensProps.R_allk{n}=sens.R_allk;
            obj.SensProps.SensorModel_allk{n}=sens.SensorModel_allk;
            
            obj.SensProps.SensorFOVpenalty_allk{n}=sens.SensorFOVpenalty_allk;
            
            obj.SensProps.SensorFOVtype_allk{n}=sens.SensorFOVtype_allk;
            obj.SensProps.SensorDynMotionModel_allk{n}=sens.SensorDynMotionModel_allk;
            obj.SensProps.SensorPos_allk{n}=sens.SensorPos_allk;
            obj.SensProps.SensorRmax_allk{n}=sens.SensorRmax_allk;
            obj.SensProps.SensorDir_allk{n}=sens.SensorDir_allk;
            obj.SensProps.SensorAlpha_allk{n}=sens.SensorAlpha_allk;
            obj.SensProps.SensorTaskTargets_allk{n}=cell(obj.SimProps.Time.nsteps,1);
            obj.SensProps.SensorModeConstraints_allk{n}=sens.SensorModeConstraints_allk;
            obj.SensProps.SensorOnOff_allk{n}=sens.SensorOnOff_allk;
            
            
            obj.SensProps.SensorModelDim_allk{n}=sens.SensorModelDim_allk;
            
            
            %             obj.SensType(obj.Ns,1)=find(strcmpi( obj.SensTypes,SensorType)==1);%(Move/Station,Sensor Model)
            %             obj.SensType(obj.Ns,2)=find(strcmpi( obj.AvailSensors,SensorModel)==1);
            %             obj.hn(obj.Ns)=obj.AvailSensors_dim{find(strcmpi( obj.AvailSensors,SensorModel)==1)};
            %
            %             
            %             obj.sk{obj.Ns}(obj.tk,:)=sk;
            %             obj.alphak{obj.Ns}(obj.tk)=alphak;
            %             obj.rmaxk{obj.Ns}(obj.tk)=rmaxk;
            %             obj.dirk{obj.Ns}(obj.tk)=dirn;
            %
            %             obj.SensState{obj.Ns}(obj.tk)=1; % has to make measurement
            %
            %             obj.R{obj.Ns}=R;
            %             obj.G(obj.Ns)=find(strcmpi( obj.FOVpenaltytypes,FOVpenaltytype)==1);
            %
            %             obj.SensConstraints{obj.Ns}=find(strcmpi( obj.SensConstraintTypes,SensorConstraints)==1);
            
        end
        %============================================================================
        function obj=SwitchOff_MovingSensors_fromTk(obj,Tk)
            
            for nsens=1:obj.SensProps.NumbSensors
                for tt=Tk:1:obj.SimProps.Time.nsteps
                    if strcmpi(obj.SensProps.SensorDynTypes_allk{nsens}{tt},'Move')
                        obj.SensProps.SensorOnOff_allk{nsens}{tt}='Off';
                    end
                end
            end
            
        end
        %============================================================================
        function obj=UpdateSensorConstraints(obj,SENSCONSTRAINTS)
            %S={[1,2];[1,3,4];...} ... the sensorids....  each row of cell is one constraint
            % use this to update the Amat, bmatm Cmat and dmat for the
            % sensor constraits used in the optimization problem.
            % These constraints remain constant over time and are copied to
            % all time steps from current time step tk to final step.
            % So in the optimization problem, we just have to append these
            % matrics with any additonal constraints
            
            %             n=obj.SimProps.Time.tk;
            % Ax<=b, Cx==d
            if max(strcmpi(fieldnames(SENSCONSTRAINTS),'Seq'))==0 || max(strcmpi(fieldnames(SENSCONSTRAINTS),'deq'))==0  || max(strcmpi(fieldnames(SENSCONSTRAINTS),'Sineq'))==0  || max(strcmpi(fieldnames(SENSCONSTRAINTS),'dineq'))==0
                error('ok dude we have an error with given field names')
            end
            Seq=SENSCONSTRAINTS.Seq;
            deq=SENSCONSTRAINTS.deq;
            Sineq=SENSCONSTRAINTS.Sineq;
            dineq=SENSCONSTRAINTS.dineq;
            
            
            ns=obj.SensProps.NumbSensors;

            for i=1:1:size(Sineq,1)
                z=zeros(1,ns);
                c=Sineq{i};
                for j=c(:)'
                    z(find(obj.SensProps.SensorsID==j) )=1;
                end
                
                obj.SensProps.Amat=vertcat(obj.SensProps.Amat,z);
                obj.SensProps.bvec=vertcat(obj.SensProps.bvec,dineq{i});
            end
            
            for i=1:1:size(Seq,1)
                z=zeros(1,ns);
                c=Seq{i};
                for j=c(:)'
                    z(find(obj.SensProps.SensorsID==j) )=1;
                end
                
                obj.SensProps.Cmat=vertcat(obj.SensProps.Cmat,z);
                obj.SensProps.dvec=vertcat(obj.SensProps.dvec,deq{i});
            end
            
            
        end
        %============================================================================
        function obj=BackUpSensorConfig(obj)
            %first back up the previous
            obj.BackUpSimProps = obj.SimProps;
            obj.BackUpTargProps = obj.TargProps;
            obj.BackUpSensProps = obj.SensProps;
           
        end
        %============================================================================
        function obj=RestoreSensorConfig(obj)

            obj.SimProps = obj.BackUpSimProps;
            obj.TargProps = obj.BackUpTargProps;
            obj.SensProps = obj.BackUpSensProps;
  
        end
        %============================================================================
        function obj=BackUpTargMeanCovs(obj)
            obj.BackUp_targMeans = obj.TargProps.MeanState_allk;
            obj.BackUp_targCovss = obj.TargProps.CovState_allk;
        end
        %============================================================================
        function obj=RestoreUpTargMeanCovs(obj)
            
            obj.TargProps.MeanState_allk=obj.BackUp_targMeans ;
            obj.TargProps.CovState_allk=obj.BackUp_targCovss ;
        end

        %============================================================================
        function obj=StoreSensProps(obj)
            obj.BackUpSensProps = obj.SensProps;
        end
        %============================================================================
        function obj=ReStoreSensProps(obj)
            obj.SensProps = obj.BackUpSensProps;
        end
        %============================================================================
        function [y,H,G]=SensorMeasurement_pts(obj,Tk,X,sensid)
            % X are the target sigma points
            % tht esimulation with sensor sensid
            % for all the points in X at time Tk
            %             AA=obj.SensProps.SensorAlpha_allk{sensid};   %  GetContCell_Nt(obj.SensProps,'SensorAlpha_allk',sensid);
            %             SS=obj.SensProps.SensorPos_allk{sensid};   %  GetContCell_Nt(obj.SensProps,'SensorPos_allk',sensid);
            %             RR=obj.SensProps.SensorRmax_allk{sensid};   %  GetContCell_Nt(obj.SensProps,'SensorRmax_allk',sensid);
            %             DD=obj.SensProps.SensorDir_allk{sensid};   %  GetContCell_Nt(obj.SensProps,'SensorDir_allk',sensid);
            %
            %             HH=obj.SensProps.SensorModelDim_allk{sensid};   %  GetContCell_Nt(obj.SensProps,'SensorModelDim_allk',sensid);
            %
            
            alpha=obj.SensProps.SensorAlpha_allk{sensid}(Tk,1);
            Rmax=obj.SensProps.SensorRmax_allk{sensid}(Tk,1);
            dirn=obj.SensProps.SensorDir_allk{sensid}(Tk,1);
            xsenspos=obj.SensProps.SensorPos_allk{sensid}(Tk,:);
            
            hnn=obj.SensProps.SensorModelDim_allk{sensid}(Tk,1);
            
            y=zeros(size(X,1),hnn);
            G=zeros(size(X,1),1);
            H=ones(size(X,1),1); % ones means inside
            
            for i=1:1:size(X,1)
                
                r=norm([X(i,1),X(i,2)]-[xsenspos(1),xsenspos(2)]);
                th=atan2(X(i,2)-xsenspos(2),X(i,1)-xsenspos(1));
                %alpha is (-pi,pi)
                %dirn is (-pi,pi)
                %th is (-pi,pi)
                
                diff=dirn-th;
                if diff>pi
                    diff=diff-2*pi;
                end
                if diff<-pi
                    diff=diff+2*pi;
                end
                Model=obj.SensProps.SensorModel_allk{sensid};
                %Model=GetContCell_Nt(obj.SensProps,'SensorModel_allk',sensid);
                if strcmpi(Model{Tk},'Range+Bearing')
                    y(i,:)=[r,th];
                elseif strcmpi(Model{Tk},'X+Y')
                    y(i,:)=[X(i,1),X(i,2)];
                elseif strcmpi(Model{Tk},'Range+X')
                    y(i,:)=[r,X(i,2)];    
                elseif strcmpi(Model{Tk},'Range')
                    y(i,:)=r;
                elseif strcmpi(Model{Tk},'Bearing')
                    y(i,:)=th;
                    
                else
                    error('Cannot find sensor model')
                end
                
                G(i)=0;
                P=obj.SensProps.SensorFOVpenalty_allk{sensid}; %  GetContCell_Nt(obj.SensProps,'SensorFOVpenalty_allk',sensid);
                %                 P=GetContCell_Nt(obj.SensProps,'SensorFOVpenalty_allk',sensid);
                
                if r>Rmax
                    if strcmpi(P{Tk},'Simple')
                        G(i)=200;
                    elseif strcmpi(P{Tk},'Quadratic')
                        % quadratic penalty
                        G(i)=5*(r-Rmax)^2;
                    elseif strcmpi(P{Tk},'None')   % no penalty
                        G(i)=1;
                    end
                    H(i)=0;  % 0 means the reading is outside the FOV
                end
                if abs(diff)>alpha
                    if strcmpi(P{Tk},'Simple')
                        % simple constant penalty
                        G(i)=G(i)+200;
                    elseif strcmpi(P{Tk},'Quadratic')
                        G(i)=G(i)+5*(abs(diff)-alpha)^2;
                    elseif strcmpi(P{Tk},'None')   % no penalty
                        G(i)=1;
                    end
                    H(i)=0;  % 0 means the reading is outside the FOV
                end
                if abs(diff)<=alpha && r<=Rmax
                    G(i)=1;
                end
            end
            
        end
        %============================================================================
        function obj=CopySensorProps2NextTime(obj,Tk,Tk1,sensid)
            % copy sensor properties at Tk to Tk1 without modifying the tk value
            
            
            
            
            sensid=sensid(:)';
            for ss=sensid
                obj.SensProps.SensorAlpha_allk{ss}(Tk1,1)=obj.SensProps.SensorAlpha_allk{ss}(Tk,1);
                obj.SensProps.SensorRmax_allk{ss}(Tk1,1)=obj.SensProps.SensorRmax_allk{ss}(Tk,1);
                obj.SensProps.SensorDir_allk{ss}(Tk1,1)=obj.SensProps.SensorDir_allk{ss}(Tk,1);
                obj.SensProps.SensorPos_allk{ss}(Tk1,1)=obj.SensProps.SensorPos_allk{ss}(Tk,1);
                obj.SensProps.SensorModelDim_allk{ss}(Tk1,1)=obj.SensProps.SensorModelDim_allk{ss}(Tk,1);
                
            end
        end
        %============================================================================
        function obj=UpdateSensorProps_Tk(obj,Tk,sensid,config)
            % Update sensor properties at Tk
            %              config= {'Move',[100,123.45];
            %                       'MoveBy',[10,10];
            %                       'RotateToTarget',Targetid;
            %                         'RotateToAngle',pi/6';
            %                       'RotateByAngle',pi/6;
            %                       'DepthFOV',100;
            %                       'ApertureFOV',p/10}
            
       
            for i=1:1:size(config,1);
                ss=config{i,1};  % string
                V=config{i,2};   % its values
                
                if strcmpi(ss,'Move')
                    obj.SensProps.SensorPos_allk{sensid}(Tk,1)=V;
                    
                elseif strcmpi(ss,'MoveBy')
                    obj.SensProps.SensorPos_allk{sensid}(Tk,1)=obj.SensProps.SensorPos_allk{sensid}(Tk,1)+V;
                    
                elseif strcmpi(ss,'RotateToTarget')
                    targpos=obj.TargProps.MeanState_allk{V}(Tk,1:2);
                    senspos=obj.SensProps.SensorPos_allk{sensid}(Tk,1:2);
                    dr=targpos-senspos;
                    dr=dr/norm(dr);
                    a=atan2(dr(2),dr(1));
                    
                    obj.SensProps.SensorDir_allk{sensid}(Tk,1)=a;
                
                elseif strcmpi(ss,'RestoreBackupDirn')
                    obj.SensProps.SensorDir_allk{sensid}(Tk,1)=obj.BackUpSensProps.SensorDir_allk{sensid}(Tk,1);
                    
                elseif strcmpi(ss,'RotateByAngle')
                    obj.SensProps.SensorDir_allk{ss}(Tk,1)=obj.SensProps.SensorDir_allk{ss}(Tk,1)+V;
                elseif strcmpi(ss,'RotateToAngle')
                    obj.SensProps.SensorDir_allk{ss}(Tk,1)=V;
                    
                end
                
                
                % Moving on a grid code
                
                
                
            end
        end
        %============================================================================
        function obj=UpdateTimek(obj,TT)
            % You have to manually ipdate the time
            obj.SimProps.Time.tk=obj.SimProps.Time.tk+1;
            if TT==obj.SimProps.Time.tk
                disp('internal time and external time are nicely in sync')
            else
               error('internal time and external time not in sync') 
            end
        end
        %============================================================================
        function obj=SimulationModeOption(obj,Mode)
            
            
            %  Then each sensor can measure all the targets 
            % no FOV constraints during tasking and measurement process
            if Mode==-1
                obj.SimProps.SimMode=-1;
                return
            end
            
            
            
            % Default mode is 0, During tasking, if 3/4th of sigma points are out of FOV: then neglect object
            %                    During Sensor reading and measurement update: hell with sensor FOV constraints
            %
            if Mode==0
                obj.SimProps.SimMode=0;
                return
            end
            % Special mode 1:  During tasking, if 3/4th of sigma points are out of FOV: then neglect object
            %                  Respect the FOV while measurememnt update i.e. 3/4 rule sigma points
            if Mode==1
                obj.SimProps.SimMode=1;
                return
            end
            % Special mode 2:  During tasking, if 3/4th of sigma points are out of FOV: then neglect object
            %                   Respect the FOV while measurememnt update i.e. 3/4 rule sigma points
            %                   Additionally respect the fact that, if
            %                   sensor cannot see the truth then do not do
            %                   measurement update
            if Mode==2
                obj.SimProps.SimMode=2;
                return
            end
            
             % Special mode 0 + use penalty  
            if Mode==3
                obj.SimProps.SimMode=2;
                return
            end
            % Special mode 1 + use penalty  
            if Mode==4
                obj.SimProps.SimMode=2;
                return
            end
            % Special mode 2 + use penalty  
            if Mode==5
                obj.SimProps.SimMode=2;
                return
            end
            
        end
        %============================================================================
        function obj=DynProp_Target(obj,targs,Tk)
            % prop grom Tk to Tk+1
            %             %the time variable is not updated
            
            %first remove the targets that have been deleted i.e. do
            %nothing to them
            %             Tk=obj.SimProps.Time.tk;
            for targid=targs(:)'
                
                if Tk==0
                    continue;
                end
                
                if strcmpi(obj.TargProps.VisibilityState_allk{targid}{Tk},'NotTracked')
                    continue;
                end

                if strcmpi(obj.TargProps.TargetDynTypes_allk{targid}{Tk},'Stationary')  || strcmpi(obj.TargProps.TargetDynTypes_allk{targid}{Tk},'VirtualStationary')  
                    Xfk=obj.TargProps.MeanState_allk{targid}(Tk,:)';
                    Pfk=reshape(obj.TargProps.CovState_allk{targid}(Tk,:),obj.TargProps.TargetModelDim_allk{targid}(Tk),obj.TargProps.TargetModelDim_allk{targid}(Tk));
                    if mean(eig(sqrtm(Pfk)))>obj.TargProps.CovLimit{targid}
                        Pk1=Pfk;
                    else
                    Pk1=Pfk+obj.TargProps.Q_allk{targid}{Tk};
                    end
                    
                    obj.TargProps.MeanState_allk{targid}(Tk+1,:)=Xfk;
                    obj.TargProps.CovState_allk{targid}(Tk+1,:)=reshape(Pk1,1,obj.TargProps.TargetModelDim_allk{targid}(Tk)^2);
                    continue;
                end
                
                Xfk=obj.TargProps.MeanState_allk{targid}(Tk,:)';
                Pfk=reshape(obj.TargProps.CovState_allk{targid}(Tk,:),obj.TargProps.TargetModelDim_allk{targid}(Tk),obj.TargProps.TargetModelDim_allk{targid}(Tk));
                
                [Xk,wk]=obj.SimProps.SigmaPntFunction(Xfk,Pfk);
                Xk1=zeros(size(Xk));
                for i=1:1:length(wk)
                    switch(obj.TargProps.TargetModel_allk{targid}{Tk}) % type of dynamics
                        case 'UM'  % UM motion
                            Xk1(i,:)=KIRB_UM_eg_dyn_disc(Xk(i,:)',obj.SimProps.Time_dt);
                        case 'CT'  % CT motion
                            Xk1(i,:)=KIRB_CT_eg_dyn_disc(Xk(i,:)',obj.SimProps.Time_dt);
                        case 'NoDyn'  % No dynamics case
                            Xk1(i,:)=Xk(i,:);
                        case 'diffdyn'
                            m0=obj.TargProps.TargetModelParas_allk{targid}{Tk,1}.m0;
                            cf=obj.TargProps.TargetModelParas_allk{targid}{Tk,1}.cf;
                            Xk1(i,:)=diffdyn(Xk(i,:),m0,cf);
                    end
                end
                
                [xk1,Pk1]=MeanCov(Xk1,wk);
                
                if mean(eig(sqrtm(Pk1)))>obj.TargProps.CovLimit{targid}
                    Pk1=Pfk;
                else
                    Pk1=Pk1+obj.TargProps.Q_allk{targid}{Tk};
                end
                
                
                obj.TargProps.MeanState_allk{targid}(Tk+1,:)=xk1;
                obj.TargProps.CovState_allk{targid}(Tk+1,:)=reshape(Pk1,1,obj.TargProps.TargetModelDim_allk{targid}(Tk)^2);
            end
        end
        %============================================================================
        function obj=MeasUpdate_Target(obj,targs,sensids,Tk)
            %             do the measurement update for specified targets at curr time
            %             tk, with out modifying the time variable
            % RMEMEBER SENS-TASKING HAS TO BE DONE
            % SensTargtask{tk}(sensors,target)
            
%             Tk=obj.SimProps.Time.tk;

            TYPE='MeasurementUpdate';
%             TYPE='TASKING';

            for targid=targs(:)'
                %Do nothing to the deleted targets
                if strcmpi(obj.TargProps.VisibilityState_allk{targid}{Tk},'NotTracked')
                    continue;
                end
                
                % obj.SensProps.SensorStationaryTasks_allk=cell(1,1); %each cell is a Matrix with rows= sensor, columns= targets, cell rows are time step
                % obj.SensProps.SensorMoveTasks_allk=cell(1,1);
                %                 STATsensind=find(obj.SensProps.SensorStationaryTasks_allk{Tk}(:,2)==targid);
                %                 MOVEsensind=find(obj.SensProps.SensorMoveTasks_allk{Tk}(:,2)==targid);
                
                
                
                
                
                Xfk=obj.TargProps.MeanState_allk{targid}(Tk,:)';
                Pfk=reshape(obj.TargProps.CovState_allk{targid}(Tk,:),obj.TargProps.TargetModelDim_allk{targid}(Tk),obj.TargProps.TargetModelDim_allk{targid}(Tk));
                
                [Xk,wk]=obj.SimProps.SigmaPntFunction(Xfk,Pfk);
                RR=[];
                YM=[];
                Y=[];
                
                
                for nsens=sensids
     
                    [yy,H,G]=SensorMeasurement_pts(obj,Tk,Xk,nsens);
                    
                    if strcmpi(obj.SensProps.SensorOnOff_allk{nsens}{Tk},'Off') 
                        continue; 
                    end
                    
                    if obj.SimProps.SimMode==-1  % no FOV anywhere
                        [ym,Hm,Gm]=SensorMeasurement_pts(obj,Tk,obj.TargProps.TruePosState_allk{targid}(Tk,:),nsens);
                        ym=ym(:);
                        ym=ym+sqrtm(obj.SensProps.R_allk{nsens}{Tk})*randn(obj.SensProps.SensorModelDim_allk{nsens}(Tk),1);
                        obj.SensProps.y_allk{Tk}{targid,nsens}=ym;
                        YM=vertcat(YM,ym);
                        Y=horzcat(Y,yy);
                        RR=blkdiag(RR,obj.SensProps.R_allk{nsens}{Tk});
                        
                    elseif obj.SimProps.SimMode==0 % During Tasking: Use FOV, During Measurement : No FOV
                        
                        [ym,Hm,Gm]=SensorMeasurement_pts(obj,Tk,obj.TargProps.TruePosState_allk{targid}(Tk,:),nsens);
                        ym=ym(:);
                        ym=ym+sqrtm(obj.SensProps.R_allk{nsens}{Tk})*randn(obj.SensProps.SensorModelDim_allk{nsens}(Tk),1);
                        obj.SensProps.y_allk{Tk}{targid,nsens}=ym;
                        YM=vertcat(YM,ym);
                        Y=horzcat(Y,yy);
                        RR=blkdiag(RR,obj.SensProps.R_allk{nsens}{Tk});
                        
                    elseif obj.SimProps.SimMode==1 % During Tasking: Use FOV, During Measurement : Use FOV, but real measurements can be outside of FOV
                        
                        if length(find(H==0))>=3/4*length(H) 
                            continue; % if no sig pt is visible just dont use this sensor
                        end
                    
                        [ym,Hm,Gm]=SensorMeasurement_pts(obj,Tk,obj.TargProps.TruePosState_allk{targid}(Tk,:),nsens);
                        ym=ym(:);
                        ym=ym+sqrtm(obj.SensProps.R_allk{nsens}{Tk})*randn(obj.SensProps.SensorModelDim_allk{nsens}(Tk),1);
                        obj.SensProps.y_allk{Tk}{targid,nsens}=ym;
                        YM=vertcat(YM,ym);
                        Y=horzcat(Y,yy);
                        RR=blkdiag(RR,obj.SensProps.R_allk{nsens}{Tk});
                        
                        
                    elseif obj.SimProps.SimMode==2  % During Tasking: Use FOV, During Measurement : Use FOV, real measurements have to be inside of FOV
                        
                        if length(find(H==0))>=3/4*length(H) 
                            continue; % if no sig pt is visible just dont use this sensor
                        end
                        
                        [ym,Hm,Gm]=SensorMeasurement_pts(obj,Tk,obj.TargProps.TruePosState_allk{targid}(Tk,:),nsens);
                        ym=ym(:);
                        ym=ym+sqrtm(obj.SensProps.R_allk{nsens}{Tk})*randn(obj.SensProps.SensorModelDim_allk{nsens}(Tk),1);
                        
                        if Hm==0  % if u cannot see the truth dont use the measuremetns===
                            disp(strcat('Target ',num2str(targid),' though tasked by sensor ',num2str(nsens),' cannot be updated for real' ))
                            continue
                        end
                        
                        obj.SensProps.y_allk{Tk}{targid,nsens}=ym;
                        YM=vertcat(YM,ym);
                        Y=horzcat(Y,yy);
                        RR=blkdiag(RR,obj.SensProps.R_allk{nsens}{Tk});
                        
                    else
                        error('Simulation mode unknown')
                    end
                    
                    
                    
                    
                    
                    
                    
                    

                    
                    
                    
                end
                
                
                
                
                
                
                % Stack all the measurmeents together for the Kalman Update
                
                if isempty(Y)==0 && isempty(YM)==0 && size(Y,2)==size(YM,1)
                    [mz,Pz]=MeanCov(Y,wk);
                    Pz=Pz+RR;
                    Pcc=CrossCov(Xk,Xfk,Y,mz,wk);
                    [xku,Pku]=KalmanUpdate(Xfk,Pfk,mz,Pz,Pcc,YM);
                else
%                     disp(' Some True measurements not takens  OR sigma points outside FOV')
                    xku=Xfk;
                    Pku=Pfk;
                end
                obj.TargProps.MeanState_allk{targid}(Tk,:)=xku;
                obj.TargProps.CovState_allk{targid}(Tk,:)=reshape(Pku,1,obj.TargProps.TargetModelDim_allk{targid}(Tk)^2);
            end
            
        end
        %============================================================================
        function [Pfk,Pku, Xfk,xku]=MeasUpdate_Pseudo(obj,targid,sensid,Tk)
            % ONLY SINGLE TARGET AS ONLY 1 OUTPUT CAN BE SENT OUT
            
            %             TYPE='MeasurementUpdate';
            TYPE='TASKING';
            
            
            Xfk=obj.TargProps.MeanState_allk{targid}(Tk,:)';
            Pfk=reshape(obj.TargProps.CovState_allk{targid}(Tk,:),obj.TargProps.TargetModelDim_allk{targid}(Tk),obj.TargProps.TargetModelDim_allk{targid}(Tk));
            
            if strcmpi(obj.TargProps.VisibilityState_allk{targid}{Tk},'NotTracked')
                Pku=Pfk;
                xku=Xfk;
                return;
            end
            
            [Xk,wk]=obj.SimProps.SigmaPntFunction(Xfk,Pfk);
            RR=[];
            Y=[];
            
            
            % if sensid is empty
            if isempty(sensid)
                SS= 1:obj.SensProps.NumbSensors;
            else
                SS=sensid(:)';
            end
            
            for nsens=SS
                
                [yy,H,G]=SensorMeasurement_pts(obj,Tk,Xk,nsens);
                
                if strcmpi(obj.SensProps.SensorOnOff_allk{nsens}{Tk},'Off')
                    continue;
                end
                
                if obj.SimProps.SimMode==-1  % no FOV anywhere
                    
                    Y=horzcat(Y,yy);
                    RR=blkdiag(RR,obj.SensProps.R_allk{nsens}{Tk});
                    
                elseif obj.SimProps.SimMode==0 % During Tasking: Use FOV, During Measurement : No FOV
                    
                    
                    Y=horzcat(Y,yy);
                    RR=blkdiag(RR,obj.SensProps.R_allk{nsens}{Tk});
                    
                elseif obj.SimProps.SimMode==1 % During Tasking: Use FOV, During Measurement : Use FOV, but real measurements can be outside of FOV
                    
                    if length(find(H==0))>=3/4*length(H)
                        continue; % if no sig pt is visible just dont use this sensor
                    end
                    
                    
                    Y=horzcat(Y,yy);
                    RR=blkdiag(RR,obj.SensProps.R_allk{nsens}{Tk});
                    
                    
                elseif obj.SimProps.SimMode==2  % During Tasking: Use FOV, During Measurement : Use FOV, real measurements have to be inside of FOV
                    
                    if length(find(H==0))>=3/4*length(H)
                        continue; % if no sig pt is visible just dont use this sensor
                    end
                    
                    
                    Y=horzcat(Y,yy);
                    RR=blkdiag(RR,obj.SensProps.R_allk{nsens}{Tk});
                    
                else
                    error('Simulation mode unknown')
                end
            end
            
            
            
            
            
            
            % Stack all the measurmeents together for the Kalman Update
            
            if isempty(Y)==0
                [mz,Pz]=MeanCov(Y,wk);
                Pz=Pz+RR;
                Pcc=CrossCov(Xk,Xfk,Y,mz,wk);
                [xku,Pku]=KalmanUpdate(Xfk,Pfk,mz,Pz,Pcc,mz);
            else
%                 disp(' Some True measurements not takens  OR sigma points outside FOV')
                xku=Xfk;
                Pku=Pfk;
            end
            
            
            %             end
        end
        %============================================================================
        function obj=MeasUpdate_Pseudo_loadit(obj,targs,sensid,Tk)
            %  Do simple pseudo measurement update to get the updated mean
            %SINGLE TARGET ONLY
            % and covarianr for tasking purposes
            % DO NOT CHANGE TIME OR UPDATE OBJ states
%             Tk=obj.SimProps.Time.tk;

%             TYPE='MeasurementUpdate';
            TYPE='TASKING';
            
            for targid=targs(:)'
                %Do nothing to the deleted targets
                
                
                % obj.SensProps.SensorStationaryTasks_allk=cell(1,1); %each cell is a Matrix with rows= sensor, columns= targets, cell rows are time step
                % obj.SensProps.SensorMoveTasks_allk=cell(1,1);
                %                 STATsensind=find(obj.SensProps.SensorStationaryTasks_allk{Tk}(:,2)==targid);
                %                 MOVEsensind=find(obj.SensProps.SensorMoveTasks_allk{Tk}(:,2)==targid);
                
                
                if strcmpi(obj.TargProps.VisibilityState_allk{targid}{Tk},'NotTracked')
                    continue
                end
                
                
                Xfk=obj.TargProps.MeanState_allk{targid}(Tk,:)';
                Pfk=reshape(obj.TargProps.CovState_allk{targid}(Tk,:),obj.TargProps.TargetModelDim_allk{targid}(Tk),obj.TargProps.TargetModelDim_allk{targid}(Tk));
                
                
                
                [Xk,wk]=obj.SimProps.SigmaPntFunction(Xfk,Pfk);
                
                
                RR=[];
                Y=[];
                
                if isempty(sensid)
                    SS= 1:obj.SensProps.NumbSensors;
                else
                    SS=sensid(:)';
                end
                
                for nsens=SS
                    
                    
                    
                   [yy,H,G]=SensorMeasurement_pts(obj,Tk,Xk,nsens);
                    
                    if strcmpi(obj.SensProps.SensorOnOff_allk{nsens}{Tk},'Off') 
                        continue; 
                    end
                    
                    if obj.SimProps.SimMode==-1  % no FOV anywhere
                       
                        Y=horzcat(Y,yy);
                        RR=blkdiag(RR,obj.SensProps.R_allk{nsens}{Tk});
                        
                    elseif obj.SimProps.SimMode==0 % During Tasking: Use FOV, During Measurement : No FOV
                        
                       
                        Y=horzcat(Y,yy);
                        RR=blkdiag(RR,obj.SensProps.R_allk{nsens}{Tk});
                        
                    elseif obj.SimProps.SimMode==1 % During Tasking: Use FOV, During Measurement : Use FOV, but real measurements can be outside of FOV
                        
                        if length(find(H==0))>=3/4*length(H) 
                            continue; % if no sig pt is visible just dont use this sensor
                        end
                    
                       
                        Y=horzcat(Y,yy);
                        RR=blkdiag(RR,obj.SensProps.R_allk{nsens}{Tk});
                        
                        
                    elseif obj.SimProps.SimMode==2  % During Tasking: Use FOV, During Measurement : Use FOV, real measurements have to be inside of FOV
                        
                        if length(find(H==0))>=3/4*length(H) 
                            continue; % if no sig pt is visible just dont use this sensor
                        end
                    
                        
                        Y=horzcat(Y,yy);
                        RR=blkdiag(RR,obj.SensProps.R_allk{nsens}{Tk});
                        
                    else
                        error('Simulation mode unknown')
                    end
                    
                    
                end
                
                
                
                
                
                
                % Stack all the measurmeents together for the Kalman Update
                
                if isempty(Y)==0
                    [mz,Pz]=MeanCov(Y,wk);
                    Pz=Pz+RR;
                    Pcc=CrossCov(Xk,Xfk,Y,mz,wk);
                    [xku,Pku]=KalmanUpdate(Xfk,Pfk,mz,Pz,Pcc,mz);
                else
%                     disp(' Some True measurements not takens  OR sigma points outside FOV')
                    xku=Xfk;
                    Pku=Pfk;
                end
                obj.TargProps.MeanState_allk{targid}(Tk,:)=xku;
                obj.TargProps.CovState_allk{targid}(Tk,:)=reshape(Pku,1,obj.TargProps.TargetModelDim_allk{targid}(Tk)^2);
                
            end
        end
        %============================================================================
        function [muk,Pk]=MeasUpdate_FDP(obj,muk,Pk,targs,sensid,Tk)
            %  Do simple pseudo measurement update to get the updated mean
            %SINGLE TARGET ONLY
            % and covarianr for tasking purposes
            % DO NOT CHANGE TIME OR UPDATE OBJ states
            
            % muk = obj.TargProps.MeanState_allk
            % Pk = obj.TargProps.CovState_allk
            
            
%             TYPE='MeasurementUpdate';
            TYPE='TASKING';
            
            for targid=targs(:)'
                if strcmpi(obj.TargProps.VisibilityState_allk{targid}{Tk},'NotTracked')
                    continue
                end
                
                
                Xfk=muk{targid}(Tk,:)';
                Pfk=reshape(Pk{targid}(Tk,:),obj.TargProps.TargetModelDim_allk{targid}(Tk),obj.TargProps.TargetModelDim_allk{targid}(Tk));
                
                
                [Xk,wk]=obj.SimProps.SigmaPntFunction(Xfk,Pfk);
                RR=[];
                Y=[];
                
                if isempty(sensid)
                    SS= 1:obj.SensProps.NumbSensors;
                else
                    SS=sensid(:)';
                end
                
                for nsens=SS
                    
                    [yy,H,G]=SensorMeasurement_pts(obj,Tk,Xk,nsens);
                    
                    if strcmpi(obj.SensProps.SensorOnOff_allk{nsens}{Tk},'Off') 
                        continue; 
                    end
                    
                    if obj.SimProps.SimMode==-1  % no FOV anywhere
                      
                        Y=horzcat(Y,yy);
                        RR=blkdiag(RR,obj.SensProps.R_allk{nsens}{Tk});
                        
                    elseif obj.SimProps.SimMode==0 % During Tasking: Use FOV, During Measurement : No FOV
                        
                       
                        Y=horzcat(Y,yy);
                        RR=blkdiag(RR,obj.SensProps.R_allk{nsens}{Tk});
                        
                    elseif obj.SimProps.SimMode==1 % During Tasking: Use FOV, During Measurement : Use FOV, but real measurements can be outside of FOV
                        
                        if length(find(H==0))>=3/4*length(H) 
                            continue; % if no sig pt is visible just dont use this sensor
                        end
                    
                        
                        Y=horzcat(Y,yy);
                        RR=blkdiag(RR,obj.SensProps.R_allk{nsens}{Tk});
                        
                        
                    elseif obj.SimProps.SimMode==2  % During Tasking: Use FOV, During Measurement : Use FOV, real measurements have to be inside of FOV
                        
                        if length(find(H==0))>=3/4*length(H) 
                            continue; % if no sig pt is visible just dont use this sensor
                        end
                    
                        
                        Y=horzcat(Y,yy);
                        RR=blkdiag(RR,obj.SensProps.R_allk{nsens}{Tk});
                        
                    else
                        error('Simulation mode unknown')
                    end
                    
                    
                end
                
                % Stack all the measurmeents together for the Kalman Update
                
                if isempty(Y)==0 
                    [mz,Pz]=MeanCov(Y,wk);
                    Pz=Pz+RR;
                    Pcc=CrossCov(Xk,Xfk,Y,mz,wk);
                    [xku,Pku]=KalmanUpdate(Xfk,Pfk,mz,Pz,Pcc,mz);
                else
%                     disp(' Some True measurements not takens  OR sigma points outside FOV')
                    xku=Xfk;
                    Pku=Pfk;
                end
                
                muk{targid}(Tk,:)=xku;
                Pk{targid}(Tk,:)=reshape(Pku,1,obj.TargProps.TargetModelDim_allk{targid}(Tk)^2);
                
            end
        end
        %============================================================================
        function [muk,Pk]=DynProp_FDP(obj,muk,Pk,targs,Tk)
            % prop grom Tk to Tk+1
            %             %the time variable is not updated
            
            %first remove the targets that have been deleted i.e. do
            %nothing to them
            %             Tk=obj.SimProps.Time.tk;
            for targid=targs(:)'
                if Tk==0
                    continue;
                end
                
                if strcmpi(obj.TargProps.VisibilityState_allk{targid}{Tk},'NotTracked')
                    continue;
                end
                
                if strcmpi(obj.TargProps.TargetDynTypes_allk{targid}{Tk},'Stationary')  || strcmpi(obj.TargProps.TargetDynTypes_allk{targid}{Tk},'VirtualStationary') 
                    Xfk=muk{targid}(Tk,:)';
                    Pfk=reshape(Pk{targid}(Tk,:),obj.TargProps.TargetModelDim_allk{targid}(Tk),obj.TargProps.TargetModelDim_allk{targid}(Tk));
                    if mean(eig(sqrtm(Pfk)))>obj.TargProps.CovLimit{targid}
                        Pk1=Pfk;
                    else
                        pp=obj.TargProps.Q_allk{targid};
                        Pk1=Pfk+pp{Tk};
                    end
                    
                    muk{targid}(Tk+1,:)=Xfk;
                    Pk{targid}(Tk+1,:)=reshape(Pk1,1,obj.TargProps.TargetModelDim_allk{targid}(Tk)^2);
                    continue;
                end
                
                Xfk=muk{targid}(Tk,:)';
                Pfk=reshape(Pk{targid}(Tk,:),obj.TargProps.TargetModelDim_allk{targid}(Tk),obj.TargProps.TargetModelDim_allk{targid}(Tk));
                
                [Xk,wk]=obj.SimProps.SigmaPntFunction(Xfk,Pfk);
                Xk1=zeros(size(Xk));
                 for i=1:1:length(wk)
                    switch(obj.TargProps.TargetModel_allk{targid}{Tk}) % type of dynamics
                        case 'UM'  % UM motion
                            Xk1(i,:)=KIRB_UM_eg_dyn_disc(Xk(i,:)',obj.SimProps.Time_dt);
                        case 'CT'  % CT motion
                            Xk1(i,:)=KIRB_CT_eg_dyn_disc(Xk(i,:)',obj.SimProps.Time_dt);
                        case 'NoDyn'  % No dynamics case
                            Xk1(i,:)=Xk(i,:);
                        case 'diffdyn'
                            m0=obj.TargProps.TargetModelParas_allk{targid}{Tk,1}.m0;
                            cf=obj.TargProps.TargetModelParas_allk{targid}{Tk,1}.cf;
                            Xk1(i,:)=diffdyn(Xk(i,:),m0,cf);
                    end
                end
                
                [xk1,Pk1]=MeanCov(Xk1,wk);
                
                if mean(eig(sqrtm(Pk1)))>obj.TargProps.CovLimit{targid}
                    Pk1=Pfk;
                else
                    Pk1=Pk1+obj.TargProps.Q_allk{targid}{Tk};
                end
                
                
                muk{targid}(Tk+1,:)=xk1;
                Pk{targid}(Tk+1,:)=reshape(Pk1,1,obj.TargProps.TargetModelDim_allk{targid}(Tk)^2);
            end
        end
        %============================================================================
        function index=Index_MUvec(obj,targid,sensid)
            % get the index number
            if size(obj.SimProps.MUindex,1)~=obj.TargProps.NumbTargets || size(obj.MUindex,2)~=obj.SensProps.NumbSensors % checking if this index is already created or not
                obj.SimProps.MUindex=zeros(obj.TargProps.NumbTargets,obj.SensProps.NumbSensors);   %[ rows are Target, columns are sensors]
                k=1;
                for j=1:1:obj.SensProps.NumbSensors
                    for i=1:1:obj.TargProps.NumbTargets
                        obj.SimProps.MUindex(i,j)=k;
                        k=k+1;
                    end
                end
            end
            if isempty(targid) %give all the indexes of the sensor sensid
                index=obj.SimProps.MUindex(:,sensid);
                return
            end
            if isempty(sensid) %give all the indexes of the sensor sensid
                index=obj.SimProps.MUindex(targid,:);
                return
            end
            index=obj.SimProps.MUindex(targid,sensid);
            
            
        end
        %============================================================================
        function obj=OptiTasking_MIUB_ALGO1(obj)
            % does for only stationary sensors
            % does not work for moving sensors
            % Do sensor tasking with upper bound i.e. binary int programmign
            obj.C=[];
            obj.A=[];
            obj.b=[];
            obj.d=[];
            II=eye(obj.TargProps.NumbTargets*obj.SensProps.NumbSensors);
            
            Tk=obj.SimProps.Time.tk;
            
            %             obj.SensProps.SensorFOVtype_allk=cell(1,1); %{'1Sensor->1Target','1Sensor->AllTarget','1Sensor->AllTarget/1Target'};
            % Ax <= b and Cx = d
            
            % adding constraints to deselect all the deleted target
            % variables
            for targid=1:1:obj.TargProps.NumbTargets
                if strcmpi(obj.TargProps.VisibilityState_allk{targid}{Tk},'NotTracked')
                    ind=Index_MUvec(obj,targid,[]);
                    c=zeros(1,obj.TargProps.NumbTargets*obj.SensProps.NumbSensors);
                    c(ind)=1;
                    obj.C=vertcat(obj.C,c);
                    obj.d=vertcat(obj.d,0);
                end
                
            end
            
            % 1 sensor can see only one target
            for sensid=1:1:obj.SensProps.NumbSensors
                if strcmpi(obj.SensProps.SensorFOVtype_allk{sensid}{Tk},'1Sensor->1Target')  % tracking sensor
                    ind=Index_MUvec(obj,[],sensid); % all the targets that are paired to this sensor
                    c=zeros(1,obj.TargProps.NumbTargets*obj.SensProps.NumbSensors);
                    c(ind)=1;
                    obj.A=vertcat(obj.A,c);
                    obj.b=vertcat(obj.b,1);
                end
                
            end
            % adding the mode cosntraints in obj.SensModeConstraintsVec
            % so that only 1 mode is selected
            %             obj.SensProps.Amat
            
            for i=1:1:size(obj.SensProps.Amat,1)
                V=find(obj.SensProps.Amat(i,:)==1); % V has the sensids that are related
                c=zeros(1,obj.TargProps.NumbTargets*obj.SensProps.NumbSensors);
                
                for sensid=V(:)'
                    ind=Index_MUvec(obj,[],sensid);
                    if strcmpi(obj.SensProps.SensorFOVtype_allk{sensid}{Tk},'1Sensor->1Target')  % tracking sensor
                        c(ind)=1;
                    elseif strcmpi(obj.SensProps.SensorFOVtype_allk{sensid}{Tk},'1Sensor->AllTarget')  % coverage sensor
                        c(ind)=1/obj.TargProps.NumbTargets;
                        % now adding that in coverage sensor all the
                        % targets have to be measured if u try to measure one
                        for j=1:1:length(ind)
                            dd=zeros(1,obj.TargProps.NumbTargets*obj.SensProps.NumbSensors);
                            dd(ind)=1/obj.TargProps.NumbTargets;
                            dd(ind(j))=dd(ind(j))-1;
                            obj.A=vertcat(obj.A,dd);
                            obj.b=vertcat(obj.b,0);
                        end
                    end
                end
                obj.A=vertcat(obj.A,c);
                obj.b=vertcat(obj.b,obj.SensProps.bvec(i));
                
                
            end
            
            %  Make BackUp of Sensor Config as we
            %             obj=RestoreSensorConfig(obj)
            %             obj=BackUpSensorConfig(obj);
            
            %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            
            disp('MIub optimization is selected')
            % first for all pairs find the MI
            
            F=zeros(obj.TargProps.NumbTargets,obj.SensProps.NumbSensors);
            for targid=1:1:obj.TargProps.NumbTargets
                for sensid=1:1:obj.SensProps.NumbSensors
                    
                    %                         obj=UpdateSensorProps_Tk(obj,sensid,{'RotateToTarget',targid});
                    
                    F(targid,sensid)=MutualInformationByChangingSensorConfig(obj,Tk,targid,sensid);
                    %                         F(targid,sensid)=MutualInformation(obj,targid,sensid);
                    %                         obj=RestoreSensorConfig(obj);
                    
                    if F(targid,sensid)==0  % remove this sensor by constraining to be 0
                        ind=Index_MUvec(obj,targid,sensid);
                        obj.C=vertcat(obj.C,II(ind,:));
                        obj.d=vertcat(obj.d,0);
                    end
                end
            end
            
            
            %now solve the bintprog problem with the constraints
        
            f=reshape(F,obj.TargProps.NumbTargets*obj.SensProps.NumbSensors,1);
            obj.FIM=F;
            obj.MUvec=bintprog(-f,obj.A,obj.b,obj.C,obj.d);
            obj.MU=reshape(obj.MUvec,obj.TargProps.NumbTargets,obj.SensProps.NumbSensors);
            
            % now load the solutoin into the SensTask variabe;
            %==============================================================
            for sensid=1:1:obj.SensProps.NumbSensors
                obj.SensProps.SensorTaskTargets_allk{sensid}{Tk}=[];
                for targid=1:1:obj.TargProps.NumbTargets
                    if obj.MU(targid,sensid)==1
                        obj.SensProps.SensorTaskTargets_allk{sensid}{Tk}=horzcat(obj.SensProps.SensorTaskTargets_allk{sensid}{Tk},targid);
                        % also update the direction of pointing of statioanry sensors
                        obj=UpdateSensorProps_Tk(obj,sensid,{'RotateToTarget',targid});
                        
                    end
                end
            end
            %==============================================================
            
            
            
        end
        %============================================================================
        function obj=CoverageSearchNewTargets(obj,P0newtarget,Qknewtarget)
            % use the list of new targets ( this taget list with its true trajectories are already know at start but not visible as they are far away)
            % once optimal tasking is done, use this optimal sensors configuration to
            % check if any new objects are there within the field of view and
            % this target to the list of known targets (till time tk append a stationary position estimate as [-100,-100])
            % the initial estimate is the current estimate, with some
            % covariance matrix (P0newtarget) same for all new targets
            % for now P0newtarget,Qknewtarget are of dim for CT model
            
            for nwtargid=1:1:length(obj.NEWxktruth)
                Y=cell(1,2);
                k=1;
                sensind=find(obj.SensTargtask{obj.tk}(:,2)==targid);
                sensids=obj.SensTargtask{obj.tk}(sensind,1);
                for nsens=sensids
                    [ym,H,G]=SensorMeasurement_pts(obj,obj.tk,obj.NEWxktruth{nwtargid}(obj.tk,:),nsens);
                    if length(find(H==0))==length(H)  %if no. of outofFOV sigpts == no. of sigpts
                        continue; % if no sig pt is visible just dont use this sensor
                    end
                    ym=ym+sqrtm(obj.R{nsens})*randn(obj.hn(nsens),1);
                    Y{k,1}=ym;
                    Y{k,2}=obj.AvailSensors(obj.SensType(nsens,2));% tell if the meas is range or bearing or both
                end
                %xx0=obj.GenerateStateFromMeasurement(Y);
                xx0=[obj.NEWxktruth{nwtargid}(obj.tk,1:2),1,1]';
                PP0=P0newtarget;
                QQk=Qknewtarget;
                obj=obj.Add_Target(obj,'Move','CT',xx0,PP0,obj.NEWxktruth{nwtargid},QQk);
            end
            
        end
        
        %============================================================================
        function MI=MutualInformation(obj,targid,sensid)
            
            [Pf,Pu, ~,~]=MeasUpdate_Pseudo(obj,targid,sensid,Tk);
            
            
            if det(Pf)/det(Pu)<1.001
                MI=0;
            else
                MI=0.5*log(det(Pf)/det(Pu));
            end
        end
        %============================================================================
        function MI=MutualInformationByChangingSensorConfig(obj,Tk,targid,sensid)
            if length(targid)>1
                error('Haha Haha !! Only 1 target at a time dude');
            end
          
            obj=BackUpSensorConfig(obj);
            for nsens=sensid
                obj=UpdateSensorProps_Tk(obj,Tk,nsens,{'RotateToTarget',targid});
            end
            [Pf,Pu, ~,~]=MeasUpdate_Pseudo(obj,targid,sensid,Tk);
            
            obj=RestoreSensorConfig(obj);
            
            
            
            
            if det(Pf)/det(Pu)<1.001
                MI=0;
            else
                MI=0.5*log(det(Pf)/det(Pu));
            end
        end
        
        %============================================================================
        function obj=BruteForceTasking(obj,Tk)
            % GOING TO FORCE TASK ALL SENSORS TO SEE EVERYTHING
            for k=Tk(:)'
                obj.SensProps.SensorTaskTargets_allk{k}=[];
                for nsens=1:obj.SensProps.NumbSensors
                    for targid=1:obj.TargProps.NumbTargets
                        obj.SensProps.SensorTaskTargets_allk{k}= vertcat(obj.SensProps.SensorTaskTargets_allk{k},[nsens,targid]);
                    end
                end
            end
        end
        %===========================================================================
        function [SensIDstat,SensIDmove]=GetMovingSensorIDs(obj,Tk)
            % GOING TO FORCE TASK ALL SENSORS TO SEE EVERYTHING
            SensIDmove=zeros(1,obj.SensProps.NumbSensors);
            for nsens=1:obj.SensProps.NumbSensors
                if strcmpi(obj.SensProps.SensorDynTypes_allk{nsens}{Tk},'Move')
                    SensIDmove(nsens)=nsens;
                end
            end
            SensIDmove(SensIDmove==0)=[];
            SensIDstat=setdiff([1:1:obj.SensProps.NumbSensors],SensIDmove);
        end
        
        %============================================================================
        function obj=SequentialExhaustiveUavPathPlanner(obj,Tk)
            % From Tk onwards do the path planning to TK1 th time step (inclduing Tk and Tk1)
            % plans all the sensor config properties for the given time steps
            % if not pplanned , the config is just copied over to the next
            % time steps
            
            Tvec=Tk:1:Tk+obj.SimProps.NoPlanningSteps ;
            
            % GOING TO FORCE TASK ALL SENSORS TO SEE EVERYTHING
            obj=BruteForceTasking(obj,Tvec);
            
            % First Copy the sensor properties from Tk all the way upto Tk1
            for k=Tvec(1:end-1)
                obj=CopySensorProps2NextTime(obj,k,k+1,1:obj.SensProps.NumbSensors);
            end
            
            % get the moving and stationary sensor ids
            [SensIDstat,SensIDmove]=GetMovingSensorIDs(obj,Tvec(1));
            
            
            %obj=BackUpTargMeanCovs(obj);  % Make MAster copy of the means and covariances
            
            % get the prior pdfs
            MMk=obj.TargProps.MeanState_allk;
            PPk=obj.TargProps.CovState_allk;
            
            for k=Tvec(1:end)  % -1 as it propagates from k to k+1 in dynprop function
                [MMk,PPk]=DynProp_FDP(obj,MMk,PPk,1:obj.TargProps.NumbTargets,k-1);
            end
            
            % get the 2-sigma grid nodes
            TargSensNodes=[];  
            for targid=1:obj.TargProps.NumbTargets
                for k=Tvec
                    Xfk=MMk{targid}(k,1:2)';
                    Pfk=reshape(PPk{targid}(k,:),obj.TargProps.TargetModelDim_allk{targid}(k),obj.TargProps.TargetModelDim_allk{targid}(k));
                    Pfk=Pfk(1:2,1:2);
                    XYellipse=sigmaEllipse(Xfk,Pfk,2);

                    IN = inpolygon(obj.SimProps.XY(:,1),obj.SimProps.XY(:,2),XYellipse(:,1),XYellipse(:,2));
                    TargSensNodes=horzcat(TargSensNodes,obj.SimProps.XYindex(IN)');
                end
            end
            TargSensNodes=unique(TargSensNodes);
            TargSensNodesPos=obj.SimProps.XY(TargSensNodes,:);
            
            %obj=RestoreUpTargMeanCovs(obj);
            
%             obj=BackUpTargMeanCovs(obj);
            

            sensids=SensIDstat;
            
            for movtarg=SensIDmove
                
                xsenspos0=obj.SensProps.SensorPos_allk{movtarg}(Tvec(1),:);
                [~,UAVnode0ind]=min(sqrt(sum((repmat(xsenspos0,size(obj.SimProps.XY,1),1)-obj.SimProps.XY).^2,2)));
                UAVnode0=obj.SimProps.XYindex(UAVnode0ind);

                StateTree=GenrateTreefromGrid(obj.SimProps.XY,obj.SimProps.XYADJindex,   obj.SimProps.NoPlanningSteps+1  ,xsenspos0);
   
                P=Generatepaths(StateTree,UAVnode0,1);
                
                TargSensNodes=setdiff(TargSensNodes,UAVnode0);
                Pp=P;
%                 Pp=PrunePaths(P,TargSensNodes);
%                 
%                 
%                 if isempty(Pp) % the nodes are not reachable
%                     
%                     %if the number of steps cannot reach the nearest ifo node,
%                     %then form a triangle [UAV, nearest node, farthest node],
%                     % prune trajectories to be in this triangle
%                     TargSensNodesPos=obj.SimProps.XY(TargSensNodes,:);
%                     dis=repmat(xsenspos0,size(TargSensNodesPos,1),1)-TargSensNodesPos;
%                     [~,mindis]=min(sqrt(sum(dis.^2,2)));
%                     [~,maxdis]=max(sqrt(sum(dis.^2,2)));
%                 
%                     nearestnode=TargSensNodes(mindis);
%                     farthestnode=TargSensNodes(maxdis);
%                     
%                     pruntriagregion=obj.SimProps.XY([UAVnode0ind,nearestnode,farthestnode],:);
%                     
%                     IN = inpolygon(obj.SimProps.XY(:,1),obj.SimProps.XY(:,2),pruntriagregion(:,1),pruntriagregion(:,2));
%                     pruntriagnodes=obj.SimProps.XYindex(IN)';
%                     
%                     Pp=PrunePaths(P,pruntriagnodes);
%                 end
                
                
                % Debug plotting
%                 keyboard 

%                 obj.plot_Dynamic_scenario(Tvec(1));
%                 hold on
%                 for jj=1:1:size(Pp,1)
%                     for ij=1:1:length(Tvec)
%                    plot(obj.SimProps.XY(Pp(jj,ij),1),obj.SimProps.XY(Pp(jj,ij),2),'ro' ) 
%                     end
%                 end
%                 for jj=1:1:length(TargSensNodes)
%                     plot(obj.SimProps.XY(TargSensNodes(jj),1),obj.SimProps.XY(TargSensNodes(jj),2),'bo','MarkerSize',8 ) 
%                 end
%                 plot3(obj.SimProps.XY(:,1),obj.SimProps.XY(:,2),obj.SimProps.XYindex,'.')
%                 
%                 
                
                
                % METHOD 1:   doing sequential estimation while taking all the static
                % sensors
                
                
                Mi_paths=zeros(size(Pp,1),1);
                
                sensids=union( sensids , movtarg); % Sequentially add the moving sensors
                
                tic
                parfor i=1:1:size(Pp,1)
                    Mi_paths(i)=UAVSeqEstimate(obj,Pp(i,:), Tvec,  movtarg,sensids);
                end
                toc
%                 keyboard
                [~,ind]=max(Mi_paths);

                p=Pp(ind,:);
                for j=1:1:length(p)
                    obj.SensProps.SensorPos_allk{movtarg}(Tvec(j),:)=obj.SimProps.XY(p(j),:);
                end
                               
                
            end
            
            
            
            
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function MI=UAVSeqEstimate(obj,Pp, Tvec,movtarg,sensids)
                MMk=obj.TargProps.MeanState_allk;
                PPk=obj.TargProps.CovState_allk;
                
                % Load the trajectory of the UAV
                for j=1:1:length(Pp)
                    obj.SensProps.SensorPos_allk{movtarg}(Tvec(j),:)=obj.SimProps.XY(Pp(j),:);
                end
                
                for k=1:1:length(Tvec)
                    [MMk,PPk]=MeasUpdate_FDP(obj,   MMk,PPk,   1:obj.TargProps.NumbTargets   ,   sensids   ,Tvec(k));
                    
                    % -1 as it propagates from k to k+1 in dynprop function
                    if k<length(Tvec)
                         [MMk,PPk]=DynProp_FDP(obj,  MMk,PPk,    1:obj.TargProps.NumbTargets,     Tvec(k));
                    end
                end
                
                MI=MI_Tree(obj,   obj.TargProps.CovState_allk   ,   PPk ,  Tvec(1),   Tvec(end),   1:obj.TargProps.NumbTargets   );

                      
            
        end
        %============================================================================
        function obj=ForwardDPUavPathPlanner(obj,Tk)
            % Implement forward DP
            %Stattic sensors aare automatically used by default and no
            %planning is done for them. Regardless of the config of the staic sensor,
            %they can see all the targets within their FOV.
            
            Tvec=Tk:1:Tk+obj.SimProps.NoPlanningSteps ;
            
            % GOING TO FORCE TASK ALL SENSORS TO SEE EVERYTHING
            obj=BruteForceTasking(obj,Tvec);
            
             % First Copy the sensor properties from Tk all the way upto Tk1
            for k=Tvec(1:end-1)
                obj=CopySensorProps2NextTime(obj,k,k+1,1:obj.SensProps.NumbSensors);
            end
            
            % get the moving and stationary sensor ids
            [SensIDstat,SensIDmove]=GetMovingSensorIDs(obj,Tk);
            
            
            % Now doing it for each moving sensor in order
            sensids=SensIDstat;
            
            
            for movtarg=SensIDmove
                
                xsenspos0=obj.SensProps.SensorPos_allk{movtarg}(Tvec(1),:);
                [~,UAVnode0ind]=min(sqrt(sum((repmat(xsenspos0,size(obj.SimProps.XY,1),1)-obj.SimProps.XY).^2,2)));
                UAVnode0=obj.SimProps.XYindex(UAVnode0ind);

                StateTree=GenrateTreefromGrid(obj.SimProps.XY,obj.SimProps.XYADJindex,   obj.SimProps.NoPlanningSteps+1  ,xsenspos0);
                

                
                % StateTree{Tk}{i,1} is node i
                % StateTree{Tk}{i,2} is node i, left parent nodes
                % StateTree{Tk}{i,3} is node i, right child nodes
                
                % Save the optimum gain
                % StateTree{Tk}{i,4} is node i, Optimum(gain) mean
                % StateTree{Tk}{i,5} is node i, Optimum(gain) Cov
                % StateTree{Tk}{i,6} is node i, Optimum(gain) Info
                % StateTree{Tk}{i,7} is node i, previous parent optimum node: to
                % follwo back and get the whole trajectory
                
                % info is calculated from the P0 i.e. at start
                % at the end nodes, find the max info na dwork your way
                % back to the starting node to get the
                
                sensids=union( sensids , movtarg); % Sequentially add the moving sensors
                
%                 keyboard
                for k=1:1:length(Tvec)
                    
                    for i=1:1:size(StateTree{k},1)
                        if k==1
                            StateTree{k}{i,4}=obj.TargProps.MeanState_allk;
                            StateTree{k}{i,5}=obj.TargProps.CovState_allk;
                            StateTree{k}{i,6}=0;
                            StateTree{k}{i,7}=[];
                        else
                            currnode=StateTree{k}{i,1};
                            MIcurr=-1e5;
                            for parentnode=StateTree{k}{i,2}
                                
                                for s=1:1:size(StateTree{k-1},1)
                                    if StateTree{k-1}{s,1}==parentnode
                                        break
                                    end
                                end
                                parent_mu=StateTree{k-1}{s,4};
                                parent_cov=StateTree{k-1}{s,5};
                                
                                
                                % Do propagation from parent node at Tvec()
                                [curr_mu,curr_Pk]=DynProp_FDP(    obj     ,parent_mu,      parent_cov,     1:obj.TargProps.NumbTargets,Tvec(k-1)       );  % from k-1 to k
                                
                                % force the UAV sensor to current node and do measurement update
                                obj.SensProps.SensorPos_allk{movtarg}(Tvec(k),:) =obj.SimProps.XY(currnode,:);
                                [curr_mu,curr_Pk]=MeasUpdate_FDP(obj,curr_mu,curr_Pk,   1:obj.TargProps.NumbTargets,  sensids,   Tvec(k));
                                
                                MI=MI_Tree(obj,obj.TargProps.CovState_allk , curr_Pk,  Tvec(1), Tvec(k) ,1:obj.TargProps.NumbTargets);
                                
                                if MI>MIcurr
                                    StateTree{k}{i,4}=curr_mu;
                                    StateTree{k}{i,5}=curr_Pk;
                                    StateTree{k}{i,6}=MI;
                                    StateTree{k}{i,7}=parentnode;
                                    MIcurr=MI;
                                end
                            end
                            
                        end
                    end
                end
                
                % Now get the Trajectory for UAV
                Trajvec=zeros(1,length(Tvec));
                for k=length(Tvec):-1:2
                    % for the final step find the maximum MI node
                    if k==length(Tvec)
                        m=-1e5;
                        pp=0;
                        for i=1:1:size(StateTree{k},1)
                            if StateTree{k}{i,6}>m
                                pp=i;
                                m=StateTree{k}{i,6};
                            end
                        end
                        Trajvec(k)=StateTree{k}{pp,1};
                        Trajvec(k-1)=StateTree{k}{pp,7};
                        continue
                    end
                    
                    for i=1:1:size(StateTree{k},1)
                            if StateTree{k}{i,1}==Trajvec(k)
                                Trajvec(k-1)=StateTree{k}{i,7};
                            end
                    end
                    
                    
                    
                end
                
                
                
                % Debug plotting
%                 keyboard     
%                 obj.plot_Dynamic_scenario(Tvec(1));     
%                 hold on      
%                 for ij=1:1:length(Trajvec)
%                     plot(obj.SimProps.XY(Trajvec(ij),1),obj.SimProps.XY(Trajvec(ij),2),'ro' )
%                 end
%                 plot3(obj.SimProps.XY(:,1),obj.SimProps.XY(:,2),obj.SimProps.XYindex,'.')



                for k=1:1:length(Tvec)
                    try
                    obj.SensProps.SensorPos_allk{movtarg}(Tvec(k),:)=obj.SimProps.XY(Trajvec(k),:);
                    catch
                        keyboard
                    end
                end
                
            end
            
        
           
        end
        
        %============================================================================
        function obj=JointCovStaticPlanner(obj,Tk,Tk1,meth)
            Tvec=Tk:1:Tk1;
            
             % First Copy the sensor properties from Tk all the way upto Tk1
            for k=Tvec(1:end-1)
                obj=CopySensorProps2NextTime(obj,k,k+1,1:obj.SensProps.NumbSensors);
            end
            
            Nt=obj.TargProps.NumbTargets;
            Ns=obj.SensProps.NumbSensors;
            Nk=length(Tvec);
            getvarind=@(j,i,k)(k-1)*Nt*Ns+(i-1)*Ns+j;
            mpol('v',Nk*Nt*Ns,1)   % Ns is repeated Nt times, then the whole thing is repated Nk times
            [COVsStateSens_Px,COVsStateSens_Py,COVsStateSens_Pxy,LocinMatJC_Pxy,LocinMatJC_Py,DeletedSensors]=GenerateStateCov_allk(obj,Tk,Tk1);
            MPx=[];
            MPy=[];
            MPxy=[];
            delrowscolsPx=cell(obj.SensProps.NumbSensors,2);  % first row and then col
            delrowscolsPy=cell(obj.SensProps.NumbSensors,2);
            delrowscolsPxy=cell(obj.SensProps.NumbSensors,2);
            
          
            
%             detJ_Px=1;
%             trJ_Px=0;
            
            detJ_Py=1;
            
            detJ_Ppy=1;
            
            for targid=1:obj.TargProps.NumbTargets
                
                delrowscolsPx{targid,2}=find(sum(abs(COVsStateSens_Px{targid}),1)==0);
                delrowscolsPx{targid,1}=find(sum(abs(COVsStateSens_Px{targid}),2)==0);
                
                delrowscolsPy{targid,2}=find(sum(abs(COVsStateSens_Py{targid}),1)==0);
                delrowscolsPy{targid,1}=find(sum(abs(COVsStateSens_Py{targid}),2)==0);
                
                delrowscolsPxy{targid,2}=find(sum(abs(COVsStateSens_Pxy{targid}),1)==0);
                delrowscolsPxy{targid,1}=find(sum(abs(COVsStateSens_Pxy{targid}),2)==0);
                
                if isequal(sort(delrowscolsPxy{targid,2}),sort( delrowscolsPy{targid,2}))==0
                    error('u have wrong cols deleted as Pxy and Py should have same cols deleted')
                end

                % get the coreners of the cross  sub bloscks
                for targk=1:1:length(Tvec)
                    for nsens=1:1:obj.SensProps.NumbSensors
                        for sensk=1:1:length(Tvec)
                            vijk=v( getvarind(nsens,targid,sensk) );
                            ff=LocinMatJC_Pxy{targid,targk,nsens,sensk};%=[ Stveclocs{targid,targk}(1)  ,  Mtveclocs{targid,sensk,nsens}(1)  ,  Stveclocs{targid,targk}(2)  , Mtveclocs{targid,sensk,nsens}(2) ];
                            rt=ff(1);
                            ct=ff(2);
                            rb=ff(3);
                            cb=ff(4);
                            e1=ones(size(COVsStateSens_Pxy{targid}));
                            e0=zeros(size(COVsStateSens_Pxy{targid}));
                            e0(rt:rb,ct:cb)=1;
                            COVsStateSens_Pxy{targid}=COVsStateSens_Pxy{targid}.*(e0*vijk+e1);
%                             COVsStateSens_Pxy{targid}(rt:rb,ct:cb)=COVsStateSens_Pxy(rt:rb,ct:cb)*vijk;
                        end
                    end
                end
                % get the coreners of the  measurement sub bloscks
                for sensk1=1:1:length(Tvec)
                    for nsens1=1:1:obj.SensProps.NumbSensors
                         vind_1=getvarind(nsens1,targid,sensk1) ;
                         vijk_1=v(vind_1 );
                        for sensk2=1:1:length(Tvec)
                            for nsens2=1:1:obj.SensProps.NumbSensors
                                
                                vind_2=getvarind(nsens2,targid,sensk2);
                                vijk_2=v( vind_2 );
                                
                                 if vind_1==vind_2
                                    continue
                                 end
                                
                                ff=LocinMatJC_Py{targid,sensk1,nsens1,sensk2,nsens2};%=[ Mtveclocs{targid,sensk1,nsens1}(1)  ,  Mtveclocs{targid,sensk2,nsens2}(1)  ,  Mtveclocs{targid,sensk1,nsens1}(2)  , Mtveclocs{targid,sensk2,nsens2}(2) ];
                                rt=ff(1);
                                ct=ff(2);
                                rb=ff(3);
                                cb=ff(4);
                                e1=ones(size(COVsStateSens_Py{targid}));
                                e0=zeros(size(COVsStateSens_Py{targid}));
                                e0(rt:rb,ct:cb)=1;
                                COVsStateSens_Py{targid}=COVsStateSens_Py{targid}.*(e0*vijk_1*vijk_2+e1);
%                                 COVsStateSens_Py{targid}(rt:rb,ct:cb)=COVsStateSens_Py(rt:rb,ct:cb)*vijk_1*vijk_2;
                            end
                        end
                    end
                end
                %removing the 0 cols and rows in the mats
                COVsStateSens_Px{targid}(delrowscolsPx{targid,1}, :)=[];
                COVsStateSens_Px{targid}( :,delrowscolsPx{targid,2})=[];
                
                COVsStateSens_Py{targid}(delrowscolsPy{targid,1}, :)=[];
                COVsStateSens_Py{targid}( :,delrowscolsPy{targid,2})=[];
                
                COVsStateSens_Pxy{targid}(delrowscolsPxy{targid,1}, :)=[];
                COVsStateSens_Pxy{targid}( :,delrowscolsPxy{targid,2})=[];
                
                %Joining all the covs of the targets
                MPx=blkdiag(  MPx,    COVsStateSens_Px{targid});
                MPy=blkdiag(  MPy,    COVsStateSens_Py{targid});
                MPxy=blkdiag( MPxy,   COVsStateSens_Pxy{targid});
                
                detJ_Py=detJ_Py*det_leibniz( COVsStateSens_Py{targid} ,length(COVsStateSens_Py{targid}) );
                
                detJ_Ppy=detJ_Ppy*det_leibniz( COVsStateSens_Py{targid}  -  COVsStateSens_Pxy{targid}'*inv( COVsStateSens_Px{targid} )*COVsStateSens_Pxy{targid} ,length(COVsStateSens_Py{targid})  );
               
            
            end
            
             
           
            

            %removing the variables that are useless
            delconsts=[];
            for kji=1:1:length(DeletedSensors) %[k,j,i]
                
                dd=DeletedSensors(kji,:);
                dk=dd(1);
                dj=dd(2);
                di=dd(3);
                delconsts=[delconsts;v( getvarind(dj,di,dk) )];
                
            end
            
            
            if length(delconsts)>0
                delconsts=[delconsts==0];
            else
                delconsts=[];
            end
            
            % Now add the constrainsts that some sensor  can see only one
            % object
            visconst=[];
            for k=1:1:length(Tvec)
                for nsens=1:1:obj.SensProps.NumbSensors
                    if strcmpi(obj.SensProps.SensorFOVtype_allk{nsens}{Tvec(k)},'1Sensor->1Target')  % it is a track sensor, plot cones
                        vv=0;
                        for targid=1:obj.TargProps.NumbTargets
                            vv=vv+v( getvarind(nsens,targid,k) );
                        end
                        visconst=[visconst;vv==1];
                    end
                end
            end
            
            % Now adding the constraints that they are binary
            binconts=[];
            for ij=1:1:Nk*Nt*Ns
                binconts=[binconts;v(ij)-v(ij)^2==0];
            end
            
             keyboard
             
             
            if strcmpi(meth,'fullmi')
                
                pnum=det_leibniz(MPy,length(MPy) );
                pdem=det_leibniz(MPy-MPxy'*inv(MPx)*MPxy,length(MPy));
                
%                 pnum=detJ_Py;
%                 pdem=detJ_Ppy;
                
                
                conts=[binconts;delconsts;visconst];
                prob = msdp(max(pnum),mom(pdem)==1,conts);
                [status,GPolyobj] = msol(prob) ;
                vsolz = round(double(v));
                
                
                
            elseif strcmpi(meth,'mitraceratio')
                
                pnum=trace(MPy);
                pdem=trace(MPy-MPxy'*inv(MPx)*MPxy);
                
                conts=[binconts;delconsts;visconst];
                prob = msdp(max(pnum),mom(pdem)==1,conts);
                [status,GPolyobj] = msol(prob) ;
                vsolz = round(double(v));
                
            elseif strcmpi(meth,'mitracediff')
                
                pdiff=trace(MPxy'*inv(MPx)*MPxy);
                conts=[binconts;delconsts;visconst];
                prob = msdp(max(pdiff),conts);
                [status,GPolyobj] = msol(prob) ;
                vsolz = round(double(v));
                
            end
            % Now loading the Sensor tasking variable at the respective
            % times
            for k=1:1:length(Tvec)
                obj.SensProps.SensorTaskTargets_allk{Tvec(k)}=[];
                for nsens=1:obj.SensProps.NumbSensors
                    for targid=1:obj.TargProps.NumbTargets
                        if vsolz( getvarind(nsens,targid,k) )==1
                            obj.SensProps.SensorTaskTargets_allk{Tvec(k)}= vertcat(obj.SensProps.SensorTaskTargets_allk{Tvec(k)},[nsens,targid]);
                            obj=UpdateSensorProps_Tk(obj,Tvec(k),nsens,{'RotateToTarget',targid});
                        end
                    end
                end
            end
            
        end
        %============================================================================
        function obj=MatrtixPolyver_JointCovStaticPlanner(obj,Tk,Tk1,meth)
            Tvec=Tk:1:Tk1;
            
             % First Copy the sensor properties from Tk all the way upto Tk1
            for k=Tvec(1:end-1)
                obj=CopySensorProps2NextTime(obj,k,k+1,1:obj.SensProps.NumbSensors);
            end
            
            Nt=obj.TargProps.NumbTargets;
            Ns=obj.SensProps.NumbSensors;
            Nk=length(Tvec);
            getvarind=@(j,i,k)(k-1)*Nt*Ns+(i-1)*Ns+j;
            
            mpol('v',Nk*Nt*Ns,1)   % Ns is repeated Nt times, then the whole thing is repated Nk times
            polydim_withoutcoeff=Nk*Nt*Ns;
            PIP=eye(polydim_withoutcoeff);
            getsinglepoly=@(ind)PIP(ind,:);
            
            [COVsStateSens_Px,COVsStateSens_Py,COVsStateSens_Pxy,LocinMatJC_Pxy,LocinMatJC_Py,DeletedSensors]=GenerateStateCov_allk(obj,Tk,Tk1);
            MPx=[];
            MPy=[];
            MPxy=[];
            delrowscolsPx=cell(obj.SensProps.NumbSensors,2);  % first row and then col
            delrowscolsPy=cell(obj.SensProps.NumbSensors,2);
            delrowscolsPxy=cell(obj.SensProps.NumbSensors,2);
            
%             keyboard
            
%             detJ_Px=1;
%             trJ_Px=0;
            
            detJ_Py=[1,zeros(1,polydim_withoutcoeff)];
            detJ_Ppy=[1,zeros(1,polydim_withoutcoeff)];
   
            tic
            for targid=1:obj.TargProps.NumbTargets
                
                delrowscolsPx{targid,2}=find(sum(abs(COVsStateSens_Px{targid}),1)==0);
                delrowscolsPx{targid,1}=find(sum(abs(COVsStateSens_Px{targid}),2)==0);
                
                delrowscolsPy{targid,2}=find(sum(abs(COVsStateSens_Py{targid}),1)==0);
                delrowscolsPy{targid,1}=find(sum(abs(COVsStateSens_Py{targid}),2)==0);
                
                delrowscolsPxy{targid,2}=find(sum(abs(COVsStateSens_Pxy{targid}),1)==0);
                delrowscolsPxy{targid,1}=find(sum(abs(COVsStateSens_Pxy{targid}),2)==0);
                
                if isequal(sort(delrowscolsPxy{targid,2}),sort( delrowscolsPy{targid,2}))==0
                    error('u have wrong cols deleted as Pxy and Py should have same cols deleted')
                end
                COVsStateSens_Pxy{targid}=ConvertMat_2_MatrixOfPolys(COVsStateSens_Pxy{targid},polydim_withoutcoeff);
                % get the coreners of the cross  sub bloscks
                for targk=1:1:length(Tvec)
                    for nsens=1:1:obj.SensProps.NumbSensors
                        for sensk=1:1:length(Tvec)
                            vind= getvarind(nsens,targid,sensk) ;
                            ff=LocinMatJC_Pxy{targid,targk,nsens,sensk};%=[ Stveclocs{targid,targk}(1)  ,  Mtveclocs{targid,sensk,nsens}(1)  ,  Stveclocs{targid,targk}(2)  , Mtveclocs{targid,sensk,nsens}(2) ];
                            rt=ff(1);
                            ct=ff(2);
                            rb=ff(3);
                            cb=ff(4);
                            
                            locs=Simple_tens_prod_vec([rt:rb]',[ct:cb]');
                            COVsStateSens_Pxy{targid}=  PolyMultiply_MatrixOfPolys(   [1,getsinglepoly(vind)]  ,  COVsStateSens_Pxy{targid}  ,locs);
                            
                        end
                    end
                end
                
                % get the coreners of the  measurement sub bloscks
                COVsStateSens_Py{targid}=ConvertMat_2_MatrixOfPolys(COVsStateSens_Py{targid},polydim_withoutcoeff);
                for sensk1=1:1:length(Tvec)
                    for nsens1=1:1:obj.SensProps.NumbSensors
                         vind_1= getvarind(nsens1,targid,sensk1) ;
                        for sensk2=1:1:length(Tvec)
                            for nsens2=1:1:obj.SensProps.NumbSensors
                                vind_2=getvarind(nsens2,targid,sensk2) ;
                                
                                if vind_1==vind_2
                                    continue
                                end
                                
                                ff=LocinMatJC_Py{targid,sensk1,nsens1,sensk2,nsens2};%=[ Mtveclocs{targid,sensk1,nsens1}(1)  ,  Mtveclocs{targid,sensk2,nsens2}(1)  ,  Mtveclocs{targid,sensk1,nsens1}(2)  , Mtveclocs{targid,sensk2,nsens2}(2) ];
                                rt=ff(1);
                                ct=ff(2);
                                rb=ff(3);
                                cb=ff(4);
                               
                                locs=Simple_tens_prod_vec([rt:rb]',[ct:cb]');
                                vprod=multiply_polyND( [1,getsinglepoly(vind_1)],[1,getsinglepoly(vind_2)] );
                                COVsStateSens_Py{targid}=PolyMultiply_MatrixOfPolys(vprod,COVsStateSens_Py{targid},locs);
                               

                            end
                        end
                    end
                end
                
                %removing the 0 cols and rows in the mats
%                 COVsStateSens_Px{targid}(delrowscolsPx{targid,1}, :)=[];
%                 COVsStateSens_Px{targid}( :,delrowscolsPx{targid,2})=[];
%                 
%                 COVsStateSens_Py{targid}(delrowscolsPy{targid,1}, :)=[];
%                 COVsStateSens_Py{targid}( :,delrowscolsPy{targid,2})=[];
%                 
%                 COVsStateSens_Pxy{targid}(delrowscolsPxy{targid,1}, :)=[];
%                 COVsStateSens_Pxy{targid}( :,delrowscolsPxy{targid,2})=[];
                
                %Joining all the covs of the targets
%                 MPx=blkdiag(  MPx,    COVsStateSens_Px{targid});
%                 MPy=blkdiag(  MPy,    COVsStateSens_Py{targid});
%                 MPxy=blkdiag( MPxy,   COVsStateSens_Pxy{targid});
                tic
                detJ_Py=multiply_polyND( detJ_Py , Det_MatrixOfPolys(COVsStateSens_Py{targid} ) ) ; 
                toc
                disp('done with Py')
                
                invPx=inv( COVsStateSens_Px{targid});
                
                PxPolyMatInv=ConvertMat_2_MatrixOfPolys(invPx,polydim_withoutcoeff);
                
                M=AddSubMultiply_MatrixOfPolys(COVsStateSens_Pxy{targid}',PxPolyMatInv,'multiply');
                M=AddSubMultiply_MatrixOfPolys(M,COVsStateSens_Pxy{targid},'multiply');
                M=AddSubMultiply_MatrixOfPolys(COVsStateSens_Py{targid},M,'sub');
                
%                 detJ_Ppy=detJ_Ppy*Det_MatrixOfPolys(M);
                tic 
                detJ_Ppy=multiply_polyND( detJ_Ppy, Det_MatrixOfPolys(M) );
                toc   
                
                disp('done with Pxy')
%                 detJ_Ppy=detJ_Ppy*det_leibniz( COVsStateSens_Py{targid}  -  COVsStateSens_Pxy{targid}'*inv( COVsStateSens_Px{targid} )*COVsStateSens_Pxy{targid} ,length(COVsStateSens_Py{targid})  );
               
                disp(strcat(num2str(targid),' is done'))
            end
            
             
            
            

            %removing the variables that are useless
            delconsts=[];
            for kji=1:1:length(DeletedSensors) %[k,j,i]
                
                dd=DeletedSensors(kji,:);
                dk=dd(1);
                dj=dd(2);
                di=dd(3);
                delconsts=[delconsts;v( getvarind(dj,di,dk) )];
                
            end
            
            
            if length(delconsts)>0
                delconsts=[delconsts==0];
            else
                delconsts=[];
            end
            
            % Now add the constrainsts that some sensor  can see only one
            % object
            visconst=[];
            for k=1:1:length(Tvec)
                for nsens=1:1:obj.SensProps.NumbSensors
                    if strcmpi(obj.SensProps.SensorFOVtype_allk{nsens}{Tvec(k)},'1Sensor->1Target')  % it is a track sensor, plot cones
                        vv=0;
                        for targid=1:obj.TargProps.NumbTargets
                            vv=vv+v( getvarind(nsens,targid,k) );
                        end
                        visconst=[visconst;vv==1];
                    end
                end
            end
            
            % Now adding the constraints that they are binary
            binconts=[];
            for ij=1:1:Nk*Nt*Ns
                binconts=[binconts;v(ij)-v(ij)^2==0];
            end
            
            toc
            
            keyboard
            
            if strcmpi(meth,'fullmi')
                pnum=evaluate_polyND_bruteforce(detJ_Py,v');
                pdem=evaluate_polyND_bruteforce(detJ_Ppy,v');
%                 pnum=detJ_Py;
%                 pdem=detJ_Ppy;
                
                
                conts=[binconts;delconsts;visconst];
                prob = msdp(max(pnum),mom(pdem)==1,conts);
                [status,GPolyobj] = msol(prob) ;
                vsolz = round(double(v));
                
                
                
            elseif strcmpi(meth,'mitraceratio')
                
                pnum=evaluate_polyND_bruteforce(trJ_Py,v');
                pdem=evaluate_polyND_bruteforce(trJ_Ppy,v');
                
                
                conts=[binconts;delconsts;visconst];
                prob = msdp(min(pdem),conts);
                [status,GPolyobj] = msol(prob) ;
                vsolz = round(double(v));
                
            elseif strcmpi(meth,'mitracediff')
                pnum=evaluate_polyND_bruteforce(trJ_Py,v');
                pdem=evaluate_polyND_bruteforce(trJ_Ppy,v');
                pdiff=pnum-pdem;
                
                conts=[binconts;delconsts;visconst];
                prob = msdp(max(pdiff),conts);
                [status,GPolyobj] = msol(prob) ;
                vsolz = round(double(v));
                
            end
            % Now loading the Sensor tasking variable at the respective
            % times
            for k=1:1:length(Tvec)
                obj.SensProps.SensorTaskTargets_allk{Tvec(k)}=[];
                for nsens=1:obj.SensProps.NumbSensors
                    for targid=1:obj.TargProps.NumbTargets
                        if vsolz( getvarind(nsens,targid,k) )==1
                            obj.SensProps.SensorTaskTargets_allk{Tvec(k)}= vertcat(obj.SensProps.SensorTaskTargets_allk{Tvec(k)},[nsens,targid]);
                            obj=UpdateSensorProps_Tk(obj,Tvec(k),nsens,{'RotateToTarget',targid});
                        end
                    end
                end
            end
            
        end
        %============================================================================
        function obj=MatrtixPolyver_JointCovStaticPlanner_Trace(obj,Tk,Tk1,meth)
            Tvec=Tk:1:Tk1;
            
             % First Copy the sensor properties from Tk all the way upto Tk1
            for k=Tvec(1:end-1)
                obj=CopySensorProps2NextTime(obj,k,k+1,1:obj.SensProps.NumbSensors);
            end
            
            Nt=obj.TargProps.NumbTargets;
            Ns=obj.SensProps.NumbSensors;
            Nk=length(Tvec);
            getvarind=@(j,i,k)(k-1)*Nt*Ns+(i-1)*Ns+j;
            
            mpol('v',Nk*Nt*Ns,1)   % Ns is repeated Nt times, then the whole thing is repated Nk times
            polydim_withoutcoeff=Nk*Nt*Ns;
            PIP=eye(polydim_withoutcoeff);
            getsinglepoly=@(ind)PIP(ind,:);
            
            [COVsStateSens_Px,COVsStateSens_Py,COVsStateSens_Pxy,LocinMatJC_Pxy,LocinMatJC_Py,DeletedSensors]=GenerateStateCov_allk(obj,Tk,Tk1);
            MPx=[];
            MPy=[];
            MPxy=[];
            delrowscolsPx=cell(obj.SensProps.NumbSensors,2);  % first row and then col
            delrowscolsPy=cell(obj.SensProps.NumbSensors,2);
            delrowscolsPxy=cell(obj.SensProps.NumbSensors,2);
            
%             keyboard
            
%             detJ_Px=1;
%             trJ_Px=0;
            
            traceJ_Py=[0,zeros(1,polydim_withoutcoeff)];
            traceJ_Ppy=[0,zeros(1,polydim_withoutcoeff)];
   
            tic
            for targid=1:obj.TargProps.NumbTargets
                
                delrowscolsPx{targid,2}=find(sum(abs(COVsStateSens_Px{targid}),1)==0);
                delrowscolsPx{targid,1}=find(sum(abs(COVsStateSens_Px{targid}),2)==0);
                
                delrowscolsPy{targid,2}=find(sum(abs(COVsStateSens_Py{targid}),1)==0);
                delrowscolsPy{targid,1}=find(sum(abs(COVsStateSens_Py{targid}),2)==0);
                
                delrowscolsPxy{targid,2}=find(sum(abs(COVsStateSens_Pxy{targid}),1)==0);
                delrowscolsPxy{targid,1}=find(sum(abs(COVsStateSens_Pxy{targid}),2)==0);
                
                if isequal(sort(delrowscolsPxy{targid,2}),sort( delrowscolsPy{targid,2}))==0
                    error('u have wrong cols deleted as Pxy and Py should have same cols deleted')
                end
                COVsStateSens_Pxy{targid}=ConvertMat_2_MatrixOfPolys(COVsStateSens_Pxy{targid},polydim_withoutcoeff);
                % get the coreners of the cross  sub bloscks
                for targk=1:1:length(Tvec)
                    for nsens=1:1:obj.SensProps.NumbSensors
                        for sensk=1:1:length(Tvec)
                            vind= getvarind(nsens,targid,sensk) ;
                            ff=LocinMatJC_Pxy{targid,targk,nsens,sensk};%=[ Stveclocs{targid,targk}(1)  ,  Mtveclocs{targid,sensk,nsens}(1)  ,  Stveclocs{targid,targk}(2)  , Mtveclocs{targid,sensk,nsens}(2) ];
                            rt=ff(1);
                            ct=ff(2);
                            rb=ff(3);
                            cb=ff(4);
                            
                            locs=Simple_tens_prod_vec([rt:rb]',[ct:cb]');
                            COVsStateSens_Pxy{targid}=  PolyMultiply_MatrixOfPolys(   [1,getsinglepoly(vind)]  ,  COVsStateSens_Pxy{targid}  ,locs);
                            
                        end
                    end
                end
                
                % get the coreners of the  measurement sub bloscks
                COVsStateSens_Py{targid}=ConvertMat_2_MatrixOfPolys(COVsStateSens_Py{targid},polydim_withoutcoeff);
                for sensk1=1:1:length(Tvec)
                    for nsens1=1:1:obj.SensProps.NumbSensors
                         vind_1= getvarind(nsens1,targid,sensk1) ;
                        for sensk2=1:1:length(Tvec)
                            for nsens2=1:1:obj.SensProps.NumbSensors
                                vind_2=getvarind(nsens2,targid,sensk2) ;
                                
%                                 if vind_1==vind_2
%                                     continue
%                                 end
                                
                                ff=LocinMatJC_Py{targid,sensk1,nsens1,sensk2,nsens2};%=[ Mtveclocs{targid,sensk1,nsens1}(1)  ,  Mtveclocs{targid,sensk2,nsens2}(1)  ,  Mtveclocs{targid,sensk1,nsens1}(2)  , Mtveclocs{targid,sensk2,nsens2}(2) ];
                                rt=ff(1);
                                ct=ff(2);
                                rb=ff(3);
                                cb=ff(4);
                               
                                locs=Simple_tens_prod_vec([rt:rb]',[ct:cb]');
                                vprod=multiply_polyND( [1,getsinglepoly(vind_1)],[1,getsinglepoly(vind_2)] );
                                COVsStateSens_Py{targid}=PolyMultiply_MatrixOfPolys(vprod,COVsStateSens_Py{targid},locs);
                               

                            end
                        end
                    end
                end
                
                %removing the 0 cols and rows in the mats
%                 COVsStateSens_Px{targid}(delrowscolsPx{targid,1}, :)=[];
%                 COVsStateSens_Px{targid}( :,delrowscolsPx{targid,2})=[];
%                 
%                 COVsStateSens_Py{targid}(delrowscolsPy{targid,1}, :)=[];
%                 COVsStateSens_Py{targid}( :,delrowscolsPy{targid,2})=[];
%                 
%                 COVsStateSens_Pxy{targid}(delrowscolsPxy{targid,1}, :)=[];
%                 COVsStateSens_Pxy{targid}( :,delrowscolsPxy{targid,2})=[];
                
                %Joining all the covs of the targets
%                 MPx=blkdiag(  MPx,    COVsStateSens_Px{targid});
%                 MPy=blkdiag(  MPy,    COVsStateSens_Py{targid});
%                 MPxy=blkdiag( MPxy,   COVsStateSens_Pxy{targid});
                tic
                traceJ_Py=add_sub_polyND( traceJ_Py , Trace_MatrixOfPolys(COVsStateSens_Py{targid} ) ,'add') ; 
                toc
                disp('done with Py')
                
                invPx=inv( COVsStateSens_Px{targid});
                
                PxPolyMatInv=ConvertMat_2_MatrixOfPolys(invPx,polydim_withoutcoeff);
                
                M=AddSubMultiply_MatrixOfPolys(COVsStateSens_Pxy{targid}',PxPolyMatInv,'multiply');
                M=AddSubMultiply_MatrixOfPolys(M,COVsStateSens_Pxy{targid},'multiply');
                M=AddSubMultiply_MatrixOfPolys(COVsStateSens_Py{targid},M,'sub');
                
%                 detJ_Ppy=detJ_Ppy*Det_MatrixOfPolys(M);
                tic 
                traceJ_Ppy=add_sub_polyND( traceJ_Ppy , Trace_MatrixOfPolys(M) ,'add') ;
%                 detJ_Ppy=multiply_polyND( detJ_Ppy, Det_MatrixOfPolys(M) );
                toc   
                
                disp('done with Pxy')
%                 detJ_Ppy=detJ_Ppy*det_leibniz( COVsStateSens_Py{targid}  -  COVsStateSens_Pxy{targid}'*inv( COVsStateSens_Px{targid} )*COVsStateSens_Pxy{targid} ,length(COVsStateSens_Py{targid})  );
               
                disp(strcat(num2str(targid),' is done'))
            end
            
             
            
            

            %removing the variables that are useless
            delconsts=[];
            for kji=1:1:length(DeletedSensors) %[k,j,i]
                
                dd=DeletedSensors(kji,:);
                dk=dd(1);
                dj=dd(2);
                di=dd(3);
                delconsts=[delconsts;v( getvarind(dj,di,dk) )];
                
            end
            
            
            if length(delconsts)>0
                delconsts=[delconsts==0];
            else
                delconsts=[];
            end
            
            % Now add the constrainsts that some sensor  can see only one
            % object
            visconst=[];
            for k=1:1:length(Tvec)
                for nsens=1:1:obj.SensProps.NumbSensors
                    if strcmpi(obj.SensProps.SensorFOVtype_allk{nsens}{Tvec(k)},'1Sensor->1Target')  % it is a track sensor, plot cones
                        vv=0;
                        for targid=1:obj.TargProps.NumbTargets
                            vv=vv+v( getvarind(nsens,targid,k) );
                        end
                        visconst=[visconst;vv==1];
                    end
                end
            end
            
            % Now adding the constraints that they are binary
            binconts=[];
            for ij=1:1:Nk*Nt*Ns
                binconts=[binconts;v(ij)-v(ij)^2==0];
            end
            
            toc
            
            keyboard
            
            if strcmpi(meth,'fullmi')
                pnum=evaluate_polyND_bruteforce(detJ_Py,v');
                pdem=evaluate_polyND_bruteforce(detJ_Ppy,v');
%                 pnum=detJ_Py;
%                 pdem=detJ_Ppy;
                
                
                conts=[binconts;delconsts;visconst];
                prob = msdp(max(pnum),mom(pdem)==1,conts);
                [status,GPolyobj] = msol(prob) ;
                vsolz = round(double(v));
                
                
                
            elseif strcmpi(meth,'mitraceratio')
                
                pnum=evaluate_polyND_bruteforce(traceJ_Py,v');
                pdem=evaluate_polyND_bruteforce(traceJ_Ppy,v');
                
                
                conts=[binconts;delconsts;visconst];
                prob = msdp(max(pnum),mom(pdem)==1,conts);
                [status,GPolyobj] = msol(prob) ;
                vsolz = round(double(v));
                
            elseif strcmpi(meth,'mitracediff')
                pnum=evaluate_polyND_bruteforce(trJ_Py,v');
                pdem=evaluate_polyND_bruteforce(trJ_Ppy,v');
                pdiff=pnum-pdem;
                
                conts=[binconts;delconsts;visconst];
                prob = msdp(max(pdiff),conts);
                [status,GPolyobj] = msol(prob) ;
                vsolz = round(double(v));
                
            end
            % Now loading the Sensor tasking variable at the respective
            % times
            for k=1:1:length(Tvec)
                obj.SensProps.SensorTaskTargets_allk{Tvec(k)}=[];
                for nsens=1:obj.SensProps.NumbSensors
                    for targid=1:obj.TargProps.NumbTargets
                        if vsolz( getvarind(nsens,targid,k) )==1
                            obj.SensProps.SensorTaskTargets_allk{Tvec(k)}= vertcat(obj.SensProps.SensorTaskTargets_allk{Tvec(k)},[nsens,targid]);
                            obj=UpdateSensorProps_Tk(obj,Tvec(k),nsens,{'RotateToTarget',targid});
                        end
                    end
                end
            end
            
        end
        %============================================================================
        function [COVsStateSens_Px,COVsStateSens_Py,COVsStateSens_Pxy,LocinMatJC_Pxy,LocinMatJC_Py,DeletedSensors]=GenerateStateCov_allk(obj,Tk,Tk1)
            % v is teh decision variables
            % Quad points from Tk time stwp itself, the tasking also starts
            % form Tk
            
            Tvec=Tk:1:Tk1;
            
            % GOING TO FORCE TASK ALL SENSORS TO SEE EVERYTHING
            obj=BruteForceTasking(obj,Tvec);
            
             % First Copy the sensor properties from Tk all the way upto Tk1
            for k=Tvec(1:end-1)
                obj=CopySensorProps2NextTime(obj,k,k+1,1:obj.SensProps.NumbSensors);
            end
            
            
            COVsStateSens_Px=cell(obj.TargProps.NumbTargets,1);
            COVsStateSens_Py=cell(obj.TargProps.NumbTargets,1);
            COVsStateSens_Pxy=cell(obj.TargProps.NumbTargets,1);
            LocinMatJC_Pxy=cell(obj.TargProps.NumbTargets,length(Tvec),obj.SensProps.NumbSensors,length(Tvec));
            LocinMatJC_Py=cell(obj.TargProps.NumbTargets,length(Tvec),obj.SensProps.NumbSensors,length(Tvec),obj.SensProps.NumbSensors);
            
            
            X=cell(length(Tvec),obj.TargProps.NumbTargets);
            Z=cell(length(Tvec),obj.SensProps.NumbSensors,obj.TargProps.NumbTargets);
            DeletedSensors=[];
            W=cell(obj.TargProps.NumbTargets,1);
            
            obj=StoreSensProps(obj);
        
            for targid=1:obj.TargProps.NumbTargets
                for k=1:1:length(Tvec)
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % state sigma points over time
                    if k==1
                        Xfk=obj.TargProps.MeanState_allk{targid}(Tk,:)';
                        Pfk=reshape(obj.TargProps.CovState_allk{targid}(Tk,:),obj.TargProps.TargetModelDim_allk{targid}(Tk),obj.TargProps.TargetModelDim_allk{targid}(Tk));
                        [XTk,wTk]=obj.SimProps.SigmaPntFunction(Xfk,Pfk);
                        X{1,targid}=XTk;
                        W{targid}=wTk;
                    else
                        
                        if strcmpi(obj.TargProps.TargetDynTypes_allk{targid}{Tvec(k)},'Stationary')  || strcmpi(obj.TargProps.TargetDynTypes_allk{targid}{Tvec(k)},'VirtualStationary')
                            X{k,targid}=X{k-1,targid};
                            
                        else
                            Xk=X{k-1,targid};
                            Xk1=zeros(size(Xk));
                            for i=1:1:length(W{targid})
                                switch(obj.TargProps.TargetModel_allk{targid}{Tvec(k)}) % type of dynamics
                                    case 'UM'  % UM motion
                                        Xk1(i,:)=KIRB_UM_eg_dyn_disc(Xk(i,:)',obj.SimProps.Time_dt);
                                    case 'CT'  % CT motion
                                        Xk1(i,:)=KIRB_CT_eg_dyn_disc(Xk(i,:)',obj.SimProps.Time_dt);
                                    case 'NoDyn'  % No dynamics case
                                        Xk1(i,:)=Xk(i,:);
                                    case 'diffdyn'
%                                         keyboard
                                        m0=obj.TargProps.TargetModelParas_allk{targid}{Tvec(k),1}.m0;
                                        cf=obj.TargProps.TargetModelParas_allk{targid}{Tvec(k),1}.cf;
                                        Xk1(i,:)=diffdyn(Xk(i,:),m0,cf);
                                end
                            end
                            X{k,targid}=Xk1;
                        end
                        [xk1,Pk1]=MeanCov(X{k,targid},W{targid});
                        obj.TargProps.MeanState_allk{targid}(Tvec(k),:)=xk1;
                        obj.TargProps.CovState_allk{targid}(Tvec(k),:)=reshape(Pk1,1,obj.TargProps.TargetModelDim_allk{targid}(Tvec(k))^2);
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % measurement sigma points at this time and target
                    Xk=X{k,targid};
                    
                    for nsens=1:obj.SensProps.NumbSensors
                        
                        
                        Y=[];
                        obj=UpdateSensorProps_Tk(obj,Tvec(k),nsens,{'RotateToTarget',targid});
                        [yy,H,G]=SensorMeasurement_pts(obj,Tvec(k),Xk,nsens);
                        
                        if strcmpi(obj.SensProps.SensorOnOff_allk{nsens}{Tvec(k)},'Off')
                            
                            yy=yy*0; % if no sig pt is visible just dont use this sensor
                            DeletedSensors=vertcat(DeletedSensors,[k,nsens,targid]);
                            Y=horzcat(Y,yy);
                            
                        elseif obj.SimProps.SimMode==-1  % no FOV anywhere
                            
                            Y=horzcat(Y,yy);
                            %                             RR=blkdiag(RR,obj.SensProps.R_allk{nsens}{Tk});
                            
                        elseif obj.SimProps.SimMode==0 % During Tasking: Use FOV, During Measurement : No FOV
                            
                            if length(find(H==0))>=3/4*length(H)
                                yy=yy*0; % if no sig pt is visible just dont use this sensor
                                DeletedSensors=vertcat(DeletedSensors,[k,nsens,targid]);
                            end
                            
                            Y=horzcat(Y,yy);
                            %                             RR=blkdiag(RR,obj.SensProps.R_allk{nsens}{Tk});
                            
                        elseif obj.SimProps.SimMode==1 % During Tasking: Use FOV, During Measurement : Use FOV, but real measurements can be outside of FOV
                            
                            if length(find(H==0))>=3/4*length(H)
                                yy=yy*0; % if no sig pt is visible just dont use this sensor
                                DeletedSensors=vertcat(DeletedSensors,[k,nsens,targid]);
                            end
                            
                            
                            Y=horzcat(Y,yy);
                            %                             RR=blkdiag(RR,obj.SensProps.R_allk{nsens}{Tk});
                            
                            
                        elseif obj.SimProps.SimMode==2  % During Tasking: Use FOV, During Measurement : Use FOV, real measurements have to be inside of FOV
                            
                            if length(find(H==0))>=3/4*length(H)
                                yy=yy*0; % if no sig pt is visible just dont use this sensor
                                DeletedSensors=vertcat(DeletedSensors,[k,nsens,targid]);
                            end
                            
                            
                            Y=horzcat(Y,yy);
                            %                             RR=blkdiag(RR,obj.SensProps.R_allk{nsens}{Tk});
                            
                        else
                            error('Simulation mode unknown')
                        end

%                             
%                             clf
%                             disp(strcat('time = ',num2str(k)))
%                             disp(strcat('nsens = ',num2str(nsens)))
%                             disp(strcat('target = ',num2str(targid)))
%                             obj.plot_Dynamic_scenario(Tvec(k))
%                             Y
%                             keyboard
                            

                        
                        Z{k,nsens,targid}=Y;
                        % Restore the rotation to the target
                        obj=UpdateSensorProps_Tk(obj,Tvec(k),nsens,{'RestoreBackupDirn',NaN});
                        
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                             
                end
            end
%             keyboard
%             LocinMatJC_Pxy=cell(obj.TargProps.NumbTargets,length(Tvec),obj.SensProps.NumbSensors,length(Tvec));
%             LocinMatJC_Py=cell(obj.TargProps.NumbTargets,length(Tvec),obj.SensProps.NumbSensors,length(Tvec),obj.SensProps.NumbSensors);
            Stveclocs=cell(obj.TargProps.NumbTargets,length(Tvec));
            Mtveclocs=cell(obj.TargProps.NumbTargets,length(Tvec),obj.SensProps.NumbSensors);
            % Now computing the cross correlations between the different combos
            % First state over time
            
            for targid=1:obj.TargProps.NumbTargets
                MasterQ=[];
                MasterSigpts=[];
                ns=0;
                pp=0;
                for k=1:1:length(Tvec)
                    MasterSigpts=horzcat(MasterSigpts,X{k,targid});
                    ns=ns+size(X{k,targid},2);
                    
                    Stveclocs{targid,k}=[pp+1,pp+size(X{k,targid},2) ]; 
                    pp=pp+size(X{k,targid},2);
                    
                    MasterQ=blkdiag(MasterQ,obj.TargProps.Q_allk{targid}{Tvec(k)});
                end
                MasterR=[];
                for k=1:1:length(Tvec)
                    for nsens=1:1:obj.SensProps.NumbSensors
                        MasterR=blkdiag(MasterR,obj.SensProps.R_allk{nsens}{Tvec(k)});
                        MasterSigpts=horzcat(MasterSigpts,Z{k,nsens,targid});
                        
                        Mtveclocs{targid,k,nsens}=[pp+1,pp+size(X{k,targid},2)];
                        pp=pp+size(Z{k,nsens,targid},2);
                    end
                end
                
                [~,JC]=MeanCov(MasterSigpts, W{targid});
                JC(1:ns,1:ns)=JC(1:ns,1:ns)+MasterQ;
                JC(ns+1:end,ns+1:end)=JC(ns+1:end,ns+1:end)+MasterR;
                
                COVsStateSens_Px{targid}=JC(1:ns,1:ns);
                COVsStateSens_Py{targid}=JC(ns+1:end,ns+1:end);
                COVsStateSens_Pxy{targid}=JC(1:ns,ns+1:end);
                
%               % get the coreners of the cross  sub bloscks  
                for targk=1:1:length(Tvec)
                    for nsens=1:1:obj.SensProps.NumbSensors
                        for sensk=1:1:length(Tvec)
                            LocinMatJC_Pxy{targid,targk,nsens,sensk}=[ Stveclocs{targid,targk}(1)  ,  Mtveclocs{targid,sensk,nsens}(1)-ns ,  Stveclocs{targid,targk}(2)  , Mtveclocs{targid,sensk,nsens}(2)-ns ];
                        end
                    end
                end
                % get the coreners of the  measurement sub bloscks  
                for sensk1=1:1:length(Tvec)
                    for nsens1=1:1:obj.SensProps.NumbSensors
                        for sensk2=1:1:length(Tvec)
                            for nsens2=1:1:obj.SensProps.NumbSensors
                            LocinMatJC_Py{targid,sensk1,nsens1,sensk2,nsens2}=[ Mtveclocs{targid,sensk1,nsens1}(1)-ns,  Mtveclocs{targid,sensk2,nsens2}(1)-ns ,  Mtveclocs{targid,sensk1,nsens1}(2)-ns  , Mtveclocs{targid,sensk2,nsens2}(2)-ns ];
                            end
                        end
                    end
                end
                
            end
            
            
            
            
        end
          %============================================================================
%         function MI=MI_Current_Backup(obj,Tk,Tk1,targs)
%             % calculate the MI between the backed up config at Tk
%             % and the current target states at TK1
%             % backed up is like prior and current is like posterior
%             % caculate the MI for the targets and add them up ( assuming the targets are independent)
%             MI=0;
%             for targid=targs(:)'
%                 
%                 Pk=reshape(obj.BackUp_targCovss{targid}(Tk,:),obj.TargProps.TargetModelDim_allk{targid}(k),obj.TargProps.TargetModelDim_allk{targid}(k));
%                 Pk1=reshape(obj.TargProps.CovState_allk{targid}(Tk1,:),obj.TargProps.TargetModelDim_allk{targid}(k),obj.TargProps.TargetModelDim_allk{targid}(k));
%                 if det(Pk)/det(Pk1)<1.001  %Pf/Pu
%                     MI=MI+0;
%                 else
%                     MI=MI+0.5*log(det(Pk)/det(Pk1));
%                 end
%             end
%         end
        %============================================================================
        function MI=MI_Tree(obj,Pk,Pk1,Tk,Tk1,targs)
            % calculate the MI between the backed up config at Tk
            % and the current target states at TK1
            % backed up is like prior and current is like posterior
            % caculate the MI for the targets and add them up ( assuming the targets are independent)
            MI=0;
            for targid=targs(:)'
                
                PPk=reshape(Pk{targid}(Tk,:),obj.TargProps.TargetModelDim_allk{targid}(Tk),obj.TargProps.TargetModelDim_allk{targid}(Tk));
                PPk1=reshape(Pk1{targid}(Tk1,:),obj.TargProps.TargetModelDim_allk{targid}(Tk1),obj.TargProps.TargetModelDim_allk{targid}(Tk1));
                if det(PPk)/det(PPk1)<1.01  %Pf/Pu
                    MI=MI+0;
                else
                    MI=MI+0.5*log(det(PPk)/det(PPk1));
                end
            end
        end
        %============================================================================
        function plot_scenario(obj) %No Tasking is considered  obj.TargProps.NumbTargets*obj.SensProps.NumbSensors
            figure(1)
            clf
            hold on
            Tk=obj.SimProps.Time.tk;
            % plotting the sensors and their FOV configurations
            for sensid=1:1:obj.SensProps.NumbSensors    % all sensors are blue
                pos=obj.SensProps.SensorPos_allk{sensid}(Tk,:);
                rmax=obj.SensProps.SensorRmax_allk{sensid}(Tk,:);
                alp=obj.SensProps.SensorAlpha_allk{sensid}(Tk,:);
                dr=obj.SensProps.SensorDir_allk{sensid}(Tk,:);
                if strcmpi(obj.SensProps.SensorFOVtype_allk{sensid}{Tk},'1Sensor->AllTarget')   % it is a coverage sensor, plot circles
                    plot_ellipse(pos(1),pos(2),rmax,rmax,0,'b','fill');
                end
                if strcmpi(obj.SensProps.SensorFOVtype_allk{sensid}{Tk},'1Sensor->1Target')  % it is a track sensor, plot cones
                    plot_2Dcone(pos(1),pos(2),alp,rmax,dr,'g','fill')
                end
                
                %                 if strcmpi(obj.SensProps.SensorDynMotionModel_allk{sensid}{Tk},'UAVgridOnly') || strcmpi(obj.SensProps.SensorDynMotionModel_allk{sensid}{Tk},'UAVgrid+smoothTraj')
                plot(pos(1),pos(2),'bo','linewidth',2,'MarkerSize',6);
                %                 end
                if strcmpi(obj.SensProps.SensorDynMotionModel_allk{sensid}{Tk},'UAVgridOnly') || strcmpi(obj.SensProps.SensorDynMotionModel_allk{sensid}{Tk},'UAVgrid+smoothTraj')
                    UAVcolors=['r','b','k','m','c'];
                    poses=obj.SensProps.SensorPos_allk{sensid}(1:Tk,:);
                    cc=obj.SensProps.SensorsID(sensid);
                    plot(poses(:,1),poses(:,2),UAVcolors(cc),'linewidth',2,'MarkerSize',6);
                end
                
            end
            %plot the target with a triangle and an ellipse for confidence
            for targid=1:1:obj.TargProps.NumbTargets
                if strcmpi(obj.TargProps.VisibilityState_allk{targid}{Tk},'NotTracked')
                    continue;
                end
                %plot covariance ellipse
                
                if strcmpi(obj.TargProps.VisibilityState_allk{targid}{Tk},'Invisible') || strcmpi(obj.TargProps.VisibilityState_allk{targid}{Tk},'Lost')
                    
                else
                    P=reshape(obj.TargProps.CovState_allk{targid}(Tk,:),obj.TargProps.TargetModelDim_allk{targid}(Tk),obj.TargProps.TargetModelDim_allk{targid}(Tk));
                    plot_1sigellip(obj.TargProps.MeanState_allk{targid}(Tk,1:2),P(1:2,1:2),'r');
                end
                
                %plot triangle as mean
                plot_triangle(obj.TargProps.MeanState_allk{targid}(Tk,1),obj.TargProps.MeanState_allk{targid}(Tk,2),2,1,0,'r','fill')
                %plot the trajectories est+truth
                plot(obj.TargProps.MeanState_allk{targid}(Tk,1),obj.TargProps.MeanState_allk{targid}(Tk,2),'r','linewidth',2)
                plot(obj.TargProps.TruePosState_allk{targid}(Tk,1),obj.TargProps.TruePosState_allk{targid}(Tk,2),'k--','linewidth',2)
            end
            
            
            hold off
        end
        %============================================================================
        function plot_Dynamic_scenario(obj,Tk) %Tasking is considered
            figure(1)
            clf
            axis([obj.SimProps.Xboundary(1)-0.2*obj.SimProps.Xboundary(2),1.5*obj.SimProps.Xboundary(2),obj.SimProps.Yboundary(1)-0.2*obj.SimProps.Yboundary(2),1.5*obj.SimProps.Yboundary(2)])
            hold on
%             keyboard
%             Tk=obj.SimProps.Time.tk;
            % plotting the sensors and their FOV configurations
            for sensid=1:1:obj.SensProps.NumbSensors   % all sensors are blue
                pos=obj.SensProps.SensorPos_allk{sensid}(Tk,:);
                rmax=obj.SensProps.SensorRmax_allk{sensid}(Tk,:);
                alp=obj.SensProps.SensorAlpha_allk{sensid}(Tk,:);
                dr=obj.SensProps.SensorDir_allk{sensid}(Tk,:);
                
                if isempty(find(obj.SensProps.SensorTaskTargets_allk{Tk}(:,1)==sensid))
                    continue;
                end
                
                plot(pos(1),pos(2),'ko','linewidth',2,'MarkerSize',6)
                if strcmpi(obj.SensProps.SensorFOVtype_allk{sensid}{Tk},'1Sensor->AllTarget')   % it is a coverage sensor, plot circles
                    plot_ellipse(pos(1),pos(2),rmax,rmax,0,'b','fill');
                end
                if strcmpi(obj.SensProps.SensorFOVtype_allk{sensid}{Tk},'1Sensor->1Target')  % it is a track sensor, plot cones
                    plot_2Dcone(pos(1),pos(2),alp,rmax,dr,'g','fill')
                end
                if strcmpi(obj.SensProps.SensorDynMotionModel_allk{sensid}{Tk},'UAVgridOnly') || strcmpi(obj.SensProps.SensorDynMotionModel_allk{sensid}{Tk},'UAVgrid+smoothTraj')
                    UAVcolors=['r','b','k','m','c'];
                    poses=obj.SensProps.SensorPos_allk{sensid}(1:Tk,:);
                    cc=obj.SensProps.SensorsID(sensid);
                    plot(poses(:,1),poses(:,2),UAVcolors(cc),'linewidth',2,'MarkerSize',6);
                end
                
            end
            %plot the target with a triangle and an ellipse for confidence
            for targid=1:1:obj.TargProps.NumbTargets
                if strcmpi(obj.TargProps.VisibilityState_allk{targid}{Tk},'NotTracked')
                    continue;
                end
                %plot covariance ellipse
                
                if strcmpi(obj.TargProps.VisibilityState_allk{targid}{Tk},'Invisible') || strcmpi(obj.TargProps.VisibilityState_allk{targid}{Tk},'Lost')
                    
                else
                    try
                    P=reshape(obj.TargProps.CovState_allk{targid}(Tk,:),obj.TargProps.TargetModelDim_allk{targid}(Tk),obj.TargProps.TargetModelDim_allk{targid}(Tk));
                    catch
                        keyboard
                    end
                    plot_1sigellip(obj.TargProps.MeanState_allk{targid}(Tk,1:2),P(1:2,1:2),'r');
                end
                
                %plot triangle as mean
                plot_triangle(obj.TargProps.MeanState_allk{targid}(Tk,1),obj.TargProps.MeanState_allk{targid}(Tk,2),2,1,0,'r','fill')
                %plot the trajectories est+truth
                plot(obj.TargProps.MeanState_allk{targid}(Tk,1),obj.TargProps.MeanState_allk{targid}(Tk,2),'r','linewidth',2)
                plot(obj.TargProps.TruePosState_allk{targid}(1:Tk,1),obj.TargProps.TruePosState_allk{targid}(1:Tk,2),'k*-','linewidth',2)
            end

            axis([obj.SimProps.Xboundary,obj.SimProps.Yboundary])
            obj.plotting_properties()
            hold off
        end
        %============================================================================
        function plot_truth_scenario(obj) %No Tasking is considered+ ploy all truth
            figure('units','normalized','outerposition',[0 0 1 1])
            axis([obj.SimProps.Xboundary(1)-0.2*obj.SimProps.Xboundary(2),1.5*obj.SimProps.Xboundary(2),obj.SimProps.Yboundary(1)-0.2*obj.SimProps.Yboundary(2),1.5*obj.SimProps.Yboundary(2)])
            clf
            hold on
            Tk=obj.SimProps.Time.tk;
            % plotting the sensors and their FOV configurations
            for sensid=1:1:obj.SensProps.NumbSensors    % all sensors are blue
                pos=obj.SensProps.SensorPos_allk{sensid}(Tk,:);
                rmax=obj.SensProps.SensorRmax_allk{sensid}(Tk,:);
                alp=obj.SensProps.SensorAlpha_allk{sensid}(Tk,:);
                dr=obj.SensProps.SensorDir_allk{sensid}(Tk,:);
                %                 keyboard
                if strcmpi(obj.SensProps.SensorFOVtype_allk{sensid}{Tk},'1Sensor->AllTarget')   % it is a coverage sensor, plot circles
                    plot_ellipse(pos(1),pos(2),rmax,rmax,0,'b','fill');
                end
                if strcmpi(obj.SensProps.SensorFOVtype_allk{sensid}{Tk},'1Sensor->1Target')  % it is a track sensor, plot cones
                    plot_2Dcone(pos(1),pos(2),alp,rmax,dr,'g','fill')
                end
                plot(pos(1),pos(2),'ko','linewidth',2,'MarkerSize',6)
                
                if strcmpi(obj.SensProps.SensorDynMotionModel_allk{sensid}{Tk},'UAVgridOnly') || strcmpi(obj.SensProps.SensorDynMotionModel_allk{sensid}{Tk},'UAVgrid+smoothTraj')
                    UAVcolors=['r','b','k','m','c'];
                    poses=obj.SensProps.SensorPos_allk{sensid}(1:Tk,:);
                    cc=obj.SensProps.SensorsID(sensid);
                    plot(poses(:,1),poses(:,2),UAVcolors(cc),'linewidth',2,'MarkerSize',6);
                end
            end
            
            %plot the target with a triangle and an ellipse for confidence
            for targid=1:1:obj.TargProps.NumbTargets
                if strcmpi(obj.TargProps.VisibilityState_allk{targid}{Tk},'Invisible') || strcmpi(obj.TargProps.VisibilityState_allk{targid}{Tk},'Lost')
                    
                else
                    P=reshape(obj.TargProps.CovState_allk{targid}(Tk,:),obj.TargProps.TargetModelDim_allk{targid}(Tk),obj.TargProps.TargetModelDim_allk{targid}(Tk));
                    plot_1sigellip(obj.TargProps.MeanState_allk{targid}(Tk,1:2),P(1:2,1:2),'r');
                end
                
                
                %                 plot_triangle(obj.TargProps.MeanState_allk{targid}(Tk,1),obj.TargProps.MeanState_allk{targid}(Tk,2),1,0.5,0,'r','fill')
                %plot the trajectories est+truth
                
                %plot triangle as mean
                if strcmpi(obj.TargProps.TargetDynTypes_allk{targid}{Tk},'VirtualMove')  || strcmpi(obj.TargProps.TargetDynTypes_allk{targid}{Tk},'VirtualStationary') % if it is a vitual target
                    plot_triangle(obj.TargProps.MeanState_allk{targid}(Tk,1),obj.TargProps.MeanState_allk{targid}(Tk,2),1,1,0,'k','fill')
                    
                else
                    plot_triangle(obj.TargProps.MeanState_allk{targid}(Tk,1),obj.TargProps.MeanState_allk{targid}(Tk,2),1,1/3,0,'r','fill')
                    
                end
                plot(obj.TargProps.MeanState_allk{targid}(Tk,1),obj.TargProps.MeanState_allk{targid}(Tk,2),'r','linewidth',2)
                plot(obj.TargProps.TruePosState_allk{targid}(:,1),obj.TargProps.TruePosState_allk{targid}(:,2),'k*-','linewidth',2)
                
            end
            
            
            for sensid=1:1:obj.SensProps.NumbSensors
                pos=obj.SensProps.SensorPos_allk{sensid}(Tk,:);
                
                
                [Xf, Yf] = ds2nfu(pos(1), pos(2));
                W=0.1;
                H=0.05;
                XY=[Xf,Yf]+0.03*randn(1,2);
                annotation('textbox',[XY,W,H],'String',strcat('Sensor : ' ,num2str(sensid)),'FontSize',12,'EdgeColor','none','Color','r');
                
            end
            for targid=1:1:obj.TargProps.NumbTargets
                pos=obj.TargProps.MeanState_allk{targid}(Tk,1:2);
                [Xf, Yf] = ds2nfu(pos(1), pos(2));
                W=0.05;
                H=0.02;
                XY=[Xf,Yf]+0.001*randn(1,2);
                annotation('textbox',[XY,W,H],'String',strcat('Target : ' ,num2str(targid)),'FontSize',7,'EdgeColor','none');
                
            end
            %
            %             axis square
            %             axis equal
            axis([obj.SimProps.Xboundary,obj.SimProps.Yboundary])
            obj.plotting_properties()
            hold off
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        function plot_tasking_debug_scenario(obj) %No Tasking is considered+ ploy all truth
            figure('units','normalized','outerposition',[0 0 1 1])
            axis([obj.SimProps.Xboundary(1)-0.2*obj.SimProps.Xboundary(2),1.5*obj.SimProps.Xboundary(2),obj.SimProps.Yboundary(1)-0.2*obj.SimProps.Yboundary(2),1.5*obj.SimProps.Yboundary(2)])
            clf
            hold on
            Tk=obj.SimProps.Time.tk;
            % plotting the sensors and their FOV configurations
            for sensid=1:1:obj.SensProps.NumbSensors    % all sensors are blue
                pos=obj.SensProps.SensorPos_allk{sensid}(Tk,:);
                rmax=obj.SensProps.SensorRmax_allk{sensid}(Tk,:);
                alp=obj.SensProps.SensorAlpha_allk{sensid}(Tk,:);
                dr=obj.SensProps.SensorDir_allk{sensid}(Tk,:);
                %                 keyboard
                if strcmpi(obj.SensProps.SensorFOVtype_allk{sensid}{Tk},'1Sensor->AllTarget')   % it is a coverage sensor, plot circles
                    plot_ellipse(pos(1),pos(2),rmax,rmax,0,'b','fill');
                end
                if strcmpi(obj.SensProps.SensorFOVtype_allk{sensid}{Tk},'1Sensor->1Target')  % it is a track sensor, plot cones
                    plot_2Dcone(pos(1),pos(2),alp,rmax,dr,'g','fill')
                end
                plot(pos(1),pos(2),'ko','linewidth',2,'MarkerSize',6)
                
                if strcmpi(obj.SensProps.SensorDynMotionModel_allk{sensid}{Tk},'UAVgridOnly') || strcmpi(obj.SensProps.SensorDynMotionModel_allk{sensid}{Tk},'UAVgrid+smoothTraj')
                    UAVcolors=['r','b','k','m','c'];
                    poses=obj.SensProps.SensorPos_allk{sensid}(1:Tk,:);
                    cc=obj.SensProps.SensorsID(sensid);
                    plot(poses(:,1),poses(:,2),UAVcolors(cc),'linewidth',2,'MarkerSize',6);
                end
            end
            
            %plot the target with a triangle and an ellipse for confidence
            for targid=1:1:obj.TargProps.NumbTargets
                if strcmpi(obj.TargProps.VisibilityState_allk{targid}{Tk},'Invisible') || strcmpi(obj.TargProps.VisibilityState_allk{targid}{Tk},'Lost')
                    
                else
                    P=reshape(obj.TargProps.CovState_allk{targid}(Tk,:),obj.TargProps.TargetModelDim_allk{targid}(Tk),obj.TargProps.TargetModelDim_allk{targid}(Tk));
                    plot_1sigellip(obj.TargProps.MeanState_allk{targid}(Tk,1:2),P(1:2,1:2),'r');
                end
                
                
                %                 plot_triangle(obj.TargProps.MeanState_allk{targid}(Tk,1),obj.TargProps.MeanState_allk{targid}(Tk,2),1,0.5,0,'r','fill')
                %plot the trajectories est+truth
                
                %plot triangle as mean
                if strcmpi(obj.TargProps.TargetDynTypes_allk{targid}{Tk},'VirtualMove')  || strcmpi(obj.TargProps.TargetDynTypes_allk{targid}{Tk},'VirtualStationary') % if it is a vitual target
                    plot_triangle(obj.TargProps.MeanState_allk{targid}(Tk,1),obj.TargProps.MeanState_allk{targid}(Tk,2),1,1,0,'k','fill')
                    
                else
                    plot_triangle(obj.TargProps.MeanState_allk{targid}(Tk,1),obj.TargProps.MeanState_allk{targid}(Tk,2),1,1/3,0,'r','fill')
                    
                end
                plot(obj.TargProps.MeanState_allk{targid}(Tk,1),obj.TargProps.MeanState_allk{targid}(Tk,2),'r','linewidth',2)
                plot(obj.TargProps.TruePosState_allk{targid}(:,1),obj.TargProps.TruePosState_allk{targid}(:,2),'k--','linewidth',2)
                
            end
            
            
            for sensid=1:1:obj.SensProps.NumbSensors
                pos=obj.SensProps.SensorPos_allk{sensid}(Tk,:);
                
                
                [Xf, Yf] = ds2nfu(pos(1), pos(2));
                W=0.1;
                H=0.05;
                XY=[Xf,Yf]+0.03*randn(1,2);
                annotation('textbox',[XY,W,H],'String',strcat('Sensor : ' ,num2str(sensid)),'FontSize',12,'EdgeColor','none','Color','r');
                
            end
            for targid=1:1:obj.TargProps.NumbTargets
                pos=obj.TargProps.MeanState_allk{targid}(Tk,1:2);
                [Xf, Yf] = ds2nfu(pos(1), pos(2));
                W=0.05;
                H=0.02;
                XY=[Xf,Yf]+0.001*randn(1,2);
                annotation('textbox',[XY,W,H],'String',strcat('Target : ' ,num2str(targid)),'FontSize',7,'EdgeColor','none');
                
            end
            %
            %             axis square
            %             axis equal
            axis([obj.SimProps.Xboundary,obj.SimProps.Yboundary])
            obj.plotting_properties()
            hold off
            
        end
        %============================================================================
        function plot_truth_dyn_scenario(obj) %plot the truth trajectories over time
            figure('units','normalized','outerposition',[0 0 1 1])
            for Tk=1:obj.SimProps.Time.nsteps
                axis([obj.SimProps.Xboundary(1)-0.2*obj.SimProps.Xboundary(2),1.5*obj.SimProps.Xboundary(2),obj.SimProps.Yboundary(1)-0.2*obj.SimProps.Yboundary(2),1.5*obj.SimProps.Yboundary(2)])
                clf
                hold on
    %             Tk=obj.SimProps.Time.tk;
                % plotting the sensors and their FOV configurations
                for sensid=1:1:obj.SensProps.NumbSensors    % all sensors are blue
                    pos=obj.SensProps.SensorPos_allk{sensid}(Tk,:);
                    rmax=obj.SensProps.SensorRmax_allk{sensid}(Tk,:);
                    alp=obj.SensProps.SensorAlpha_allk{sensid}(Tk,:);
                    dr=obj.SensProps.SensorDir_allk{sensid}(Tk,:);
                    %                 keyboard
                    if strcmpi(obj.SensProps.SensorFOVtype_allk{sensid}{Tk},'1Sensor->AllTarget')   % it is a coverage sensor, plot circles
                        plot_ellipse(pos(1),pos(2),rmax,rmax,0,'b','fill');
                    end
                    if strcmpi(obj.SensProps.SensorFOVtype_allk{sensid}{Tk},'1Sensor->1Target')  % it is a track sensor, plot cones
                        plot_2Dcone(pos(1),pos(2),alp,rmax,dr,'g','fill')
                    end
                    plot(pos(1),pos(2),'ko','linewidth',2,'MarkerSize',6)

                    if strcmpi(obj.SensProps.SensorDynMotionModel_allk{sensid}{Tk},'UAVgridOnly') || strcmpi(obj.SensProps.SensorDynMotionModel_allk{sensid}{Tk},'UAVgrid+smoothTraj')
                        UAVcolors=['r','b','k','m','c'];
                        poses=obj.SensProps.SensorPos_allk{sensid}(1:Tk,:);
                        cc=obj.SensProps.SensorsID(sensid);
                        plot(poses(:,1),poses(:,2),UAVcolors(cc),'linewidth',2,'MarkerSize',6);
                    end
                end

                %plot the target with a triangle and an ellipse for confidence
                for targid=1:1:obj.TargProps.NumbTargets
                    if strcmpi(obj.TargProps.VisibilityState_allk{targid}{Tk},'Invisible') || strcmpi(obj.TargProps.VisibilityState_allk{targid}{Tk},'Lost')

                    else
                        P=reshape(obj.TargProps.CovState_allk{targid}(Tk,:),obj.TargProps.TargetModelDim_allk{targid}(Tk),obj.TargProps.TargetModelDim_allk{targid}(Tk));
                        plot_1sigellip(obj.TargProps.MeanState_allk{targid}(Tk,1:2),P(1:2,1:2),'r');
                    end


                    %                 plot_triangle(obj.TargProps.MeanState_allk{targid}(Tk,1),obj.TargProps.MeanState_allk{targid}(Tk,2),1,0.5,0,'r','fill')
                    %plot the trajectories est+truth

                    %plot triangle as mean
                    if strcmpi(obj.TargProps.TargetDynTypes_allk{targid}{Tk},'VirtualMove')  || strcmpi(obj.TargProps.TargetDynTypes_allk{targid}{Tk},'VirtualStationary') % if it is a vitual target
                        plot_triangle(obj.TargProps.MeanState_allk{targid}(Tk,1),obj.TargProps.MeanState_allk{targid}(Tk,2),1,1,0,'k','fill')

                    else
                        plot_triangle(obj.TargProps.MeanState_allk{targid}(Tk,1),obj.TargProps.MeanState_allk{targid}(Tk,2),1,1/3,0,'r','fill')

                    end
                    plot(obj.TargProps.MeanState_allk{targid}(Tk,1),obj.TargProps.MeanState_allk{targid}(Tk,2),'r','linewidth',2)
                    plot(obj.TargProps.TruePosState_allk{targid}(:,1),obj.TargProps.TruePosState_allk{targid}(:,2),'k*-','linewidth',2)

                end


                for sensid=1:1:obj.SensProps.NumbSensors
                    pos=obj.SensProps.SensorPos_allk{sensid}(Tk,:);


                    [Xf, Yf] = ds2nfu(pos(1), pos(2));
                    W=0.1;
                    H=0.05;
                    XY=[Xf,Yf]+0.03*randn(1,2);
                    annotation('textbox',[XY,W,H],'String',strcat('Sensor : ' ,num2str(sensid)),'FontSize',12,'EdgeColor','none','Color','r');

                end
                for targid=1:1:obj.TargProps.NumbTargets
                    pos=obj.TargProps.MeanState_allk{targid}(Tk,1:2);
                    [Xf, Yf] = ds2nfu(pos(1), pos(2));
                    W=0.05;
                    H=0.02;
                    XY=[Xf,Yf]+0.001*randn(1,2);
                    annotation('textbox',[XY,W,H],'String',strcat('Target : ' ,num2str(targid)),'FontSize',7,'EdgeColor','none');

                end
                %
                %             axis square
                %             axis equal
                axis([obj.SimProps.Xboundary,obj.SimProps.Yboundary])
                obj.plotting_properties()
                hold off
                pause(0.5)
                
                for targid=1:1:obj.TargProps.NumbTargets
                    obj=DynProp_Target(obj,targid,Tk);
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        function plotting_properties(obj)
            
            txtsz=15;
            set(gca,'FontSize',txtsz)
            screen_size = get(0, 'ScreenSize');
            h_plotszz = get(gca, 'title');
            kpp_plotszz = get(gca, 'xlabel');
            l_plotszz = get(gca, 'ylabel');
            m_plotszz = get(gca, 'zlabel');
            set(h_plotszz, 'FontName', 'Helvetica', 'FontSize', txtsz*1.5)
            set(kpp_plotszz, 'FontName', 'Helvetica', 'FontSize', txtsz)
            set(l_plotszz, 'FontName', 'Helvetica', 'FontSize', txtsz)
            set(m_plotszz, 'FontName', 'Helvetica', 'FontSize', txtsz)
            set(gcf, 'PaperOrientation', 'landscape');
            set(gcf, 'PaperUnits', 'centimeters');
            set(gcf, 'PaperType', 'A4');
            set(gcf, 'PaperPositionMode', 'manual');
            set(gcf, 'PaperUnits', 'inches');
            set(gcf,'PaperOrientation','landscape');
            set(gcf,'PaperUnits','normalized');
            set(gcf,'PaperPosition', [0 0 1 1]);
            
        end
        %============================================================================
        %         function PrintDebug(obj)
        %
        %
        %             SensConst=zeros(obj.Ns,1);
        %             dirs=zeros(obj.Ns,1);
        %             pos=zeros(obj.Ns,2);
        %             for i=1:1:obj.Ns
        %                 SensConst(i,1)=obj.SensConstraints{i};
        %                 dirs(i,1)=obj.dirk{i}(obj.tk);
        %                 pos(i,:)=obj.sk{i}(obj.tk,1:2);
        %             end
        %             format shortg
        %             disp([[1:1:obj.Ns]', obj.SensType , obj.hn(:) , SensConst , dirs,  pos ])
        %         end
        %============================================================================
    end
end