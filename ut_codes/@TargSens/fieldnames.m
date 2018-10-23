% fileds in targsens

            obj.SimProps.SimMode=[];
            obj.SimProps.SensorMsgsAvail={'OutOfFOV','PseudoUpdate','SimpleTask'};
            obj.SimProps.SensorConfigStatesAvail={'HasToTakeMeas','AlreadyTakenMeas'};
            obj.SimProps.TargetConfigStatesAvail={'CurrAt_tk','aPrioriAt_tk','aPostAt_tk'};
            obj.SimProps.SensorDynTypesAvail={'Move','Stationary'};
            obj.SimProps.TargtDynTypesAvail={'Move','Stationary','VirtualMove','VirtualStationary'};
            obj.SimProps.SigmaPntsTypeAvail={'UT','CUT4','CUT6','CUT8'};
            obj.SimProps.SensorTasksAvail={'MIUB','FIM'};
            obj.SimProps.SensorModelTypesAvail={'Range+Bearing','Range','Bearing'};
            obj.SimProps.TargetModelTypesAvail={'UM','CT','NoDyn'};
            obj.SimProps.TargetVisibilityStateAvail={'Lost','Tracked','Invisible','NotTracked'};
            obj.SimProps.SensorDynMotionModelTypes={'UAV','None'};
            obj.SimProps.SensorModeConstraintTypesAvail={'1ModeOperation'};
            
            
            obj.SimProps.Xboundary=[];
            obj.SimProps.Yboundary=[];
            
            obj.SimProps.SigmaPntFunction=[];
            
            obj.SimProps.Time.t0=t0;
            obj.SimProps.Time.tf=tf;
            obj.SimProps.Time_dt=dt;
            obj.SimProps.Time.tvec=t0:dt:tf;
            obj.SimProps.Time.nsteps=length(obj.Tvec);
            obj.SimProps.Time.tk=1; %set the curr time to 0 or Tvec(1)
            
            
            obj.TargProps.NumbTargets=0;
            obj.TargProps.TargetID=[];
            obj.TargProps.MeanState_allk=cell(1,1); %one cell for each each target
            obj.TargProps.CovState_allk=cell(1,1); %one cell for each each target
            obj.TargProps.VisibilityState_allk=cell(1,1); %one cell for each each target
            obj.TargProps.TruePosState_allk=cell(1,1); %one cell for each each target
            obj.TargProps.TargetDynTypes_allk=cell(1,1); %one cell for each each target
            obj.TargProps.Q_allk=cell(1,1); %one cell for each each target
            obj.TargProps.TargetModel_allk=cell(1,1); %one cell for each each target
            obj.TargProps.TargetModelDim_allk=cell(1,1); %one cell for each each target
            obj.TargProps.TargetsID_deleted=[];
            
            
            
            obj.SensProps.NumbSensors=0;
            obj.SensProps.SensorsID=[];
            obj.SensProps.SensorDynTypes_allk=cell(1,1); %one cell for each each sensor at each time step
            obj.SensProps.R_allk=cell(1,1); %one cell for each each sensor at each time step
            obj.SensProps.y_allk=cell(1,1,1); %rows= sensor, columns= targets, depth is time step
            obj.SensProps.SensorNmodes_allk=cell(1,1);
            obj.SensProps.SensorPlannedConfigTasks_allk=cell(1,1); %each cell for 1 sensor. 
            obj.SensProps.SensorModel_allk=cell(1,1); %one cell for each each target
            obj.SensProps.SensorModelDim_allk=cell(1,1);
            obj.SensProps.SensorFOVtype_allk=cell(1,1);
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
            obj.SensProps.SensorDynMotionModel_allk=cell(1,1);
            
            obj.SensProps.SensorPos_allk=cell(1,1);
            obj.SensProps.SensorRmax_allk=cell(1,1);
            obj.SensProps.SensorDir_allk=cell(1,1);
            obj.SensProps.SensorAlpha_allk=cell(1,1);

            
            
            