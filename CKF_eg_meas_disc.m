%% discrete mesuremtn function
    function yk=CKF_eg_meas_disc(x)
    zi=x(1)-0;
    ni=x(3)+0;
    
    yk(1,1)=sqrt((zi)^2+(ni)^2);
    yk(2,1)=atan2(ni,zi);
    end
