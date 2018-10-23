%MPTEST
% This script contains various tests of the mp*.m files
%  Note that MPDER, MPSHIFT, MPNEXT  are tested incidentally,
%  with the major stress on MPVAL, MPMAK, MPAPI.

% cb aug97; 24jan99 (added test for MPDIR)

% try trivial stuff first

% univariate polynomial  p(x) = 1-x;
mp = mpmak([1 -1],1);
if(any(mpval(mp,[-1:2])-[2:-1:-1]))
   fprintf('error in test1.\n'), end
if any(mpbrk(mpder(mp,1),'co')-(-1))
   fprintf('error in test1der.\n'), end

% univariate polynomial  p(x) = (1-x)^2;
mp = mpmak([1 -2  1],1);
if(any(mpval(mp,[-1:2])-(1-[-1:2]).^2))
   fprintf('error in test2.\n'), end

% univariate polynomial  p(x) = (1-x)^3;
mp = mpmak([1 -3 3 -1],1);
if(any(mpval(mp,[-1:2])-(1-[-1:2]).^3))
   fprintf('error in test3.\n'), end
if any(mpbrk(mpder(mp,-2),'co')-[6 -12 6])
   fprintf('error in test3der.\n'), end

% a bivariate polynomial  p(x,y) = (1-x)^2;
mp = mpmak([1 0 -2 0 0 1],2);
if(any(mpval(mp,dcube(2))-[0 0 4 4]))
   fprintf('error in test4.\n'), end

% a trivariate quadratic that vanishes at all the corners of the cube
mp = mpmak([-3 0 0 0 1 0 1 0 0 1],3);
if(any(mpval(mp,dcube(3))))
   fprintf('error in test5.\n'), end
p = [1;-2;3]; xi = [2;1;-2]; h = 1.e-7;
if abs(mpval(mp,[p+h*xi p-h*xi])*[1;-1]/(2*h) - mpval(mpder(mp,xi),p))>h
   fprintf('error in test5der.\n'), end

% test also change of variable, on the same example, using rotation by 2pi/3
% three times around a randomly chosen axis:
a = rand(3,1); a = a/norm(a); phi = (2/3)*pi; c = cos(phi); s = sin(phi);
A = c*eye(3) + (1-c)*(a*a') + s*[0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0];
mpn = mpcov(mpcov(mpcov(mp,A),A),A);
if any(abs(mpbrk(mpn,'coefs')-mpbrk(mp,'coefs'))>3e-15)
   fprintf('error in test5c.\n'), end

% test also change of form, on the same example, comparing values at
% some randomly chosen points:
points = (rand(3,10)-.5)*4;
if any(abs(mpval(mp,points)-mpval(mp2fm(mp),points))>3e-15)
   fprintf('error in test5v.\n'), end

% finally, check it out when there is less symmetry
% i.e., the polynomial (x-1)(y-1)(z-1) which vanishes at all but one corner
mp = mpmak([-1 1 1 1 0 -.5 0 -.5 -.5 0 0 0 0 0 0 1/6 0 0 0 0],3);
if(abs(mpval(mp,dcube(3))-[0 0 0 0 0 0 0 -8])>1.e-15)
   fprintf('error in test6.\n'), end
p = [1;-2;3]; xi = [2;1;-2]; h = 1.e-5;
if abs(mpval(mp,[p+h*xi p p-h*xi])*[1;-2;1]/(h*h) - ...
                                   mpval(mpder(mp,[xi xi]),p))>10*h
   fprintf('error in test6der.\n'), end
if mpbrk(mpder(mp,[xi -2*xi xi -xi 3*xi]),'co')
   fprintf('error in test6ader.\n'), end
% Now check the vector-valued calculation:
%  (x-1)(y-1)(z-1)
%  (x+1)(y-1)(z-1)
%  (x-1)(y+1)(z-1)
%  (x-1)(y-1)(z+1)

cc=  [-1 +1 +1 +1 0 -.5 0 -.5 -.5 0 0 0 0 0 0 1/6 0 0 0 0
  1 -1 -1 +1 0 +.5 0 -.5 -.5 0 0 0 0 0 0 1/6 0 0 0 0
  1 -1 +1 -1 0 -.5 0 +.5 -.5 0 0 0 0 0 0 1/6 0 0 0 0
  1 +1 -1 -1 0 -.5 0 -.5 +.5 0 0 0 0 0 0 1/6 0 0 0 0];
mp = mpmak(cc,3);
dd = [0 0 0 0 0 0 0 -8
 0 0 0 8 0 0 0 0
 0 0 0 0 0 8 0 0
 0 0 0 0 0 0 8 0];
if(any(abs(mpval(mp,dcube(3))-dd)>1.e-15))
   fprintf('error in test7.\n'), end

% Now we are ready to tackle the least interpolant...
% first, the simplest case
mp = mpapi(0,0);
if(mpval(mp,0))
   fprintf('error in test8.\n'), end

% try bivariate, linear interpolation
points = [0 1 0;0 0 1]; values = [0 1 1];
mp = mpapi(points, values);
if(mpval(mp,points) - [0 1 1])
   fprintf('error in test9.\n'), end

% try bivariate, quadratic, still simple
points = [0 1 2 0 1 0; 0 0 0 1 1 2];
coefs = [1 -1 1 -1 1 -1];
mp0 = mpmak(coefs,2);
if any(abs(mpbrk(mpapi(points,mpval(mp0,points)),'co')-coefs)>1.e-14)
   fprintf('error in test10.\n'), end

% try bivariate, cubic, still simple
points = [0 1 2 3 0 1 2 0 1 0;0 0 0 0 1 1 1 2 2 3];
coefs = [1 -1 1 -1 1 -1 1 -1 1 -1];
mp0 = mpmak(coefs,2);
if any(abs(mpbrk(mpapi(points,mpval(mp0,points)),'co')-coefs)>1.e-14)
   fprintf('error in test11.\n'), end

% try trivariate interpolation to a linear function
points = [0 1 0 0 1; 0 0 1 0 1; 0 0 0 1 1];
coefs = [1 -2 3 -1];
mp0 = mpmak(coefs,3);
values = mpval(mp0, points);
if any(abs(mpbrk(mpapi(points, values, 1.e-14, [0;0;0]),'co')- ...
                                              [coefs 0 0 0 0 0 0])>1.e-14)
   fprintf('error in test12.\n'), end

% try trivariate, quadratic, interpolation to a quadratic function
points = [0 1 2 0 1 0 0 1 0 0;0 0 0 1 1 2 0 0 1 0;0 0 0 0 0 0 1 1 1 2];
coefs = [1 -1 1 -1 1 -1 1 -1 1 -1];
mp0 = mpmak(coefs,3);
values = mpval(mp0, points);
mp = mpapi(points, values);
mps = mpshft(mp,[0;0;0]);
if any(abs(mpbrk(mps,'co')-coefs)>1.e-14)
   fprintf('error in test13.\n'), end
if any(abs(mpbrk(mp,'co') - mpbrk(mpshft(mps,mpbrk(mp,'ce')),'co'))>1.e-14)
   fprintf('error in test13shift.\n'), end
 
% try trivariate, cubic, interpolation to a quadratic function
points = dcube(3);
coefs = [1 -1 1 -1 1 -1 1 -1 1 -1];
mp0 = mpmak(coefs,3);
values = mpval(mp0, points);
% In this example, the coefficients are not reproduced, but the interpolant is
% in Pi_2. To check it
cc = mpbrk(mpapi(points, values),'co');
cc(1:10) = cc(1:10)-coefs;
mpd = mpmak(cc,3);
if any(abs(mpval(mpd,points))>1.e-14)
   fprintf('error in test14.\n'), end

% test the re-use of factorization information:
points = rand(3,10); values = [1 -2 1]*(points.*(points - 1));
[mp, mpfactor] = mpapi(points,values,1.e-12,[0;0;0]);
mpf = mpapi(mpfactor,values);
if (max(abs(mpbrk(mp,'co')-mpbrk(mpf,'co'))))  
   fprintf('error in test15.\n'), end

% test the vector capability:
points = rand(3,10); values = (points.*(points - 1));
mp = mpapi(points,values);
if any(abs(values - mpval(mp,points))>1.e-14)
   fprintf('error in test16.\n'), end

% test the noise level
%fprintf(['The next test checks interpolation at 10 data sites along a',...
%     ' straight line in R^3 to quartic values.\n',...
%  'The interpolant should formally be of degree 9 but its coefficients\n',...
%  'of degree > 4 should be zero. Also, the derivative in any direction\n',...
%  'perpendicular to that straight line should be zero.\n',...
%'The next three numbers give the size of the supposedly zero coefficients\n',...
%'relative to the size of the coefficients of the interpolant.\n'])
t = 6*rand(1,10)-3; xi = [1;2;3]; points = xi*t; values = t.^4;
mp = mpapi(points,values);
% coefficients of degree > 4 should be zero:
coefs = mpbrk(mp,'co'); blaise = mpbrk(mp,'b');
mc = max(abs(coefs)); mc5 = max(abs(coefs(blaise(end,6)+1:blaise(end,11))));
if mc5/mc>5e-15
fprintf('in test17, noise ratio = %8g \n', mc5/mc), end
% Also, the interpolant should be constant in directions orthogonal to  xi :
noise1 = max(abs(mpbrk(mpder(mp,[1;-2;1]),'co')));
noise2 = max(abs(mpbrk(mpder(mp,[3;0;-1]),'co')));
if any([noise1,noise2]>5e-15)
   fprintf('in test17, derivative noise = %8g, %8g \n',noise1/mc, noise2/mc)
end

% % Now add a speed test, to compare the old bivariate routine with the present
% % one, on a large number of points:
% nn = 8;
% for n=1:4
%    nn = nn*2
%    points = rand(2,nn); values = mpval(mpmak([1 -1 1 -1 1 -1],2),points);
%    tic
%    mp = mpapi(points,values,1.e-8);
%    toc
%    tic
%    bp = bpapi(points,values);
%    toc
% end

% some tests of MPDIR, ending in a nice picture

%  test 18
% make up a linear interpolant
p = mpapi([0 1 0;0 0 1],eye(3));
% this gives the vector valued interpolant with linear terms
%   -y -x
%       x
%    y
%  hence with linear coefficients
linco = [-1 -1
   0  1
   1  0];

% mpbrk(p);
dp = mpdir(p,[0 1;1 0]);
if any(linco - reshape(mpbrk(dp,'coefs'),3,2)),
   fprintf('trouble with the test 18\n')
end

% test 19 
% compare it with MPDER, using the trivariate example used in test7
mp = mpmak(cc,3); y = [1;2;5]; derp = mpder(mp,y); dirp = mpdir(mp,y);
if max(max(abs(mpbrk(derp,'coe')-mpbrk(dirp,'coe')))>1e-6)
   fprintf('trouble with the test 19\n')
end

% test 20 
% make up the monkey saddle get a contour plot, then superpose some gradients

[1:6]*(pi/3); p = mpapi([cos(ans);sin(ans)],[-1 1 -1 1 -1 1]);
% p = mpmak([0 0 0 1 0 1],2);
x = -1:.1:1; [xx,yy] = meshgrid(x,x); z = zeros(size(xx));
z(:) = mpval(p,[xx(:) yy(:)].');
contour(x,x,z,21), hold on, axis equal
title('level lines and gradients for the Monkey Saddle')
xlabel('test 20 for MPDIR')
dp = mpdir(p,eye(2));
grads = reshape(mpval(dp,[xx(:) yy(:)].'),[2,length(xx),length(yy)]);
quiver(x,x,squeeze(grads(1,:,:)),squeeze(grads(2,:,:))), hold off
