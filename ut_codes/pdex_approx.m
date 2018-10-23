function uend = pdex_approx(time, domain, model, Ttmp, GF_px)

global modelQ dom_x val_x
modelQ = model.Q;
val_x = GF_px;

% grid space and time
dom_x = domain.x;
t = linspace(Ttmp(1),Ttmp(2),time.dt_N/10);

% solve pde
m = 0;
sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,dom_x,t);

% Extract the first solution component as u.
u = sol(:,:,1);
uend = u(end,:) ./ trapz(dom_x, u(end,:));


% --------------------------------------------------------------
function [c,f,s] = pdex1pde(x,t,u,DuDx)

global modelQ

c = 1;
f = 1/2*modelQ*DuDx;
% f = 0;
s = -sin(x)*DuDx - u*cos(x);

% --------------------------------------------------------------
function p0 = pdex1ic(x)

global dom_x val_x

p0 = interp1(dom_x, val_x, x);

% --------------------------------------------------------------
function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t)
pl = ul;
ql = 0;
pr = ur;
qr = 0;