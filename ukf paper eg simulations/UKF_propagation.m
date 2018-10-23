function dx = UKF_propagation(t, x, model)

meu = feval(model.fx, t, x);
dx = meu;