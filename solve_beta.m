syms x

b=fac*25*1.71^4;

eqn =  (b*x^4 +x^1) == 1;
solx = solve(eqn,x)

syms z

root(b*z^4 + z - 1, z, 1)