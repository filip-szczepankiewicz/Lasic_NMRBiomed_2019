function [V, L] = mcw19_fB(q, dt)

[V, L] = eig(q'*q*dt);