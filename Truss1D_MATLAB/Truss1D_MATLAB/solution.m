function disp = solution(nDof,fixedDof,K,force)
%%%%% This function is to solve: K*U=F to obtain U
activeDof = setdiff((1:nDof)',fixedDof);
U = K(activeDof,activeDof)\force(activeDof);
disp = zeros(nDof,1);
disp(activeDof) = U;
end