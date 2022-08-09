function [FMAT,JMAT] = calFJ(u)
FMAT = calF(u);
JMAT = calJ(u);
end