function [vectAngle] = vectorAngle(vect1,vect2)
%%%% vectorAngle %%%%

%   If there are insufficient markers to create rotation matrices then this
%   will provide simple angles between two separate vectors created from two different markers

dotProd = dot(vect1',vect2');
vectAngle = 90 - acosd(dotProd);
end

