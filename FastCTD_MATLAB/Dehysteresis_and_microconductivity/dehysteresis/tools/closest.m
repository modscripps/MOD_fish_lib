function [closestval,ind] = closest(vec,val,numvals)
% Find the values in one vector (vec) closest to a set of target values (val)
%
% Inputs:
%   vec: A one-dimensional vector of values
%   val: A one-dimensional vector of target values
%  Optional:
%   numvals: If specified, gives the N closest values to the target value
%       val.
%
% Outputs:
%   closestval: A vector of values representing the closest value in vec to
%       the requested value in val. If numvals>1, closestval will be size
%       length(val) x length(numvals); Each row of closestval will correspond
%       to the numval closest numbers to val.
%   ind: The indices in val corresponding to the values in closestval
%
% Alex Andriatis
% 10-02-20
  if ~exist('numvals','var')
      numvals=1;
  end
  if length(numvals)<length(val)
      ind = NaN(length(val),numvals);
  else
      ind = NaN(length(val));
  end
  closestval=ind;

  for i=1:length(val)
    minvec = abs(vec-val(i));
    [~,I]=sort(minvec);
    tmp = vec(I);
    for j=1:numvals
        closestval(i,j) = tmp(j);
        ind(i,j) = I(j);
    end
  end
end

