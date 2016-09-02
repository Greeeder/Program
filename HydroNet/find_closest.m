function [ difference, index ] = find_closest( value, vector )
% Find closest element to 'value' inside 'vector'.
% Return difference between 'value' and closest element and index of that
% element.

[difference, index] = min(value - vector);

end
