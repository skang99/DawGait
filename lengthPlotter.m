function [dyn_data] = lengthPlotter(~,M1,M2)

% Returns an array of lengths between the dynamic marker M1, M2

[M1_x, M1_y, M1_z] = extract_XYZ(M1);
[M2_x, M2_y, M2_z] = extract_XYZ(M2);

for i = 1:length(M1_x)
    x = M1_x(i) - M2_x(i);
    z = M1_z(i) - M2_z(i);
    y(i) = sqrt(x^2 + z^2);
end

dyn_data = y;

end

