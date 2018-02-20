function S = rotBlochVector(gamma,theta)
%   evaluates spin-up Bloch vector rotated around torque in XY-plane
%
%   S: 3x1 Bloch vector (x,y,z); z-axis is quantisation   
%
%   gamma: torque vector tilt angle from X-axis
%   theta: rotation angle (rabi freq * duration)
%

if isnan(gamma) || isnan(theta)
    S=NaN([3,1]);
    return
end

S0 = [0;0;1];     % spin-up vector

phi_y=rad2deg(-(pi/2-gamma));       % y-rot angle to align torque to z (degrees)

Stemp=roty(phi_y)*S0;       % y-rot align torque to z
Stemp=rotz(rad2deg(theta))*Stemp;           % precession
S=roty(-phi_y)*Stemp;       % back to original coords

end