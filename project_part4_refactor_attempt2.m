% Problem Statement: Design a 2-D inlet for a turbojet engine with 3
% oblique shocks and a normal shock
% Procedure: Set the 3 deflection angles to be equal and find the one that
% results in the greatest stagnation pressure ratio. Then tune the first
% angle so that the stagnation pressure loss across the first two shocks is
% equal. Then do the same for the third deflection angle (d1 < d2 < d3)

% In order to reduce runtime,  the optimal angle is calculated over a range
% with a relatively large resolution. The resolution is then decreased and
% used to create a new range around the previously calculated angle.
% The optimal angle is recalculated over the new range. This is repeated
% until the desired resolution is met.

% Reference for obliqueshock()
% Mohamed Hassan (2022). Oblique Shock relations solver function
% (https://www.mathworks.com/matlabcentral/fileexchange/100898-oblique-shock-relations-solver-function),
% MATLAB Central File Exchange. Retrieved December 8, 2022.

clear all;
clc;

gamma = 1.4;
m_flight = 2.1;

angle_res = 0.5;
current_res = deg2rad(angle_res);
desired_res = deg2rad(0.01);
deflection_range = [8:current_res:9.5];
deflection_range = deg2rad(deflection_range);

[~, ~, p1_po1, ~, ~] = flowisentropic(gamma, m_flight);


% Finding optimal value for second deflection angle
while current_res > desired_res
    % Preallocate memory
    [theta, mn_pre, mn_post] = deal(zeros(3, length(deflection_range)));
    p_ratio = zeros(4, length(deflection_range));
    stagp_ratio = ones(1, length(deflection_range));
    p5_po5 = zeros(1, length(deflection_range));
    m = zeros(5, length(deflection_range));

    for i=1:length(deflection_range)
        m(1, i) = m_flight;
    end

    % calculate flow properties after each shock
    for i=1:length(deflection_range)%for d_num=1:3
        for d_num=1:3%for i=1:length(deflection_range)
            % deflection oblique shock relations
            [~,theta(d_num, i),~,~,~,~,~,~,~] = obliqueshock(3,1,m(d_num, i),deflection_range(i));
            theta(d_num, i) = double(theta(d_num, i));
            mn_pre(d_num, i) = m(d_num, i) * sin(theta(d_num, i));
    
            % deflection normal shock relations
            [~, ~, p_ratio(d_num, i), ~, mn_post(d_num, i), ~, ~] = flownormalshock(gamma, mn_pre(d_num, i));

            % calc mach number after shock
            m(d_num+1, i) = mn_post(d_num, i)/sin(theta(d_num, i) - deflection_range(i));
        end
        % calculate properties after terminating normal shock
        [~, ~, p_ratio(4, i), ~, m(5, i), ~, ~] = flownormalshock(gamma, m(4, i));
        % calculate isentropic relations for final flows
        [~, ~, p5_po5(i), ~, ~] = flowisentropic(gamma, m(5, i));

        % calculate stagnation pressure ratios for each delfelction angle
        for j=1:size(p_ratio, 1)
            stagp_ratio(i) = stagp_ratio(i) * p_ratio(j, i);
        end
        stagp_ratio(i) = stagp_ratio(i) * p1_po1 * (1/p5_po5(i));
    end
    % find max stagnation pressure ratio
    k = find(stagp_ratio==max(stagp_ratio));
    delta2 = deflection_range(k);
    new_res = current_res / 10;
    deflection_range = [deflection_range(k)-current_res:new_res:deflection_range(k)+current_res];
    current_res = new_res;
end
%%
% Finding optimal value for first deflection angle by distibuting
% stagnation pressure loss equally across deflections 1 and 2
current_res = deg2rad(angle_res);
deflection_range = [delta2-(2*current_res):current_res:delta2+(2*current_res)];
while current_res > desired_res
    % Preallocate memory
    m = zeros(3, length(deflection_range));
    [p_ratio, stagp_ratio, theta, mn_pre, mn_post] = deal(zeros(2, length(deflection_range)));
    [p2_po2, p3_po3, po2_po1, po3_po2, stagp_diff] = deal(zeros(1, length(deflection_range)));
    
    for i=1:length(deflection_range)
        m(1, i) = m_flight;
    end

    for i=1:length(deflection_range)
        for j=1:2
            if j==1
                angle = deflection_range(i);
            else
                angle = delta2;
            end
            % oblique shock to normal shock relations
            [~,theta(j, i),~,~,~,~,~,~,~] = obliqueshock(3,1,m(j, i),angle);
            theta(j, i) = double(theta(j, i));
            mn_pre(j, i)  = m(j, i) * sin(theta(j, i));
    
            % deflection normal shock relations
            [~, ~, p_ratio(j, i), ~, mn_post(j, i), ~, ~] = flownormalshock(gamma, mn_pre(j, i));
    
            % calc mach number after shock
            m(j+1, i) = mn_post(j, i)/sin(theta(j, i) - angle);
        end
        [~, ~, p3_po3(i), ~, ~] = flowisentropic(gamma, m(3, i));
        [~, ~, p2_po2(i), ~, ~] = flowisentropic(gamma, m(2, i));
        po2_po1(i) = p1_po1 * p_ratio(1, i) * (1/p2_po2(i));
        po3_po2(i) = p2_po2(i) * p_ratio(2, i) * (1/p3_po3(i));
        stagp_diff(i) = abs(po2_po1(i) - po3_po2(i));
    end
    k = find(stagp_diff == min(stagp_diff));
    delta1 = deflection_range(k);
    new_res = current_res / 10;
    deflection_range = [deflection_range(k)-current_res:new_res:deflection_range(k)+current_res];
    current_res = new_res;
end
m2 = m(2, k);
m3 = m(3, k);
po3_po2_optimal = po3_po2(k);
p3_po3_optimal = p3_po3(k);
%%
% Finding optimal value for third deflection angle by distributing
% stagnation pressure loss equally across all deflections
current_res = deg2rad(angle_res);
deflection_range = [delta2:current_res:delta2+(2*current_res)];
while current_res > desired_res
    % Preallocate memory
    [m, mn_pre, mn_post, p4_po4, p_ratio, theta, po4_po3, stagp_diff] = deal(zeros(1, length(deflection_range)));

    for i=1:length(deflection_range)
        % oblique shock relations
        [~,theta(i),~,~,~,~,~,~,~] = obliqueshock(3,1,m3,deflection_range(i));
        theta(i) = double(theta(i));
        mn_pre(i)  = m3 * sin(theta(i));

        % deflection normal shock relations
        [~, ~, p_ratio(i), ~, mn_post(i), ~, ~] = flownormalshock(gamma, mn_pre(i));

        % calc mach number after shock
        m(i) = mn_post(i)/sin(theta(i) - deflection_range(i));

        [~, ~, p4_po4(i), ~, ~] = flowisentropic(gamma, m(i));
        po4_po3(i) = p3_po3_optimal * p_ratio(i) * (1/p4_po4(i));
        stagp_diff(i) = abs(po3_po2_optimal - po4_po3(i));
        k = find(stagp_diff == min(stagp_diff));
        delta3 = deflection_range(k);
        new_res = current_res / 10;
        deflection_range = [deflection_range(k)-current_res:new_res:deflection_range(k)+current_res];
        current_res = new_res;
    end
end
%%
fprintf('delta1 = %.2f degrees\ndelta2 = %.2f degrees\ndelta3 = %.2f degrees', rad2deg(delta1), rad2deg(delta2), rad2deg(delta3));