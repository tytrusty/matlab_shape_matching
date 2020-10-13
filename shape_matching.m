function shape_matching

    d=3;                % dimension
    dt = 0.01;          % timestep
    alpha = 1;          % stiffness for shape matching
    beta = 0.95;        % blending parameter for A & R in shape matching
    k_spring = 10;      % spring stiffness
    gravity = 0;        % gravity force magnitude
    
    % How much to pad rotation matrix when using quadratic motion
    padnum = 3;         
    
    % Spring positions (two springs)
    spring_a = [0.0 2];
    spring_b = [1 0];
    
    % Gravity force vector
    f_gravity = [0 gravity];
    
    if d == 3
        padnum=6;
        spring_a = [spring_a 0];
        spring_b = [spring_b 0];
        f_gravity = [f_gravity 0];
    end
    
    % Load mesh
    [V,F] = readOBJ('sphere.obj');
    
    % Transform shape a little
    V = V / 4;
    V = V + [0 1 0];
    V = V(:,1:d);
    
    % Plot mesh
    figure(2);
    p1 = tsurf(F, V);
    view([0 90])     
    xlim([-0.7 2.5]);
    ylim([0 2.5]);
    %axis equal
    %shading interp
    lighting gouraud
    lightangle(gca,-130,90)

    x = V';
    v = zeros(size(x));
    
    x0_com = mean(x,2); % assuming constant point masses
    q = x - x0_com;
    
    % Computing quadratic terms (Section 4.3)
    if d == 2
        q_ = zeros(5, size(q,2));
        q_(1:2,:) = q;
        q_(3:4,:) = q.^2;
        q_(5,:) = q(1,:).*q(2,:);
    else
        q_ = zeros(9, size(q,2));
        q_(1:3,:) = q;
        q_(4:6,:) = q.^2;
        q_(7,:) = q(1,:).*q(2,:);        
        q_(8,:) = q(2,:).*q(3,:);        
        q_(9,:) = q(3,:).*q(1,:);        
    end
    
    Aqq_ = inv(q_ * q_');
 
    disp('Press ANY key');
    waitforbuttonpress
    
    % Simulation loop
    for t=0:dt:30
        x_com = mean(x, 2);
        p = x - x_com;
       
        Apq = p * q';
        Apq_ = p * q_';
        
        % Polar decomposition (Section 3.5)
        [U, e] = eig(Apq'*Apq);
        e = 1 ./ sqrt(diag(e));
        Sinv = repmat(e',size(d,1),1).*U*U';
        R = Apq*Sinv;
        
        % Computing A & R matrix for quadratic deformation
        R_ = padarray(R,[0,padnum],0,'post');
        A_ = Apq_ * Aqq_; 
        A_ = A_ / sqrt(det(A_*A_')^(1/d)); % preserve volume
        AR_ = (1-beta)*R_ + beta*A_;       % blending R & A matrices
        
        % Computing goal positions
        g = AR_*q_ + x_com;
        
        % Computing spring forces
        f_spring = zeros(size(v));
        f_spring(:,F(5,:)) = -k_spring * (x(:,F(5,:)) - spring_a);
        f_spring(:,F(7,:)) = -k_spring * (x(:,F(7,:)) - spring_b);
        
        % Time integration
        v = v + alpha*(g-x) + dt*f_gravity' + dt*f_spring;
        x = x + dt*v;
        
        % Update plot
        p1.Vertices = x';
        drawnow;
    end
    
end