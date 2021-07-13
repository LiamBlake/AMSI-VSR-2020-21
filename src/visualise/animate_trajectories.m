function animate_trajectories(traj,style,delay,omit_const)

    nT = size(traj,3);
    
    if omit_const
        % Do not plot constant trajectories
        constx = sum(diff(traj(:,:,:,1)) == 0, 3) == nT;
        consty = sum(diff(traj(:,:,:,2)) == 0, 3) == nT;
        
        
        land = logical(repmat(constx.*consty, 1,1,1,2));
        
        traj(land) = nan;
        
    end

    figure;
    for t = 1:nT
        plot(squeeze(traj(:,:,t,1)), squeeze(traj(:,:,t,2)), style);
        pause(delay)
    end
end