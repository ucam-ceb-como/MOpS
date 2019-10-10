%{  
 Author:    Casper Lindberg (csl37)

 Purpose:   Draws a single TEM frame
%}
function TEMframe()

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Input parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% location of data
    coords_file = '.\test\Z1-video(3).csv';
    
    %% draw scale bar (width metres)?
    scale_bar_width = 10e-9;
    
    %% frame number to draw ( < 1 for final frame) 
    n_frame = 0; 
    
    %% half frame sizes in metres
    %% box will be [-x, x] 
    half_size_x =  50e-9;
    half_size_y =  50e-9;
    half_size_z =  50e-9;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    addpath('./matlab')
    
    %% read time, primary coordinates, and frame coordinates
    frames = read_coords(coords_file);
    n = length(frames);
    
    %% frame to draw
    if (n_frame > n || n_frame < 1) 
        n_frame = n;
    end
    disp(strcat('Drawing frame: ',num2str(n_frame),' of ',num2str(n)));
    
    %reorient frame
    new_frames = reorientate(frames, n_frame);
    
    %%draw particle
    [fig, ax] = draw_particle(new_frames(n_frame).coords,1);
    
    set(fig,'visible','on');
    
    ax.Color = [0.7 0.7 0.7];
    
    xlim([-half_size_x, half_size_x])
    ylim([-half_size_y, half_size_y])
    zlim([-half_size_z, half_size_z])
    
    %%light object
    light('Position',[-1,-1,0]);

    %draw scale bar in lower right corner
    if (scale_bar_width > 0)
        sh = 0.8;
        plot([-0.9*half_size_x, -0.9*half_size_x+scale_bar_width],[-sh*half_size_y,-sh*half_size_y],'-k','linewidth',6)
        text(-0.9*half_size_x,-(sh+0.05)*half_size_y,strcat(num2str(scale_bar_width*1e9),' nm'),'fontsize',20)
    end
end

