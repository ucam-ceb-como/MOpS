%{  
 Author:    Casper Lindberg (csl37)

 Purpose:   Creates a video from MOpS "-videos().csv" file format

 To do:     Add composition change
%}
function TEMvideo()

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Input parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% location of data
    coords_file = '.\test\Z1-video(3).csv';
    
    %% video parameters
    frame_rate = 20;
    
    %% half frame sizes in metres
    %  box will be [-x, x] 
    half_size_x =  50e-9;
    half_size_y =  50e-9;
    half_size_z =  50e-9;
    
    %% draw scale bar (width in metres)
    scale_bar_width = 10e-9;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    addpath('./matlab');
    close all;
    
    %% read time, primary coordinates, and frame coordinates
    frames = read_coords(coords_file);
    n_frames = length(frames);
    disp(strcat('Frames:',num2str(n_frames)));
    
    %reorientate frame
    new_frames = reorientate(frames, [1:n_frames]);
    
    %%create video structure
    v = VideoWriter('video.mp4','MPEG-4');
    v.FrameRate = frame_rate;
    open(v);

    %%loop over frames
    for k = 1:n_frames
        
        %%draw particle
        [fig, ax] = draw_particle(new_frames(k).coords);

        ax.Color = [0.7 0.7 0.7];
            
        %%light object
        light('Position',[0,0,-1]); %light from underneath

        xlim([-half_size_x, half_size_x])
        ylim([-half_size_y, half_size_y])
        zlim([-half_size_z, half_size_z])

        %draw scale bar in lower right corner
        if (scale_bar_width > 0)
            sh = 0.8;
            plot([-0.9*half_size_x, -0.9*half_size_x+scale_bar_width],[-sh*half_size_y,-sh*half_size_y],'-k','linewidth',6)
            text(-0.9*half_size_x,-(sh+0.05)*half_size_y,strcat(num2str(scale_bar_width*1e9),' nm'),'fontsize',20)
        end

        %add clock
        text(-0.9*half_size_x,0.9*half_size_y,sprintf('%0.6f s',new_frames(k).t),'fontsize',20)
        
        %%write movie frame
        writeVideo(v,getframe(gca))
        close(fig); %close frame
    end
    
    close(v);   %close video object
end