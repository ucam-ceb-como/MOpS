
function mov = video()

    close all;
    
    %% location of data
    mops_path = '.\';           %%root folder
    coords_file = strcat(mops_path,'\Part-Coords().csv');
    
    %%read coordinates into 3D array
    [t, coords] = read_coords(coords_file);
    [L,m,n] = size(coords);
    
    %%create movie structure

    mov(:) = struct('cdata',[],'colormap',[]);
    set(gca,'nextplot','replacechildren');
    
    v = VideoWriter('movie.avi');
    v.FrameRate = 15;
    open(v);

    nframe = 1;
    box_size =  [-120e-9 120e-9];
    k=200; %test
    %for k = 60:n
        %if (or(k==1, isequal(coords(:,:,k),coords(:,:,k-1))==0))
     %   if (mod(k,2)==0)
            %%draw particle
            fig = draw_particle(coords,k);
            xlim(box_size)
            ylim(box_size)
            zlim(box_size)
            %%light object
            light('Position',[-1,-1,1]);
            %%get movie frame
            mov(nframe) = getframe(fig);
            writeVideo(v,getframe(gca))
            nframe = nframe + 1;
      %  end
    %end
        
    close(v)
end

%%draws a particle
function [fig] = draw_particle(coords,k)
    
    figure()
    hold on
    
    [L,m,n] = size(coords);
    for i = 1:L
        %%loop over primaries
        s = draw_primary(coords(i,:,k));
        set(s,'FaceColor',[0.7 0.7 0.7],'FaceAlpha',1,'EdgeColor' ,'none');
    end
    colormap gray
    fig = gcf;
end

%%draws a primary particle
function [s] = draw_primary(coords)

[x0,y0,z0] = sphere(25); %coordinates of 20x20 face unit sphere at (0,0,0) %sphere(N) specifies N faces

%rescale and recentre sphere 
x1 = x0*coords(4) + coords(1);
y1 = y0*coords(4) + coords(2);
z1 = z0*coords(4) + coords(3);

s = surf(x1,y1,z1);

end

%%read primary coordinates time series 
%%coords is a 3D array with i: primary, j: x,y,z,r coordinates, k: time
%%In: file_name
%%Out: times(k), coordinates(i,j,k)
function [t, coords] = read_coords(file_name)

    %%read data 
    data = csvread(file_name,1,0);
    [m,n] = size(data);
    
    %%first time coordinate
    t(1) = data(1,1);   
    %%initialize counters
    i=1;
    k=1;
    %%loop over rows of data (2D) and transfer to coords (3D)
    for r = 1:m
        if data(r,1) ~= t(k)
            i = 1;
            k = k+1;
            t(k) = data(r,1);
        end
        coords(i,1:4,k) = data(r,2:5);
        i=i+1;
    end
    
end
