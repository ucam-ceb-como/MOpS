%{  
 Author:    Casper Lindberg (csl37)

 Purpose:   Create TEM-style images from a MOpS "-primary.csv" file
            
            Expects  a csv file with primary data in columns
                1:  Particle Index
                7:  Position x
                8:  Position y
                9:  Position z
                10: Primary radius (m)
%}
function TEMimage()
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  Input parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% location of data
    coords_file = './videos/test/Z1-primary.csv';
    
    %% output
    show_image = 'on';          %make image visible
    output = 1;                 %save TEM?: yes(1)/no(0)
    folder_out = './test-TEM';  %output folder name
    
    %%
    scale_bar_width = 50e-9;    %scale bar width (metres)
    frames = 5;                 %number of frames
    max_coverage = 0.05;        %maximum coverage fraction
    particles = 200;            %max particles per frame
    
    %% frame size
    x_box = 200e-9;             %image half width and height (m)
    y_box = 150e-9;
    x_pix = 1200;               %image width and height (pixels)
    y_pix = 800;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    close all;
        
    %%read coordinates into 3D array
    [coords] = read_coords(coords_file);
    [~,~,n_part] = size(coords);
    
    %%create output folder
    if (output == 1)
        mkdir(folder_out);
    end
    
    %generate a random list of particle indices without replacement
    index_list = datasample([1:n_part],n_part,'replace',false);
    
    %initialise some variables
    stop = false;
    index_count = 1;
    
    %loop over frames
    for i = 1:frames
        
        %create figure
        fig = figure('units','pixels','Position', [0 0 20+x_pix 20+y_pix],'visible',show_image);
        hold on
        xlim([-x_box x_box]);
        ylim([-y_box y_box]);
        zlim([-x_box x_box]);
               
        %loop over particles in figure
        for j = 1:particles
             
            %%particle index
            if(index_count > n_part)
               %%not enough particles: stop 
               disp(strcat('Not enough particles, stopping after ',num2str(i-1),' frames'));
               stop = true;
               break;
            else
                k = index_list(index_count);
                index_count = index_count + 1;
            end
            
            %%generate particle position
            position = [(rand-1/2)*2*x_box, (rand-1/2)*2*y_box,0] ; 
            %%draw particle at given x,y,z position
            s = draw_particle(position,coords,k);
            
            %%check coverage
            ax = gca;
            ax.Color = [1 1 1];
            fi = getframe(ax);
            fiBW = imbinarize(fi.cdata,0.99);
            coverage = 1-sum(sum(fiBW(:,:,1)))/numel(fiBW(:,:,1));
           
            %if coverage exceeds threshold remove the last particle and
            %break
            if (coverage > max_coverage)
                for it = i:numel(s)
                   delete(s(it));
                end 
                break;
            end            
        end
        
        %if not enough particle stop without saving
        if (stop)
            break;
        end
        
        %%set background colour
        ax.Color = [0.6 0.6 0.6];
                
        %draw scale bar in lower left corner
        if (scale_bar_width > 0)
            sh = 0.8;
            plot([-0.9*x_box, -0.9*x_box+scale_bar_width],[-sh*y_box,-sh*y_box],'-k','linewidth',6)
            text(-0.9*x_box,-(sh+0.05)*y_box,strcat(num2str(scale_bar_width*1e9),' nm'),'fontsize',20)
        end
        
        ax.Color = [0.6 0.6 0.6];
        
        %scale image
        set(ax,'units','pixels','position',[0 0 x_pix y_pix]);
        
        set(ax,'xtick',[]);
        set(ax,'ytick',[]);
        
        %%light object
        light('Position',[0,0,-1]); %light from underneath
   
        if (output == 1)
            %%save frame
            fr=getframe; 
            imwrite(fr.cdata,strcat(folder_out,'/TEM-',num2str(i),'.png'));
        end
        
    end
    %close figure
    close(fig);
end

%%draws a particle
function s = draw_particle(position,coords,k)
     
    hold on
    
    [L,m,n] = size(coords);
    
    %construct rotation matrix
    %Implementation of: Arvo, J., Fast Random Rotation Matrices in Graphics
    %Gems III, Academic Press, 1992
    %M = (2VV^T - I)R
    %R is a rotation about the vertical axis 
    %(2VV^T - I) rotates the north pole to a random position
    
    %random numbers
    x1 = 2*pi*rand;
    x2 = 2*pi*rand;
    x3 = rand;
    
    %identity matrix
    I = eye(3);
  
    %random rotation matrix about z axis
    R(1,1) = cos(x1);
    R(1,2) = sin(x1);
    R(1,3) = 0;
    R(2,1) = -sin(x1);
    R(2,2) = cos(x1);
    R(2,3) = 0;
    R(3,1) = 0;
    R(3,2) = 0;
    R(3,3) = 1;
    
    %now rotate north pole to random position
    V(1,1)= sqrt(x3)*cos(x2);
    V(2,1)= sqrt(x3)*sin(x2);
    V(3,1)= sqrt(1-x3);
    
    G = 2*mtimes(V,transpose(V)) - I;
    
    %final rotation matrix
    M = mtimes(G,R);
    
    %%Apply rotation to primaries 
    %%loop over primaries
    for i = 1:L
        if (coords(i,4,k) > 0 )  
            %%rotate coordinates
            coords(i,1:3,k) = transpose(mtimes(M,transpose(coords(i,1:3,k))));
            %%translate coordinates    
            coords(i,1:3,k) = translate_coords(coords(i,1:3,k),position);
        end
    end
    
    %vector of primary object handles
    s = [];
    
    for i = 1:L      
        %check if primary (radius > 0)
        if (coords(i,4,k) > 0 )          
            
            %%draw primary 
            p = draw_primary(coords(i,:,k));
            
            %%check for overlap and truncate sphere
%            p = check_overlap(p,i,coords(:,:,k));
            
            set(p,'FaceColor',[0.6,0.6,0.6],'FaceAlpha',0.5,'EdgeColor' ,'none');             
            s = [s p];
        end
    end
end

%% translate coordinates to new centre
function new_coords = translate_coords(coords,centre)
    
    new_coords(1,1) = coords(1) + centre(1);
    new_coords(1,2) = coords(2) + centre(2);
    new_coords(1,3) = coords(3) + centre(3);
end 

%% draws a primary particle
function s = draw_primary(coords)

    [x0,y0,z0] = sphere(50); %coordinates of 50x50 face unit sphere at (0,0,0)

    %rescale and recentre sphere and translate position 
    x1 = x0*coords(4) + coords(1);
    y1 = y0*coords(4) + coords(2);
    z1 = z0*coords(4) + coords(3);

    s = surf(x1,y1,z1);

end

%% read primary coordinates time series 
%{
 In:   file_name
 Out:  3D coords array with: primary; x,y,z,r coordinates; particle
%}
function [coords] = read_coords(file_name)

    %%read data 
    data = csvread(file_name,1,0);
    [m,~] = size(data);
       
    %%initialize counters
    i=1;    %particle number counter
    k=1;    %primary number counter
    %%loop over rows of data (2D) and transfer to coords (3D)
    for r = 1:m
        ii = data(r,1); %particle number counter (in file)
        coords(k,1,i) = data(r,7);
        coords(k,2,i) = data(r,8);
        coords(k,3,i) = data(r,9);
        coords(k,4,i) = data(r,10);

        if r < m
            if data(r+1,1)~=ii %next particle
                k = 1;      %reset primary counter
                i = i+1;    %increment particle counter
            else
                k = k+1; %increment primary counter
            end
        end
    end
end

%%function to check for overlaps and truncate sphere
%%In: surface coordinates x,y,z, primary index, list of all primaries
function s=check_overlap(s,primary,primaries)

    [L,~] = size(primaries);
    %%loop through list of primaries checking for overlap
   
    x1 = primaries(primary,1);
    y1 = primaries(primary,2);
    z1 = primaries(primary,3);
    r1 = primaries(primary,4);
    
    for i = 1:L
        if i ~= primary && r1 ~= 0.0
            %%check overlap
            x2 = primaries(i,1);
            y2 = primaries(i,2);
            z2 = primaries(i,3);
            r2 = primaries(i,4);
            %vector along centre-centre line from primary to neighbour
            centre = [x1,y1,z1];
            d = [x2-x1,y2-y1,z2-z1];
            d_2 = (d(1))^2+(d(2))^2+(d(3))^2;    %separation squared
            if  d_2 < (r1+r2)^2 && r2 ~= 0.0
                %truncate sphere
                s=truncate(s,r1,centre,r2,d_2,d);
            end
        end
    end
    
end

%%function to truncate sphere
%%In: primary surface coordinates; primary radius; neighbour radius;
%%separation squared; vector joining primaries; centre of primary
function s_trunc = truncate(s,r1,centre,r2,d_2,d)

    %%calculate size of cap
    x_ij = (d_2 - r2^2 + r1^2) / (2 * sqrt(d_2));
    
    %%rotate coordinates such that line joining primaries lies along the z-axis
    %%axis of rotation perpendicular to z-axis and d
    ax_rot = cross(d,[0,0,1]);  %%cross product of z-axis and d
    %%angle of rotation
    ang_rot = (180/pi)*acos(dot([0,0,1],d)/sqrt(d_2));  %% z.d = |z||d|cos(theta)
    %%rotate s about ax_rot, centred on centre by ang_rot
    rotate(s,ax_rot,ang_rot,centre);
    
    x1 = get(s,'XData');
    y1 = get(s,'YData');
    z1 = get(s,'ZData');
    
    %%truncate cap
    cap = z1-centre(3) > x_ij; 
    z1(cap) = centre(3)+x_ij;

    %%recreate surface
    delete(s);
    s_trunc = surf(x1,y1,z1);
    %%rotate primary back
    rotate(s_trunc,-ax_rot,ang_rot,centre);
    
end