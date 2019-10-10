%{  
 Author:    Casper Lindberg (csl37)

 Purpose:   Functions to draw a single particle

 In:        List of primary particle coordinates, 
            facealpha: draw solid or TEM style (optional)
 Out:       Figure and axes handles
%}
function [fig, ax] = draw_particle(coords,solid)
       
    %check for solid (optional)
    %in not supplied reorientate all coordinates
    if (~exist('solid', 'var'))
        solid = 0.5;    %assume TEM style, set alpha to 0.5
    end

    figure('Position', [ 0 0 1000 1000],'visible','off')
    hold on
    
    [L,m,n] = size(coords);
    for i = 1:L
        %check if primary (radius > 0)
        if (coords(i,4) > 0 )
            %%loop over primaries
            s = draw_primary(coords(i,:));
            %%check for overlap and truncate sphere
            s = check_overlap(s,i,coords(:,:));
            %draw sphere
            %set(s,'FaceColor','interp','FaceAlpha',0.65,'EdgeColor' ,'none'); 
            set(s,'FaceColor',[0.7,0.7,0.7],'FaceAlpha',solid,'EdgeColor' ,'none'); 
        end
    end

    %% pass figure and axes handles
    ax= gca;
    fig = gcf;
    
end

%% Draws a spherical primary particle
function s = draw_primary(coords)

    [x0,y0,z0] = sphere(50); %coordinates of 50x50 face unit sphere at (0,0,0) %sphere(N) specifies N faces

    %rescale and recentre sphere 
    x1 = x0*coords(4) + coords(1);
    y1 = y0*coords(4) + coords(2);
    z1 = z0*coords(4) + coords(3);

    s = surf(x1,y1,z1);

end

%% Check for overlaps between primaries and truncate sphere
%% In:  surface coordinates of primary being drawn 
%%      index of primary being drawn
%%      list of primary coordinates
function s=check_overlap(s,primary_index,primaries)

    [L,~] = size(primaries);
    %%loop through list of primaries checking for overlap
   
    x1 = primaries(primary_index,1);
    y1 = primaries(primary_index,2);
    z1 = primaries(primary_index,3);
    r1 = primaries(primary_index,4);
    
    for i = 1:L
        if i ~= primary_index && r1 ~= 0.0
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
                %  Not absolutely necessary but may improve a TEM style image
                s=truncate(s,r1,centre,r2,d_2,d);
            end
        end
    end
    
end

%%  function to truncate sphere
%%  Not absolutely necessary but may improve a TEM style image
%%
%% In:  s:      surface coordinates of primary being drawn 
%%      r1:     radius of primary being drawn
%%      centre: coordinates of primary being drawn
%%      r2:     neighbour separation squared
%%      d:      x_neighbour - x_primary
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