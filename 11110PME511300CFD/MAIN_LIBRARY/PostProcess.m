function [] = PostProcess(xc, yc, VelocityU, VelocityV, Pressure, VelocityMag)
    %% To save in folder.
    % If not exist, create folder
    current_folder = pwd;
    mkdir pics
    mkdir pics CA3_2022Fall
    % This is the place to put your directory folder
    saved_folder = [current_folder,'\pics\CA3_2022Fall'];
    %% Import globals from 'GlobalsSIMPLE.m'.
    GlobalsSIMPLE;
    % Figure index
    fig = 1;
    %% To plot the figure of velocity quiver.
    figure(fig);
    contourf(xc,yc,VelocityMag, 200,'linecolor','none');
    hold on;
    %% To plot the figure of velocity contour.
    figure(fig);
    quiver(xc, yc, VelocityU', VelocityV', 2, 'w');
    hold off;
    title(sprintf('Velocity Magnitude/Quiver Distribution Using %s', string(sch)));
    subtitle(sprintf('$N_{x} = N_{y} = %d, Re = %2.0e$', nx, Re));
    % To setup the axes
    xlabel('x'); 
    ylabel('y',rotation=0);
    colormap(jet(256));
    colorbar('TickLabelInterpreter', 'latex');
    axis equal;
    set(get(gca,'XLabel'),'FontSize',16);
    set(get(gca,'YLabel'),'FontSize',16);
    set(gca,'FontSize',10);
    % To save the contour
    filename = sprintf('Contour(%s)_N(%d)_Re(%2.0e)_VelocityMagnitude.png', string(sch), nx, Re);
    file     = fullfile(saved_folder, filename);
    exportgraphics(gcf, file);
    %% To plot the figure of pressure contour.
    figure(fig+1);
    contourf(xc,yc,Pressure',20,'linecolor','none');
    title(sprintf('Pressure Magnitude Distribution Using %s', string(sch)));
    subtitle(sprintf('$N_{x} = N_{y} = %d, Re = %2.0e$', nx, Re));
    % To setup the axes
    xlabel('x'); 
    ylabel('y',rotation=0);
    axis equal;
    axis([0 Lx 0 Ly]);
    colormap(jet(256));
    colorbar('TickLabelInterpreter', 'latex');
    set(get(gca,'XLabel'),'FontSize',16);
    set(get(gca,'YLabel'),'FontSize',16);
    set(gca,'FontSize',10);
    % To save the contour
    filename = sprintf('Contour(%s)_N(%d)_Re(%2.0e)_Pressure.png', string(sch), nx, Re);
    file     = fullfile(saved_folder, filename);
    exportgraphics(gcf, file);

    writematrix(VelocityU, sprintf('VelocityU(%s)_N(%d)_Re(%2.0e).csv', string(sch), nx, Re));
    writematrix(VelocityV, sprintf('VelocityV(%s)_N(%d)_Re(%2.0e).csv', string(sch), nx, Re));
end