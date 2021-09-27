clear
for l=1:20
    
    WLM
    axes.SortMethod = 'ChildOrder';
    addpath(genpath('X:\Documents\Project\Matlab\export_fig-master\export_fig-master'));
    export_fig(['./figures/WLM_a2_DP100_',num2str(l),'.png'],'-r192','-a2','-opengl','-transparent','-nocrop')
    
    clear
end