clear
for l=1:40
    %micelle
%     PNIPAM_chain
%     axes.SortMethod = 'ChildOrder';
%     addpath(genpath('D:\Documents\Project\Matlab\export_fig-master\export_fig-master'));
%     export_fig(['PNIPAM_rb',num2str(l),'.png'],'-r300','-a2','-opengl','-transparent')
    
    %unimer
    PNIPAM_chain_unimer
    axes.SortMethod = 'ChildOrder';
    addpath(genpath('D:\Documents\Project\Matlab\export_fig-master\export_fig-master'));
    export_fig(['PNIPAM_rb_u',num2str(l),'.png'],'-r300','-a2','-opengl','-transparent')
    
    clear
end