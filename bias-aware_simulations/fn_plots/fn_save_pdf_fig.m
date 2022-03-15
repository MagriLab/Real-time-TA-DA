
function fn_save_pdf_fig(folder, name)
% folder = ['C:\Users\an553\OneDrive - University of Cambridge\',...
%             'PhD\Papers\Real-time TA DA\Figures_paper\Chaotic\'];
    exportgraphics(gcf,[folder,name,'.pdf'],'ContentType','vector')
end