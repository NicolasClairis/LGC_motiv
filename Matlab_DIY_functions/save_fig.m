function [] = save_fig( fig_hdl, dirPath, img_nm, overwrite_allowed )
%[] = save_fig( fig_hdl, dirPath, img_nm, overwrite_allowed )
% save_fig will save the figure stored in the figure handle fig_hdl inside
% the dirPath directory with the img_nm name.
%
% INPUTS
% fig_hdl: figure handle (what you get if you write fig_hdl = figure;)
%
% dirPath: path of the directory where you want to store your figure
%
% img_nm: name you want to give to your new figure
%
% overwrite_allowed:
% (0) if a file already exists with the specified name do not overwrite it
% (1) overwrite any file with the same name
%

%% set default parameters
if ~exist('fig_hdl','var') || isempty(fig_hdl)
    fig_hdl = gcf; % by default current figure (= last opened)
end
if ~exist(dirPath,'dir') || isempty(dirPath)
    dirPath = pwd; % by default current folder if no path is specified
end
if ~exist('img_nm','var') || isempty(img_nm)
    img_nm = 'figure.png';
end
if ~exist('overwrite_allowed','var') || isempty(overwrite_allowed)
    overwrite_allowed = 0; % do not allow overwritting by default
end

%% save the image
image_full_path = [dirPath, filesep, img_nm];
if ~exist([image_full_path,'.png'],'file') ||...
        ~exist(image_full_path,'file') ||...
        overwrite_allowed == 1
    
    % maximize figure size
    set(fig_hdl,'PaperPosition',[0 0 1 1]);
    set(fig_hdl,'PaperPositionMode','auto');
    
    if strcmp(image_full_path(end-3:end),'.png')
        % if name of the file already includes the extension of the file,
        % save it as it is
        saveas(fig_hdl, image_full_path);
    else % otherwise add '.png' at the end of the name
        saveas(fig_hdl, [image_full_path,'.png']);
    end
    
else % if overwritting not allowed
    save_nevertheless_yn = questdlg(['A figure with the name ',img_nm,' already exists in the specified location',...
        ' Do you want to overwrite it?']);
    if strcmp(save_nevertheless_yn,'Yes')
        set(fig_hdl,'PaperPosition',[0 0 1 1]);
        set(fig_hdl,'PaperPositionMode','auto');
        
        if strcmp(image_full_path(end-3:end),'.png')
            % if name of the file already includes the extension of the file,
            % save it as it is
            saveas(fig_hdl, image_full_path);
        else % otherwise add '.png' at the end of the name
            saveas(fig_hdl, [image_full_path,'.png']);
        end
        
    else
        save_but_with_different_name = questdlg('Are you willing to save it with a _bis extension at the end?');
        if strcmp(save_but_with_different_name,'Yes')
            while exist(image_full_path,'file')
                image_full_path = [image_full_path, '_bis'];
            end
            set(fig_hdl,'PaperPosition',[0 0 1 1]);
            set(fig_hdl,'PaperPositionMode','auto');
            
            
            if strcmp(image_full_path(end-3:end),'.png')
                % if name of the file already includes the extension of the file,
                % save it as it is
                saveas(fig_hdl, image_full_path);
            else % otherwise add '.png' at the end of the name
                saveas(fig_hdl, [image_full_path,'.png']);
            end
            
        end
    end
end

end % function