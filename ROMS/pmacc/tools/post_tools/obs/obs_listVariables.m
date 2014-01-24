function [FilesUsed,AllVars,VarMatrix]=obs_listVariables(files, dim)

% [FilesUsed,AllVars,VarMatrix]=obs_listVariables(files, dim)
%
%takes the file input (used for obs_extract) and creates a figure and
%text file showing all variables listed in these files.  This program does
%not differentiate between a variable existing and it being full of nans.
%
%dim is the dimension of variables in the netcdf files
%mooring files have dim=2, all other files have dim=1.  If dim is not assigned then
%variables from both types of files are returned.
% built to work with standard netcdf files created for obs_extract package

%C. Bassin Dec 2010


if nargin==1
    dim=0;
end




[netcdfpaths,folderarray]=dirwalker(files);

if strcmp(netcdfpaths,'No Files Found')
    disp(['no netcdf files found in ' files])
    FilesUsed=nan;
    AllVars=nan;
    VarMatrix=nan;
    
    return
end

    
x=1;
for i=1:length(netcdfpaths);

    Info_file=nc_info(netcdfpaths{i});
         switch dim
             case 0
                  if size(Info_file.Dimension,1)==1 | size(Info_file.Dimension,1)==2
                     V{x,1}={Info_file.Dataset.Name};
                     V{x,2}=size(Info_file.Dataset,1);
                     V{x,3}=netcdfpaths{i};
                     x=x+1;

                  else 
                     disp(['skipping ' netcdfpaths{i} ' due to incorrect dimension size'])
                  end
             case {1,2}
                if size(Info_file.Dimension,1)==dim
                     V{x,1}={Info_file.Dataset.Name};
                     V{x,2}=size(Info_file.Dataset,1);
                     V{x,3}=netcdfpaths{i};
                     x=x+1;

                else 
                disp(['skipping ' netcdfpaths{i} ' due to incorrect dimension size'])
                end
         end

       
  

end

     if exist('V','var')==0
     disp(['no netcdf files found in with dimension size = ' num2str(dim)])
    FilesUsed=nan;
    AllVars=nan;
    VarMatrix=nan;
     return
     end
    
    AllVars=unique([V{:,1}]);
    p=size(AllVars,2);
    VarMatrix=nan(x-1,p);
    for y=1:x-1;
        for z=1:p;
            Findvar=strcmp(AllVars{z},V{y,1});
            if sum(Findvar)>0
               % FV=find(Findvar==1);
             VarMatrix(y,z)=1;
            end
           
        end
 
    scatter(1:p ,ones(p,1)*y,100,VarMatrix(y,:),'s','filled')
    hold on
    
    end
    
    set(gca,'xticklabel',' ')
    ax = axis; % Current axis limits
    axis(axis); % Set the axis limit modes (e.g. XLimMode) to manual
    Yl = ax(3:4);
    set(gca,'xtick',1:p,'FontSize',8)


    t = text(1:p,Yl(1)*ones(1,p),AllVars);
    set(t,'HorizontalAlignment','right','VerticalAlignment','top', ...
    'Rotation',45);

%truncate file names
FilesUsed=V(:,3);
 a=regexp(FilesUsed,'\','split');
 for l=1:size(a,1)
     F{l,1}=a{l}{end};
 end
     set(gca,'ytick',1:x-1 ,'yticklabel',F(:,1) ,'FontSize',8)
if x > 2, ylim([1 x-1]); end
    xlim([1 p])    
    
    


% outputfilename=['Vars_in_nc_files_with_dim1' num2str(dim) '.txt'];
% 
% fid = fopen(outputfilename, 'wt');
% 
% fprintf(fid,'Filename \t');
% fprintf(fid,'\t');
% for i=1:p; fprintf(fid,[AllVars{i} '\t']); end
% fprintf(fid,'\n');
% 
%     for j=1:x-1;
% 
%     a=regexp(FilesUsed{j},'\');
%         fprintf(fid,[FilesUsed{j}(a(end)+1:end-3) '\t \t']);
%          for q=1:p-1;
%         fprintf(fid,'%d\t',VarMatrix(j,q));
%          end
%       fprintf(fid,'%d\n',VarMatrix(j,end));
%     end
% 
% fclose(fid)



