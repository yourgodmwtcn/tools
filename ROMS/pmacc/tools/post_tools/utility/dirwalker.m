function [netcdfpaths,folderarray]=dirwalker(folderarray)
%%[netcdfpaths]=dirwalker(folderarray) returns the full netcdf file paths
%%of all netcdf files within the all of the folders in folderarray
%The argument folderarray should be a cell array of one or more folder
%names or a string of a single folder name
%To verify which folders and subfolders were searched, return the second
%argument 'folderarray' which is the original input argument of folder
%names, with any subfolders found appended to the end of the array
if ~isa(folderarray,'cell')
    folderarray={folderarray};
end

slash = '/'; % '\' on a PC

%now add on any .nc files to the output that were in the original folderarray directly
%specified- 
count=0;
for i=1:length(folderarray)
    h=dir(char(folderarray(i)));
    if length(h)==1&&h.bytes>0
        count=count+1;
        netcdfpaths{count}=[char(folderarray(i)) ];
    end
    
end

fol.isdir=[1 1 1];
start=1;
i=0;
while i<length(folderarray)
    for i=start:length(folderarray)
%        if ~strcmp(folderarray{length(folderarray)}(end),'\')
        if ~strcmp(folderarray{i}(end),slash)
            folderarray{i}(end+1)=slash; %changed from length(folderarray)
        end

        h=dir(char(folderarray(i)));
        p=struct2cell(h);
        fol.isdir=cell2mat(p(4,:));
        if length(h)>2
            if any(fol.isdir(3:end))%if there are any subdirectories then append them onto the list of folderarray
                inds=find(fol.isdir(3:end));
                for j=1:length(inds)
                    
                    folderarray{length(folderarray)+1}=[char(folderarray(i)) h(inds(j)+2).name];
                    
                end
            end
         

        end
      
           
            nccell{i}=dir([char(folderarray(i)) '*.nc']);


    end
    start=i+1;
end



for i=1:length(folderarray)
    for j=1:length(nccell{i})
        count=count+1;
        netcdfpaths{count}=[char(folderarray(i)) char(nccell{i}(j).name)];
    end
end

if ~count %if no netcdf files were found
    netcdfpaths='No Files Found';
    
end



