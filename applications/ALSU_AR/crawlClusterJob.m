function allResults = crawlClusterJob(pathToJobID,saveFlag)
% Collects the data from the individual 'results.mat' files from a cluster run into one structure

	files = dir(pathToJobID);
	% Remove current and main
	dirFlags = [files.isdir] & ~strcmp({files.name},'.') & ~strcmp({files.name},'..');
	
	subFolders = files(dirFlags);
	
	if isempty(subFolders)
		error('Nothing found.')
	end
	
	% Print folder names to command window.
	i=1;
	for k = 1 : length(subFolders)
		fprintf('Sub folder #%d = %s\n', k, subFolders(k).name);
		
		if exist(fullfile(subFolders(k).folder,subFolders(k).name,'results.mat'))
			
			locResults = load(fullfile(subFolders(k).folder,subFolders(k).name,'results.mat'));
			
			fieldNames = fieldnames(locResults.results);
			fieldValues = struct2cell(locResults.results);
			
			for nF = 1:length(fieldNames)
				allResults.(fieldNames{nF}){i} = fieldValues{nF};
			end
			
			allResults.RUNID(i) = str2double(subFolders(k).name);
			
			i = i+1;
		else
			fprintf('WARNING: no results found in folder ''%s'' \n',subFolders(k).name)
		end
	end
	if saveFlag
		save(['results' pathToJobID(end-7:end)],'allResults','-v7.3')
	end
end
