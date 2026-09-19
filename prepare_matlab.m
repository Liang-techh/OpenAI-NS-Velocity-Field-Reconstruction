function prepare_matlab
%PREPARE_MATLAB Restore missing pinned MATLAB sources/data, once per checkout.
% Complete offline bundles already contain these files. Downloads are public,
% SHA-256 checked, and restricted to the declared MATLAB destination folder.
% Existing files with different bytes are never silently overwritten.
root=fileparts(mfilename('fullpath'));
M=jsondecode(fileread(fullfile(root,'tools','source_manifest.json')));
entries=M.entries;
for k=1:numel(entries)
    e=entries(k);rel=char(e.path);
    if ~startsWith(rel,'visualization/matlab/'),continue,end
    if contains(rel,'..')||contains(rel,'\')||~strcmp(e.source_repository,'Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1')
        error('ns:ImportPath','Unexpected pinned MATLAB import');
    end
    dest=fullfile(root,strrep(rel,'/',filesep));
    if isfile(dest)
        if ~strcmp(filehash(dest),char(e.sha256)),error('ns:ModifiedFile','Existing file differs from the pinned source: %s. Preserve your edits before replacing it.',rel);end
        continue
    end
    folder=fileparts(dest);if ~isfolder(folder),mkdir(folder);end
    tmp=tempname(folder);cleanup=onCleanup(@()removeTemp(tmp));
    url=['https://raw.githubusercontent.com/' char(e.source_repository) '/' char(e.source_commit) '/' char(e.source_path)];
    fprintf('Restoring pinned MATLAB asset: %s\n',rel);
    websave(tmp,url,weboptions('Timeout',90,'ContentType','binary'));
    if ~strcmp(filehash(tmp),char(e.sha256)),error('ns:Checksum','Downloaded bytes failed SHA-256 verification: %s',rel);end
    movefile(tmp,dest);clear cleanup
end
end
function value=filehash(path)
fid=fopen(path,'rb');if fid<0,error('ns:File','Cannot read %s',path);end
c=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');clear c
md=java.security.MessageDigest.getInstance('SHA-256');md.update(typecast(bytes(:),'int8'));
value=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
function removeTemp(path)
if isfile(path),delete(path);end
end
