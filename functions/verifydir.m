% this function verify the existence of the directory specified in 
% "targetpath" and creates it if missing

function verifydir(targetpath)
dirchk=dir(targetpath);
if isempty(dirchk)
       mkdir(targetpath)
       disp('% directory created ...')
else
    disp('% directory already existing ...')
end
end