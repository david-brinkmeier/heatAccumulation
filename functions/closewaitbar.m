function [] = closewaitbar()
progressbar = findall(0,'type','figure','tag','TMWWaitbar');
if ~isempty(progressbar)
    delete(progressbar);
end
end

