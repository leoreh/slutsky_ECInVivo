function sep = getOSsep()

if ispc
    sep = '\';
else
    sep = '/';
end

return