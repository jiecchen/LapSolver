function ijvwrite(filename, a)
%IJVWRITE Exports a graph in ijv adjacency list form.
    [ai,aj,av] = find(tril(a,-1));
    dlmwrite(filename, [length(a), length(ai)], ' ');
    dlmwrite(filename, [ai-1 aj-1 av], 'delimiter', ' ', '-append');
end

