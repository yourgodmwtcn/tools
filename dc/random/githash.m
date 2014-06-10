function [hash] = githash()

    [~, hash] = system('TERM=xterm-256color git log -n 1 --pretty=format:''%H''');
    % remove bash escape characters
    hash = hash(9:48)