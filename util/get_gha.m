function gha = get_gha(time)
    secs_in_day = 86400;
    time_of_day = mod(time, secs_in_day); % seconds in current day
    
    % calculate the greenwich hour angle
    % positive means the prime meridian is east of the vernal equinox
    gha = time_of_day/secs_in_day*360; % degrees
end