function str=datetimefilename
%DATETIMEFILENAME Date time filename. Generate a string suitable for
%filename based on the date/time.
clockvec=clock;
str=[num2str(clockvec(1)),'-',num2str(clockvec(2)),'-',num2str(clockvec(3)),...
    '-',num2str(clockvec(4)),'-',num2str(clockvec(5)),'-',num2str(clockvec(6))];
return