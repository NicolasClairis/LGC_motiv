
      
a = 0;
timelimit = false;
t_max = 10;
keySpace = 65;

onset = GetSecs;
timeNow = onset;

KbReleaseWait();
while (a < 5) &&...
        ( ( (timelimit == true) && (timeNow < onset + t_max) ) || (timelimit == false) )
    
    timeNow = GetSecs;
    
    [keyisdown, timeAnswer, keyCode] = KbCheck();
    if keyisdown == 1 && keyCode(keySpace) == 1
        a = a + 1;
        disp(num2str(a));
    end
    %     WaitSecs(0.001);
    KbReleaseWait;
end % loop