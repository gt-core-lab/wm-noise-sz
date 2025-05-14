function [result, endSecs, keyCode] = get_response(answer,keyboardnum)
FlushEvents('keyDown');

firstKey = KbName('LeftArrow');
secondKey = KbName('RightArrow');

t1 = GetSecs;
keyState = 0;
while keyState == 0
	[keyIsDown, endSecs, keyCode] = KbCheck(keyboardnum);
    if answer == -1
        if keyCode(firstKey)
            result = 1;
            keyState = 1;
            break;
        elseif keyCode(secondKey)
            result = 0;
            keyState = 1;
            break;
        elseif keyCode(KbName('q'))
            Screen('CloseAll');
            ShowCursor;
            break;
        end
    elseif answer == 1
        if keyCode(firstKey)
            result = 0;
            keyState = 1;
            break;
        elseif keyCode(secondKey)
            result = 1;
            keyState = 1;
            break;
        elseif keyCode(KbName('q'))
            Screen('CloseAll');
            ShowCursor;
            break;
        end
    end
end
         

end

