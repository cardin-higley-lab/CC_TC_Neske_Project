function [trim_chan1, trim_chan4] = tracetrimmer(dataofinterest) 
% find first frame then trim vis_anal to same start as Ca2+ signal
% chan1 is the framerate, chan4 is the analog visual stimulus

chan1= dataofinterest.chan1;
chan4= dataofinterest.chan4;

frametrig = diff(chan1);                                                    % derivates chan1 frame data
xIndex = find (frametrig == max(frametrig), 1);                             % finds time(x) of first frame

trim_chan1 = repmat(chan1,1);
trim_chan1(1:xIndex) = [];

trim_chan4 = repmat(chan4,1);
trim_chan4(1:xIndex) = [];                                                  % deletes chan4 values before Ca2+ trace starts

end

