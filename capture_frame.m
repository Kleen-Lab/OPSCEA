function frame = capture_frame(fig, reset)
% CAPTURE_FRAME  getframe wrapper that enforces a consistent cdata size
% across an entire video capture session.
%
% VideoWriter requires every frame's 'cdata' to be the same size. A stray
% display sleep, window move/resize, or DPI change mid-run can make a
% single getframe() call return a differently-sized cdata, which
% otherwise only surfaces as a crash at writeVideo() after the full
% render is done. This wrapper remembers the size of the first captured
% frame and resizes any later mismatch back to it, warning instead of
% failing.
%
% Call capture_frame(fig, true) once to reset/start a new session before
% the capture loop begins.

persistent targetSize
if nargin > 1 && reset
    targetSize = [];
end

frame = getframe(fig);
sz = size(frame.cdata);

if isempty(targetSize)
    targetSize = sz;
elseif ~isequal(sz, targetSize)
    warning('OPSCEA:frameSizeMismatch', ...
        'Captured frame size %s does not match session target %s (likely display sleep/resize) - resizing to match.', ...
        mat2str(sz), mat2str(targetSize));
    frame.cdata = imresize(frame.cdata, targetSize(1:2));
end
end
