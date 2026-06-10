function [y, areas] = evalVoigt(params,freq)

% initialize outputs
y = zeros(size(freq));
areas = zeros(size(params,1),1);

% loop through each fit
for i = 1:size(params,1)
    amp = params(i,1);
    shift = params(i,2);
    fwhmL = params(i,3);
    fwhmG = params(i,3);
    eta = params(i,4);
    phase = params(i,5);
    phaserec = params(i,6);
    offsetR = params(i,7);
    offsetI = params(i,8);

    df = freq - shift;
    offset = complex(offsetR,offsetI);

    lorR = fwhmL ./ (df.^2 + fwhmL.^2);
    lorI = -df ./ (df.^2 + fwhmL.^2);
    gauss = exp(-(df .* sqrt(log(2)) / fwhmG).^2);

    v = amp .* exp(1i*phase) .* ...
        (eta * (complex(lorR,0) + complex(0,lorI) * exp(1i*phaserec)) +...
        (1 - eta) * gauss) + offset;

    areas(i) = sum(abs(v))./length(v);

    if size(params,2) == 9
        areas(i) = areas(i) ./ params(i,9);
        v = v ./ params(i,9);
    end
    
    y = y + v;
end
end