% Constrained fitting
function GParameters = Fit4LorentzianCplxPh_Dixon_Con_20230126(Freq, DataCplx, StartPoint, LowerBounds, UpperBounds)

warning('off', 'all')

% Ensure sufficient amplitude for accurate fitting. Otherwise fminunc will terminate fitting process too early
ScaleFac = 1000 / max(abs(DataCplx));
% ScaleFac = 1;
DataCplx_Scaled = ScaleFac * DataCplx;
StartPoint(1) = ScaleFac * StartPoint(1);
StartPoint(5) = ScaleFac * StartPoint(5);
StartPoint(9) = ScaleFac * StartPoint(9);
StartPoint(13) = ScaleFac * StartPoint(13);
StartPoint(17) = ScaleFac * StartPoint(17);
StartPoint(18) = ScaleFac * StartPoint(18);
LowerBounds(1) = ScaleFac * LowerBounds(1);
LowerBounds(5) = ScaleFac * LowerBounds(5);
LowerBounds(9) = ScaleFac * LowerBounds(9);
LowerBounds(13) = ScaleFac * LowerBounds(13);
LowerBounds(17) = ScaleFac * LowerBounds(17);
LowerBounds(18) = ScaleFac * LowerBounds(18);
UpperBounds(1) = ScaleFac * UpperBounds(1);
UpperBounds(5) = ScaleFac * UpperBounds(5);
UpperBounds(9) = ScaleFac * UpperBounds(9);
UpperBounds(13) = ScaleFac * UpperBounds(13);
UpperBounds(17) = ScaleFac * UpperBounds(17);
UpperBounds(18) = ScaleFac * UpperBounds(18);

% Use fitting parameters from an initial fit with fminsearchbnd as starting parameters for fmincon for improved fitting
FMS_Parameters = fminsearchbnd(@lor4fun, StartPoint, LowerBounds, UpperBounds, optimset('MaxFunEvals',10000,'MaxIter',10000));
options = optimoptions(@fmincon,'MaxFunctionEvaluations',10000,'Display','off');
GParameters = fmincon(@lor4fun, FMS_Parameters, [], [], [], [], LowerBounds, UpperBounds, [], options);
% GParameters = fmincon(@lor4fun, StartPoint, [], [], [], [], LowerBounds, UpperBounds, [], options);
% options = optimoptions(@fminunc,'MaxFunctionEvaluations',10000,'Display','off');
% GParameters = fminunc(@lor4fun, StartPoint);
% GParameters = fminunc(@lor4fun, FMS_Parameters);

% Undo scaling of fitting parameters
GParameters(1) = GParameters(1)/ScaleFac;
GParameters(5) = GParameters(5)/ScaleFac;
GParameters(9) = GParameters(9)/ScaleFac;
GParameters(13) = GParameters(13)/ScaleFac;
GParameters(17) = GParameters(17)/ScaleFac;
GParameters(18) = GParameters(18)/ScaleFac;

warning('on', 'all')

    function sse = lor4fun(params)
        Amp1        = params(1);
        HalfWidthL1 = params(2);
        Shift1      = params(3);
        Phase1      = params(4);
        Amp2        = params(5);
        HalfWidthL2 = params(6);
        Shift2      = params(7);
        Phase2      = params(8);
        Amp3        = params(9);
        HalfWidthL3 = params(10);
        Shift3      = params(11);
        Phase3      = params(12);
        Amp4        = params(13);
        HalfWidthL4 = params(14);
        Shift4      = params(15);
        Phase4      = params(16);
        Offset_R    = params(17);
        Offset_I    = params(18);
        Offset  = complex(Offset_R, Offset_I);

        FitPeakLorR1 = HalfWidthL1 ./ ((Freq - Shift1) .* (Freq - Shift1) + HalfWidthL1 * HalfWidthL1);
        % FitPeakLorI1 = -(Freq - Shift1) ./ ((Freq - Shift1) .* (Freq - Shift1) + HalfWidthL1 * HalfWidthL1);
        FitPeakLorI1 = (Freq - Shift1) ./ ((Freq - Shift1) .* (Freq - Shift1) + HalfWidthL1 * HalfWidthL1);
        FitPeakC1    = Amp1 * complex(FitPeakLorR1, FitPeakLorI1) * exp(complex(0,Phase1));
        FitPeakLorR2 = HalfWidthL2 ./ ((Freq - Shift2) .* (Freq - Shift2) + HalfWidthL2 * HalfWidthL2);
        % FitPeakLorI2 = -(Freq - Shift2) ./ ((Freq - Shift2) .* (Freq - Shift2) + HalfWidthL2 * HalfWidthL2);
        FitPeakLorI2 = (Freq - Shift2) ./ ((Freq - Shift2) .* (Freq - Shift2) + HalfWidthL2 * HalfWidthL2);
        FitPeakC2    = Amp2 * complex(FitPeakLorR2, FitPeakLorI2) * exp(complex(0,Phase2));
        FitPeakLorR3 = HalfWidthL3 ./ ((Freq - Shift3) .* (Freq - Shift3) + HalfWidthL3 * HalfWidthL3);
        % FitPeakLorI3 = -(Freq - Shift3) ./ ((Freq - Shift3) .* (Freq - Shift3) + HalfWidthL3 * HalfWidthL3);
        FitPeakLorI3 = (Freq - Shift3) ./ ((Freq - Shift3) .* (Freq - Shift3) + HalfWidthL3 * HalfWidthL3);
        FitPeakC3    = Amp3 * complex(FitPeakLorR3, FitPeakLorI3) * exp(complex(0,Phase3));
        FitPeakLorR4 = HalfWidthL4 ./ ((Freq - Shift4) .* (Freq - Shift4) + HalfWidthL4 * HalfWidthL4);
        % FitPeakLorI4 = -(Freq - Shift4) ./ ((Freq - Shift4) .* (Freq - Shift4) + HalfWidthL4 * HalfWidthL4);
        FitPeakLorI4 = (Freq - Shift4) ./ ((Freq - Shift4) .* (Freq - Shift4) + HalfWidthL4 * HalfWidthL4);
        FitPeakC4    = Amp4 * complex(FitPeakLorR4, FitPeakLorI4) * exp(complex(0,Phase4));
        FittedCurve = FitPeakC1 + FitPeakC2 + FitPeakC3 + FitPeakC4 + Offset;
        ErrorVector = FittedCurve - DataCplx_Scaled;
        sse = double(sum(abs(ErrorVector) .^ 2));
    end
end
