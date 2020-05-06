function quantilePlot(Time, Quantiles, baseColor)

if nargin<3
       baseColor = [44,127,184]./255;  % Blue
end

if size(Quantiles, 2) == 5
    OuterBot = Quantiles(:, 1);
    InnerBot = Quantiles(:, 2);
    Center   = Quantiles(:, 3);
    InnerTop = Quantiles(:, 4);
    OuterTop = Quantiles(:, 5);

    plot(Time, Center, 'LineWidth', 2, 'Color', baseColor)
    hold on

    InnerIdx = ~isnan(InnerBot + InnerTop);
    fill([    Time(InnerIdx)    ; flip(Time(InnerIdx))], ...
         [InnerBot(InnerIdx); flip(InnerTop(InnerIdx))], baseColor,...
        'LineStyle','none',...
        'FaceAlpha', .4);


    OuterIdx = ~isnan(OuterBot + OuterTop);
    fill([   Time(OuterIdx); flip(    Time(OuterIdx))], ...
        [OuterBot(OuterIdx); flip(OuterTop(InnerIdx))], baseColor,...
        'LineStyle','none',...
        'FaceAlpha', .15);

elseif size(Quantiles, 2) == 7
    OuterBot  = Quantiles(:, 1);
    MiddleBot = Quantiles(:, 2);
    InnerBot  = Quantiles(:, 3);
    Center    = Quantiles(:, 4);
    InnerTop  = Quantiles(:, 5);
    MiddleTop =  Quantiles(:, 6);
    OuterTop  =  Quantiles(:, 7);
    baseColor = [44,127,184]./255;

    plot(Time, Center, 'LineWidth', 2, 'Color', baseColor * 1.2)
    hold on

    InnerIdx = ~isnan(InnerBot + InnerTop);
    fill([   Time(InnerIdx); flip(    Time(InnerIdx))], ...
        [InnerBot(InnerIdx); flip(InnerTop(InnerIdx))], baseColor,...
        'LineStyle','none',...
        'FaceAlpha', .5);

    MiddleIdx = ~isnan(MiddleBot + MiddleTop);
    fill([    Time(MiddleIdx); flip(     Time(MiddleIdx))], ...
        [MiddleBot(MiddleIdx); flip(MiddleTop(InnerIdx))], baseColor,...
        'LineStyle','none',...
        'FaceAlpha', .25);



    OuterIdx = ~isnan(OuterBot + OuterTop);
    fill([   Time(OuterIdx); flip(    Time(OuterIdx))], ...
        [OuterBot(OuterIdx); flip(OuterTop(InnerIdx))], baseColor,...
        'LineStyle','none',...
        'FaceAlpha', .15);

else
    error('Enter a valid quantile matrix')

end

end
