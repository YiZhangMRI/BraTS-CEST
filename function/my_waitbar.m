% def the class of customized waitbar
classdef my_waitbar
    %%
    properties (SetAccess=private)
        Fig
        Txt
        maxIter
        tstart
        Line1
        Line2
    end
    %% 
    methods
        function WB = my_waitbar(maxIter, FigTitle)
            WB.tstart = tic;
            WB.maxIter = maxIter;
            %% canvas
            WS = get(0,'ScreenSize');
            wt = 500;
            ht = 100;
            WB.Fig = figure('Position',[(WS(3)-wt)/2,(WS(4)-ht)/2,wt,ht],...
                'Name',FigTitle,'NumberTitle','off','menu','none',...
                'Color','white','Resize','off');
            %% icon
            barH = 35;
            AxesIcon = axes(WB.Fig,'Units','pixels','Position',[1,1,ht,ht]);
            try
                icon = imread('./icon.jpg');
            catch
                icon = ind2rgb(round(255*rescale(peaks(100))+1),[1-hot(128); hot(128)]);
            end
            imshow(icon,'Parent',AxesIcon)
            %% text display
            PnlInfo = uipanel(WB.Fig,'Units','pixels','Position', [1+ht,barH,wt-ht,ht-barH]);
            strlist = {'current',0,'total',WB.maxIter,...
                       'past','0','remian','º∆À„÷–'};
            WB.Txt = cell(8,1);
            for i = 0:1
                for j = 0:3
                    n = i*4+j+1;
                    WB.Txt{n} = uicontrol(PnlInfo,'style','edit','Enable','off',...
                        'Units','normalized','Position',[j/4,i/2,1/4,1/2],...
                        'String',strlist{n},'Fontsize',16);
                end
            end
            %% plot waitbar
            AxesBar = axes(WB.Fig,'Units','pixels','Position',[1+ht,1,wt-ht,barH]);
            axis(AxesBar,[-0.05,1.05,-0.2,0.2])
            axis(AxesBar,'off')
            hold(AxesBar,'on')
            WB.Line1 = plot(AxesBar,[0,1],[0,0],'-','LineWidth',15,'Color',[0.9,0.9,0.9]);
            WB.Line2 = plot(AxesBar,[0,0],[0,0],'-','LineWidth',15,'Color',[0.1,0.9,0.1]);            
            hold(AxesBar,'off')
            drawnow

        end
        
        function updata(WB,iter) 
            %% num iter
            set(WB.Txt{2},'string',iter)
            %% time
            tnow = toc(WB.tstart);
            trem = tnow/iter*(WB.maxIter-iter);
            set(WB.Txt{6},'string',...
                [int2str(floor(tnow/60)),':',int2str(floor(rem(tnow,60)))])
            set(WB.Txt{8},'string',...
                [int2str(floor(trem/60)),':',int2str(floor(rem(trem,60)))])
            %% length for waitbar
            WB.Line1.XData = [iter/WB.maxIter,1];
            WB.Line2.XData = [0,iter/WB.maxIter];
            drawnow

        end

        function closeWaitBar(WB)
            close(WB.Fig)
            clear("WB")
        end
        
    end

end
















