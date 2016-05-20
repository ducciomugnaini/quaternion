classdef coosReader
    
    methods (Static)
        
        function XYZ = Coo3DReader(filePath)            
            
            fileID = fopen(filePath);
            
            idCurv = 1;
            
            XYZ{1} = [];
            
            % disp([ '>> reading curve #' num2str(idCurv)]);
            while (~feof(fileID))
                
                currLine = textscan(fileID,'%s',1,'Delimiter','\n');
                currRow = char(currLine{1});
                splittedRow = strsplit(currRow,' ');
                
                numCoo = size(splittedRow,2);
                
                switch(numCoo)
                    
                    case 2
                        
                        XYZ{idCurv} = [XYZ{idCurv}; str2num(splittedRow{1}) str2num(splittedRow{2})];
                        
                    case 3
                        
                        XYZ{idCurv} = [XYZ{idCurv}; str2num(splittedRow{1}) str2num(splittedRow{2}) str2num(splittedRow{3})];
                    
                    otherwise
                        idCurv = idCurv + 1;
                        XYZ{idCurv} = [];
                end
                
                % disp([ '>> reading curve #' num2str(idCurv)]);
            end 
            % disp([ '>> curve #' num2str(idCurv) ' is empty']);
            
        end%_K3DReader
        
        function plotRef(arrowLength, idFig)
            figure(idFig);
            quiver3(0,0,0, arrowLength,0,0, 'Color','k');
            quiver3(0,0,0, 0,arrowLength,0, 'Color','k');
            quiver3(0,0,0, 0,0,arrowLength, 'Color','k');
        end
        
    end
    
end

