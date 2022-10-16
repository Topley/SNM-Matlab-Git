function saveVars4MatFile(Vars)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
if contains(fullFilename, 'TA')
    
    TACOHT = m1COHT;
    TAF = m1F;
    TACOF = m1COF;
    TACOF2 = m1COF2;
    TAZ = m1Z;
    TACOFz = m1COFz;
    TACOF2z = m1COF2z;
    
    MGCOHT = m2COHT;
    MGF = m2F;
    MGCOF = m2COF;
    MGCOF2 = m2COF2;
    MGZ = m2Z;
    MGCOFz = m2COFz;
    MGCOF2z = m2COF2z;
    
    SolCOHT = m3COHT;
    SolF = m3F;
    SolCOF = m3COF;
    SolCOF2 = m3COF2;
    SolGZ = m3Z;
    SolCOFz = m3COFz;
    SolCOF2z = m3COF2z;
    
%     LGCOHT = m4COHT;
%     LGF = m4F;
%     LGCOF = m4COF;
%     LGCOF2 = m4COF2;
%     LGZ = m4Z;
%     LGCOFz = m4COFz;
%     LGCOF2z = m4COF2z;
    
    try
        filename_save = [fullFilename(1:end-4) '_pooled_coher_', num2str(stax), '_',num2str(endax)];
        save([filename_save '.mat'],'TACOHT','TAF','TACOF','TACOF2','TAZ','TACOFz','TACOF2z',...
            'MGCOHT','MGF','MGCOF','MGCOF2','MGZ','MGCOFz','MGCOF2z',...
            'SolCOHT','SolF','SolCOF','SolCOF2','SolZ','SolCOFz','SolCOF2z',...
            'cstCOHT','cstF','cstCOF','cstCOF2','cstZ','cstCOFz','cstCOF2z',...
            'poolCOHT','poolF','poolCOF','poolCOF2','poolZ','poolCOFz','poolCOF2z',...
            'fsamp','MUNumber','MULength', 'WinLen','rSpike1','a','meanDR','CoVISI');
        print (gcf,'-dpdf',[filename_save '.pdf']);
        print(figure(99),'-dpdf',[filename_save '_comp.pdf']);
        
    catch
        try
            filename_save = [fullFilename(1:end-4) '_pooled_coher_', num2str(stax), '_',num2str(endax)];
            save([filename_save '.mat'],'TACOHT','TAF','TACOF','TACOF2','TAZ','TACOFz','TACOF2z',...
                'MGCOHT','MGF','MGCOF','MGCOF2','MGZ','MGCOFz','MGCOF2z',...
                'cstCOHT','cstF','cstCOF','cstCOF2','cstZ','cstCOFz','cstCOF2z',...
                'poolCOHT','poolF','poolCOF','poolCOF2','poolZ','poolCOFz','poolCOF2z',...
                'fsamp','MUNumber','MULength', 'WinLen','rSpike1','a','meanDR','CoVISI');
            print (gcf,'-dpdf',[filename_save '.pdf']);
            print(figure(99),'-dpdf',[filename_save '_comp.pdf']);
        catch
            filename_save = [fullFilename(1:end-4) '_pooled_coher_', num2str(stax), '_',num2str(endax)];
            save([filename_save '.mat'],'TACOHT','TAF','TACOF','TACOF2','TAZ','TACOFz','TACOF2z',...
                'cstCOHT','cstF','cstCOF','cstCOF2','cstZ','cstCOFz','cstCOF2z',...
                'poolCOHT','poolF','poolCOF','poolCOF2','poolZ','poolCOFz','poolCOF2z',...
                'fsamp','MUNumber','MULength', 'WinLen','rSpike1','a','meanDR','CoVISI');
            print (gcf,'-dpdf',[filename_save '.pdf']);
            print(figure(99),'-dpdf',[filename_save '_comp.pdf']);
        end
    end
    
elseif contains(fullFilename, 'MG')
    
    TACOHT = m2COHT;
    TAF = m2F;
    TACOF = m2COF;
    TACOF2 = m2COF2;
    TAZ = m2Z;
    TACOFz = m2COFz;
    TACOF2z = m2COF2z;
    
    MGCOHT = m1COHT;
    MGF = m1F;
    MGCOF = m1COF;
    MGCOF2 = m1COF2;
    MGZ = m1Z;
    MGCOFz = m1COFz;
    MGCOF2z = m1COF2z;
    
    SolCOHT = m3COHT;
    SolF = m3F;
    SolCOF = m3COF;
    SolCOF2 = m3COF2;
    SolGZ = m3Z;
    SolCOFz = m3COFz;
    SolCOF2z = m3COF2z;
    
    LGCOHT = m4COHT;
    LGF = m4F;
    LGCOF = m4COF;
    LGCOF2 = m4COF2;
    LGZ = m4Z;
    LGCOFz = m4COFz;
    LGCOF2z = m4COF2z;
    
    try
        filename_save = [fullFilename(1:end-4) '_pooled_coher_', num2str(stax), '_',num2str(endax)];
        save([filename_save '.mat'],'TACOHT','TAF','TACOF','TACOF2','TAZ','TACOFz','TACOF2z',...
            'MGCOHT','MGF','MGCOF','MGCOF2','MGZ','MGCOFz','MGCOF2z',...
            'SolCOHT','SolF','SolCOF','SolCOF2','SolZ','SolCOFz','SolCOF2z',...
            'cstCOHT','cstF','cstCOF','cstCOF2','cstZ','cstCOFz','cstCOF2z',...
            'poolCOHT','poolF','poolCOF','poolCOF2','poolZ','poolCOFz','poolCOF2z',...
            'fsamp','MUNumber','MULength', 'WinLen','rSpike1','a','meanDR','CoVISI');
        print (gcf,'-dpdf',[filename_save '.pdf']);
        print(figure(99),'-dpdf',[filename_save '_comp.pdf']);
        
    catch
        try
            filename_save = [fullFilename(1:end-4) '_pooled_coher_', num2str(stax), '_',num2str(endax)];
            save([filename_save '.mat'],'TACOHT','TAF','TACOF','TACOF2','TAZ','TACOFz','TACOF2z',...
                'MGCOHT','MGF','MGCOF','MGCOF2','MGZ','MGCOFz','MGCOF2z',...
                'cstCOHT','cstF','cstCOF','cstCOF2','cstZ','cstCOFz','cstCOF2z',...
                'poolCOHT','poolF','poolCOF','poolCOF2','poolZ','poolCOFz','poolCOF2z',...
                'fsamp','MUNumber','MULength', 'WinLen','rSpike1','a','meanDR','CoVISI');
            print (gcf,'-dpdf',[filename_save '.pdf']);
            print(figure(99),'-dpdf',[filename_save '_comp.pdf']);
            
        catch
            try
                filename_save = [fullFilename(1:end-4) '_pooled_coher_', num2str(stax), '_',num2str(endax)];
                save([filename_save '.mat'],'MGCOHT','MGF','MGCOF','MGCOF2','MGZ','MGCOFz','MGCOF2z',...
                    'SolCOHT','SolF','SolCOF','SolCOF2','SolZ','SolCOFz','SolCOF2z',...
                    'cstCOHT','cstF','cstCOF','cstCOF2','cstZ','cstCOFz','cstCOF2z',...
                    'poolCOHT','poolF','poolCOF','poolCOF2','poolZ','poolCOFz','poolCOF2z',...
                    'fsamp','MUNumber','MULength', 'WinLen','rSpike1','a','meanDR','CoVISI');
                print (gcf,'-dpdf',[filename_save '.pdf']);
                print(figure(99),'-dpdf',[filename_save '_comp.pdf']);
                
            catch
                filename_save = [fullFilename(1:end-4) '_pooled_coher_', num2str(stax), '_',num2str(endax)];
                save([filename_save '.mat'],'MGCOHT','MGF','MGCOF','MGCOF2','MGZ','MGCOFz','MGCOF2z',...
                    'cstCOHT','cstF','cstCOF','cstCOF2','cstZ','cstCOFz','cstCOF2z',...
                    'poolCOHT','poolF','poolCOF','poolCOF2','poolZ','poolCOFz','poolCOF2z',...
                    'fsamp','MUNumber','MULength', 'WinLen','rSpike1','a','meanDR','CoVISI');
                print (gcf,'-dpdf',[filename_save '.pdf']);
                print(figure(99),'-dpdf',[filename_save '_comp.pdf']);
            end
        end
    end
    
elseif contains(fullFilename, 'Sol')
    
    TACOHT = m2COHT;
    TAF = m2F;
    TACOF = m2COF;
    TACOF2 = m2COF2;
    TAZ = m2Z;
    TACOFz = m2COFz;
    TACOF2z = m2COF2z;
    
    MGCOHT = m3COHT;
    MGF = m3F;
    MGCOF = m3COF;
    MGCOF2 = m3COF2;
    MGZ = m3Z;
    MGCOFz = m3COFz;
    MGCOF2z = m3COF2z;
    
    SolCOHT = m1COHT;
    SolF = m1F;
    SolCOF = m1COF;
    SolCOF2 = m1COF2;
    SolGZ = m1Z;
    SolCOFz = m1COFz;
    SolCOF2z = m1COF2z;
    
    LGCOHT = m4COHT;
    LGF = m4F;
    LGCOF = m4COF;
    LGCOF2 = m4COF2;
    LGZ = m4Z;
    LGCOFz = m4COFz;
    LGCOF2z = m4COF2z;
    
    try
        filename_save = [fullFilename(1:end-4) '_pooled_coher_', num2str(stax), '_',num2str(endax)];
        save([filename_save '.mat'],'TACOHT','TAF','TACOF','TACOF2','TAZ','TACOFz','TACOF2z',...
            'MGCOHT','MGF','MGCOF','MGCOF2','MGZ','MGCOFz','MGCOF2z',...
            'SolCOHT','SolF','SolCOF','SolCOF2','SolZ','SolCOFz','SolCOF2z',...
            'cstCOHT','cstF','cstCOF','cstCOF2','cstZ','cstCOFz','cstCOF2z',...
            'poolCOHT','poolF','poolCOF','poolCOF2','poolZ','poolCOFz','poolCOF2z',...
            'fsamp','MUNumber','MULength', 'WinLen','rSpike1','a','meanDR','CoVISI');
        print (gcf,'-dpdf',[filename_save '.pdf']);
        print(figure(99),'-dpdf',[filename_save '_comp.pdf']);
        
    catch
        try
            filename_save = [fullFilename(1:end-4) '_pooled_coher_', num2str(stax), '_',num2str(endax)];
            save([filename_save '.mat'],'TACOHT','TAF','TACOF','TACOF2','TAZ','TACOFz','TACOF2z',...
                'MGCOHT','MGF','MGCOF','MGCOF2','MGZ','MGCOFz','MGCOF2z',...
                'cstCOHT','cstF','cstCOF','cstCOF2','cstZ','cstCOFz','cstCOF2z',...
                'poolCOHT','poolF','poolCOF','poolCOF2','poolZ','poolCOFz','poolCOF2z',...
                'fsamp','MUNumber','MULength', 'WinLen','rSpike1','a','meanDR','CoVISI');
            print (gcf,'-dpdf',[filename_save '.pdf']);
            print(figure(99),'-dpdf',[filename_save '_comp.pdf']);
            
        catch
            try
                filename_save = [fullFilename(1:end-4) '_pooled_coher_', num2str(stax), '_',num2str(endax)];
                save([filename_save '.mat'],'MGCOHT','MGF','MGCOF','MGCOF2','MGZ','MGCOFz','MGCOF2z',...
                    'SolCOHT','SolF','SolCOF','SolCOF2','SolZ','SolCOFz','SolCOF2z',...
                    'cstCOHT','cstF','cstCOF','cstCOF2','cstZ','cstCOFz','cstCOF2z',...
                    'poolCOHT','poolF','poolCOF','poolCOF2','poolZ','poolCOFz','poolCOF2z',...
                    'fsamp','MUNumber','MULength', 'WinLen','rSpike1','a','meanDR','CoVISI');
                print (gcf,'-dpdf',[filename_save '.pdf']);
                print(figure(99),'-dpdf',[filename_save '_comp.pdf']);
                
            catch
                filename_save = [fullFilename(1:end-4) '_pooled_coher_', num2str(stax), '_',num2str(endax)];
                save([filename_save '.mat'], 'SolCOHT','SolF','SolCOF','SolCOF2','SolZ','SolCOFz','SolCOF2z',...
                    'cstCOHT','cstF','cstCOF','cstCOF2','cstZ','cstCOFz','cstCOF2z',...
                    'poolCOHT','poolF','poolCOF','poolCOF2','poolZ','poolCOFz','poolCOF2z',...
                    'fsamp','MUNumber','MULength', 'WinLen','rSpike1','a','meanDR','CoVISI');
                print (gcf,'-dpdf',[filename_save '.pdf']);
                print(figure(99),'-dpdf',[filename_save '_comp.pdf']);
            end
        end
    end
    
end

end

