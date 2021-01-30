function [] = plot_clinical_cmp(db_rep,db0,db1,dateFrom,dateTo,varargin)

ip = inputParser;
addParamValue(ip, 'mm', true, @islogical);%#ok<*NVREPL>
addParamValue(ip, 'reported',true, @islogical);
addParamValue(ip, 'reduced',true, @islogical);
parse(ip, varargin{:});
results = ip.Results;
mm = results.mm;
reported = results.reported;
reduced = results.reduced;

figure('Name','Situation in Hospitals I.');
if reduced
    subplot(2,1,1)
    if mm
        bar(resize(db0.H,dateFrom:dateTo),'FaceAlpha',0.5,'EdgeAlpha',0.5);hold on;
        bar(resize(db1.H,dateFrom:dateTo),'FaceAlpha',0.5,'EdgeAlpha',0.5);hold on;
    end
    pp1=plot(resize(db0.H_smooth,dateFrom:dateTo),'linewidth',2,'Color','b');hold on;
    pp2=plot(resize(db1.H_smooth,dateFrom:dateTo),'linewidth',2,'Color','r');hold on;
    if reported
        pp3=plot(resize(db_rep.H_smooth,dateFrom:dateTo),'linewidth',1,'Color','k','linestyle','--');hold on;
        legend([pp1,pp2,pp3],{'observed','implied by reported daily new cases','reconstructed from implied cases'});
    else
        legend([pp1,pp2],{'observed','implied by reported daily new cases'});
    end
    grid on;
    title('Hospitalisations (total)');    
      
    subplot(2,1,2)
    if mm
        bar(resize(db0.D,dateFrom:dateTo),'FaceAlpha',0.5,'EdgeAlpha',0.5);hold on;
        bar(resize(db1.D,dateFrom:dateTo),'FaceAlpha',0.5,'EdgeAlpha',0.5);hold on;
    end
    pp1=plot(resize(db0.D_smooth,dateFrom:dateTo),'linewidth',2,'Color','b');hold on;
    pp2=plot(resize(db1.D_smooth,dateFrom:dateTo),'linewidth',2,'Color','r');hold on;
    if reported
        pp3=plot(resize(db_rep.D_smooth,dateFrom:dateTo),'linewidth',1,'Color','k','linestyle','--');hold on;
        legend([pp1,pp2,pp3],{'observed','implied by reported daily new cases','reconstructed from implied cases'});
    else
        legend([pp1,pp2],{'observed','implied by reported daily new cases'});
    end
    grid on;
    title('Cummulative deaths (total)');    
else
    subplot(2,2,1);
    if mm
        bar(resize(db0.H,dateFrom:dateTo),'FaceAlpha',0.5,'EdgeAlpha',0.5);hold on;
        bar(resize(db1.H,dateFrom:dateTo),'FaceAlpha',0.5,'EdgeAlpha',0.5);hold on;
    end
    pp1=plot(resize(db0.H_smooth,dateFrom:dateTo),'linewidth',2,'Color','b');hold on;
    pp2=plot(resize(db1.H_smooth,dateFrom:dateTo),'linewidth',2,'Color','r');hold on;
    if reported
        pp3=plot(resize(db_rep.H_smooth,dateFrom:dateTo),'linewidth',1,'Color','k','linestyle','--');hold on;
        legend([pp1,pp2,pp3],{'observed','implied by reported daily new cases','reconstructed from implied cases'});
    else
        legend([pp1,pp2],{'observed','implied by reported daily new cases'});
    end
    grid on;
    title('Hospitalisations (total)');    

    subplot(2,2,2);
    if mm
        bar(resize(db0.C,dateFrom:dateTo),'FaceAlpha',0.5,'EdgeAlpha',0.5);hold on;
        bar(resize(db1.C,dateFrom:dateTo),'FaceAlpha',0.5,'EdgeAlpha',0.5);hold on;
    end
    pp1=plot(resize(db0.C_smooth,dateFrom:dateTo),'linewidth',2,'Color','b');hold on;
    pp2=plot(resize(db1.C_smooth,dateFrom:dateTo),'linewidth',2,'Color','r');hold on;
    if reported
        pp3=plot(resize(db_rep.C_smooth,dateFrom:dateTo),'linewidth',1,'Color','k','linestyle','--');hold on;
        legend([pp1,pp2,pp3],{'observed','implied by reported daily new cases','reconstructed from implied cases'});
    else
        legend([pp1,pp2],{'observed','implied by reported daily new cases'});
    end
    grid on;
    title('ICU hospitalisations');    

    subplot(2,2,3)
    if mm
        bar(resize(db0.V,dateFrom:dateTo),'FaceAlpha',0.5,'EdgeAlpha',0.5);hold on;
        bar(resize(db1.V,dateFrom:dateTo),'FaceAlpha',0.5,'EdgeAlpha',0.5);hold on;
    end
    pp1=plot(resize(db0.V_smooth,dateFrom:dateTo),'linewidth',2,'Color','b');hold on;
    pp2=plot(resize(db1.V_smooth,dateFrom:dateTo),'linewidth',2,'Color','r');hold on;
    if reported
        pp3=plot(resize(db_rep.V_smooth,dateFrom:dateTo),'linewidth',1,'Color','k','linestyle','--');hold on;
        legend([pp1,pp2,pp3],{'observed','implied by reported daily new cases','reconstructed from implied cases'});
    else
        legend([pp1,pp2],{'observed','implied by reported daily new cases'});
    end
    grid on;
    title('Ventilations');    
    
    subplot(2,2,4)
    if mm
        bar(resize(db0.D,dateFrom:dateTo),'FaceAlpha',0.5,'EdgeAlpha',0.5);hold on;
        bar(resize(db1.D,dateFrom:dateTo),'FaceAlpha',0.5,'EdgeAlpha',0.5);hold on;
    end
    pp1=plot(resize(db0.D_smooth,dateFrom:dateTo),'linewidth',2,'Color','b');hold on;
    pp2=plot(resize(db1.D_smooth,dateFrom:dateTo),'linewidth',2,'Color','r');hold on;
    if reported
        pp3=plot(resize(db_rep.D_smooth,dateFrom:dateTo),'linewidth',1,'Color','k','linestyle','--');hold on;
        legend([pp1,pp2,pp3],{'observed','implied by reported daily new cases','reconstructed from implied cases'});
    else
        legend([pp1,pp2],{'observed','implied by reported daily new cases'});
    end
    grid on;
    title('Cummulative deaths (total)');    
end

% % %
% % figure('Name','Situation in Hospitals: Comparison I')
% % if idx_fun==1
% %     subplot(2,2,1)
% %     plot(resize(hospit_smooth,disp_from:t1),'linewidth',1);hold on;
% %     plot(resize(out.H,disp_from:t1),'linewidth',1);hold on;
% %     legend({'observed','implied by reported daily new cases'});
% %     grid on;
% %     title('Hospitalisations (total)');
% %     subplot(2,2,2)
% %     plot(resize(icu_smooth,disp_from:t1),'linewidth',1);hold on;
% %     plot(resize(out.C,disp_from:t1),'linewidth',1);hold on;
% %     grid on;
% %     title('ICU');
% %     subplot(2,2,3)
% %     plot(resize(vent_smooth,disp_from:t1),'linewidth',1);hold on;
% %     plot(resize(out.V,disp_from:t1),'linewidth',1);hold on;
% %     grid on;
% %     title('Ventilations');
% %     subplot(2,2,4)
% %     plot(resize(deaths_total_smooth,disp_from:t1),'linewidth',1);hold on;
% %     plot(resize(out.D,disp_from:t1),'linewidth',1);hold on;
% %     grid on;
% %     title('Deaths');
% % else
% %     subplot(2,1,1)
% %     ratio_h = resize(hospit_smooth,disp_from:t1)./resize(out_check.H,disp_from:t1);
% %     plot(resize(hospit_smooth,disp_from:t1),'linewidth',2);hold on;
% %     plot(ratio_h.*resize(out.H,disp_from:t1),'linewidth',2);hold on;
% %     % plot(ratio.*resize(out_check.H,disp_from:t1),'k--','linewidth',1);hold on;
% %     legend({'observed','implied by reported daily new cases'});%,'reconstructed from implied cases'});
% %     grid on;
% %     title('Hospitalisations (total)');    
% %     xls_out.H_imp = ratio_h.*resize(out.H,disp_from:t1);
% %     xls_out.H_rep = resize(hospit_smooth,disp_from:t1);
% %       
% %     subplot(2,1,2)
% %     ratio_d = resize(deaths_total_smooth,disp_from:t1)./resize(out_check.D,disp_from:t1);
% %     plot(resize(deaths_total_smooth,disp_from:t1),'linewidth',2);hold on;
% %     plot(ratio_d.*resize(out.D,disp_from:t1),'linewidth',2);hold on;
% %     % plot(resize(out_check.D,disp_from:t1),'k--','linewidth',1);hold on;
% %     grid on;
% %     title('Deaths');
% %     xls_out.D_imp = ratio_d.*resize(out.D,disp_from:t1);
% %     xls_out.D_rep = resize(deaths_total_smooth,disp_from:t1);
% %   
% %     legend({'observed','implied by reported daily new cases'}); %,'reconstructed from implied cases'}); 
% %     figure('Name','Situation in Hospitals: Comparison II');
% %     subplot(2,1,1)
% %     ratio_c = 1/3*ratio_d+2/3*ratio_h;
% %     plot(resize(icu_smooth,disp_from:t1),'linewidth',2);hold on;
% %     plot(resize(icu_smooth,disp_from:t1)./(2/3*resize(hospit_smooth,disp_from:t1)./(ratio_h.*resize(out.H,disp_from:t1))+1/3*resize(deaths_total_smooth,disp_from:t1)./(ratio_d.*resize(out.D,disp_from:t1))),'linewidth',2);hold on;
% %     legend({'observed','implied by reported daily new cases'});%,'reconstructed from implied cases'});
% %     grid on;
% %     title('ICU');
% %     subplot(2,1,2)
% %     plot(resize(vent_smooth,disp_from:t1),'linewidth',2);hold on;
% %     plot(resize(vent_smooth,disp_from:t1)./(1/3*resize(hospit_smooth,disp_from:t1)./(ratio_h.*resize(out.H,disp_from:t1))+2/3*resize(deaths_total_smooth,disp_from:t1)./(ratio_d.*resize(out.D,disp_from:t1))),'linewidth',2);hold on;
% %     grid on;
% %     title('Ventilations');
% %     legend({'observed','implied by reported daily new cases'}); %,'reconstructed from implied cases'}); 
% % end    
% % 
% % 
% % 
% % figure('Name','Clinical statistics');
% % %
% % subplot(2,2,1)
% % if mm
% %     bar(resize(data.H,dateFrom:dateTo));hold on;
% % end
% % if smooth
% %     plot(resize(data.H_smooth,dateFrom:dateTo),'linewidth',2);hold on;
% % end
% % if raw
% %     plot(resize(data.H_raw,dateFrom:dateTo),'linestyle','-.','Color',[0.5 0.5 0.5]);hold on;
% % end
% % grid on;
% % title('Total hospitalizations');
% % %
% % subplot(2,2,2)
% % if mm
% %     bar(resize(data.C,dateFrom:dateTo));hold on;
% % end
% % if smooth
% %     plot(resize(data.C_smooth,dateFrom:dateTo),'linewidth',2);hold on;
% % end
% % if raw
% %     plot(resize(data.C_raw,dateFrom:dateTo),'linestyle','-.','Color',[0.5 0.5 0.5]);hold on;
% % end
% % grid on;
% % title('ICU hospitalizations');
% % %
% % subplot(2,2,3)
% % if mm
% %     bar(resize(data.V,dateFrom:dateTo));hold on;
% % end
% % if smooth
% %     plot(resize(data.V_smooth,dateFrom:dateTo),'linewidth',2);hold on;
% % end
% % if raw
% %     plot(resize(data.V_raw,dateFrom:dateTo),'linestyle','-.','Color',[0.5 0.5 0.5]);hold on;
% % end
% % grid on;
% % title('Ventilations');
% % %
% % subplot(2,2,4)
% % if mm
% %     bar(resize(data.D,dateFrom:dateTo));hold on;
% % end
% % if smooth
% %     plot(resize(data.D_smooth,dateFrom:dateTo),'linewidth',2);hold on;
% % end
% % if raw
% %     plot(resize(data.D_raw,dateFrom:dateTo),'linestyle','-.','Color',[0.5 0.5 0.5]);hold on;
% % end
% % grid on;
% % title('Cummulative deaths (on+with)');

end