function [pdf_data_smooth,cdf_data_smooth] = apply_censoring(cdf_data,k,s_width,s_type,s_ends)

%
pdf_data_new = cdf_data(2:end)-cdf_data(1:end-1); 
pdf_data_new = [pdf_data_new;0];
pdf_data_new(k+1:end) = 0;
pdf_data_new = pdf_data_new./sum(pdf_data_new);

%
cdf_data = cumsum(pdf_data_new);
cdf_data_smooth = smooth_series(cdf_data,s_width,s_type,s_ends);
cdf_data_smooth(1) = 0;

%
pdf_data_smooth = cdf_data_smooth(2:end)-cdf_data_smooth(1:end-1);
pdf_data_smooth = smooth_series([pdf_data_smooth;0],ceil(s_width/2),s_type,ends);
pdf_data_smooth = pdf_data_smooth/sum(pdf_data_smooth);

end
