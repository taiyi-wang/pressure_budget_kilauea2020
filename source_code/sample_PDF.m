function samples = sample_PDF(PDF, Nsamples, replacement_flag)
% Sample posterior distributions
% Input: 
% PDF              = a N X 1 vector of a PDF distribution or a matrix of N X Ndistr PDFs
% Nsamples         = number of samples to take from the distribution
% replacement_flag = whether to replace after sampling

% Output:
% samples          = a Nsamples X 1 vector of samples for one distribution
%                    or a Nsamples X Ndistributions matrix of samples

% Taiyi Wang 08/26/21

Ndistr = length(PDF(1,:));
samples = zeros(Nsamples, Ndistr);
for i = 1:Ndistr
    PDFi = PDF(:,i);
    samples(:, i) = randsample(PDFi, Nsamples, replacement_flag);
end

end