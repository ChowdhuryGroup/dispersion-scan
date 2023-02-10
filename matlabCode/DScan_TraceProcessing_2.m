clear all
close all

% Process D-scan measurement into a common D-scan trace format to input
% into retrieval algorithm codes. 
% 
% Version 2 is designed to take d-scan traces output by the automated
% LabVIEW code. The expected input is a single matrix with the
% following format:
% 
% 0  l1  l2  l3  l4 ... (one zero followed by the wavelength array)
% 0  f1  f2  f3  f4 ... (one zero followed by the fundamental spectrum)
% p1 s11 s12 s13 s14 ... (initial stage pos; first row of raw d-scan)
% p2 s21 s22 s23 s24 ... (second stage pos; second row of d-scan)
% .  .
% .    .
% .      .
%
% NOTE THAT THIS VERSION DOES NOT CALIBRATE INPUT SPECTRA
%
% Each collected spectrum needs to be 
% 1) background-subtracted
% 2) plotted against frequency (instead of wavelength) and resampled for
% uniform sampling in frequency space
% 3) truncated to include only the SH, centering the frequency on that
% which corresponds to the central frequency of the fundamental spectrum
%
% The final output will be a .txt file with a matrix for the body,
% representing the D-scan trace (mm of dispersive material/compressor
% grating distance vs. SH frequency).  The first row of the matrix will
% correspond to the fundamental spectrum.
% The header will specify the central SH frequency (f0SH), frequency spacing (df), 
% starting position of dispersive material/stage position (Pi), spacing of 
% positions (dP):
%
% output file format:
% header:      f0SH, df, Pi, dP,
% fund spec:    s1, s2, s3, ... ,   sN
% Dscan trace: S11, S12, S13, ... , S1N
%              S21, S22, S23, ... , S2N
%               .    .    .     .    .
%               .    .    .     .    .
%              SM1, SM2, SM3, ... , SMN
%
% NxM matrix with N spectral samples and M stage position samples.
%
% The output file is TAB-DELIMITED

path = 'C:\Users\Noah\Documents\OhioStateU\ChowdhuryResearchGroup\data\DScan\GratingCompressor_Dscan\2020_07_27\';
fn='2020_07_27_64pos_3.txt';

% specify what to append to the end of the filename for the processed data
% file (e.g. '_Processed.txt')
procFN='_Processed_NoShift.txt';

dat=dlmread([path,fn],'\t');

% extract stage position, wavelength, fundamental spectrum, and dscan
% arrays
pos=dat(3:end,1);
lambda=dat(1,2:end);
fundspec=dat(2,2:end);
dscanRaw=dat(3:end,2:end);

% Compute corresponding frequency array:
rawFreqs=300./lambda; % PHz

% calculate average frequency spacing to be used as the sampling rate for
% the uniformly resampled frequency axis
df=mean(rawFreqs(1:end-1)-rawFreqs(2:end));
uniFreqs=min(rawFreqs):df:max(rawFreqs);

% background subtract by averaging over empty parts of the spectra:
BGfund=mean([fundspec(1000:1200),fundspec(1550:1600)]);
fundspec=fundspec-BGfund;

BGtrace=mean(dscanRaw(:,600:1000),2);
dscanRaw=dscanRaw-BGtrace;

% Resample over uniform frequency grid:
RStrace=zeros(length(pos),length(lambda));
for i=1:length(pos)
    F=griddedInterpolant(flip(rawFreqs),flip(dscanRaw(i,:)));
    RStrace(i,:)=300./(uniFreqs.^2).*F(uniFreqs);
end
F=griddedInterpolant(flip(rawFreqs),flip(fundspec));
RSfundSpec=300./(uniFreqs.^2).*F(uniFreqs);

% Find weighted average fundamental frequency:
freq0fund=sum(uniFreqs(60:300).*RSfundSpec(60:300))/sum(RSfundSpec(60:300));
[~, indf0]=min(abs(uniFreqs-freq0fund));

% Use twice the central frequency of the fundamental spectrum as the
% central frequency of the SH spectrum
[~, ind] = min(abs(uniFreqs-(freq0fund*2)));
f0SH=uniFreqs(ind);

% choose number of samples to truncate about the central SH frequency:
Nsamps=2^7;
finalFreqArray=uniFreqs(ind-Nsamps/2+1:ind+Nsamps/2)-f0SH;
finalTrace=RStrace(:,ind-Nsamps/2+1:ind+Nsamps/2);
finalFundSpec=RSfundSpec(indf0-Nsamps/2+1:indf0+Nsamps/2);

% Plot d-scan trace:
figure('Position',[250,200,500,400])
imagesc(finalFreqArray,pos,finalTrace);
xlabel('SH frequency (PHz)'); ylabel('Compressor Stage Position (mm)'); zlabel('SH signal');
colormap('jet');
shading interp
set(gca,'YDir','normal')

% Plot Fundamental Spectum:
figure('Position',[780,200,500,400])
plot(finalFreqArray,finalFundSpec)


% Save data to output file (SEE TOP OF CODE FOR FILE FORMAT DETAILS)
initialStagePos=pos(1);
ScanIncrement=pos(2)-pos(1);

str=input('Save data to disc? y/n: ','s');
if str == 'n'
    disp('rerun program with new settings');
elseif str == 'y'
    header=[f0SH, df, initialStagePos, ScanIncrement];
    fid=fopen([path,fn(1:end-4),procFN],'w');
    fprintf(fid, '%f\t',header);
    tline=fgets(fid);
    fclose(fid);
    
    body=cat(1,finalFundSpec,finalTrace);
    dlmwrite([path,fn(1:end-4),procFN],body,'-append','delimiter','\t','roffset',1);
    
else
    disp('incorrect input.  Rerun the program.');
end


