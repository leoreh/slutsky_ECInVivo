function sig = accurateResampling(original_sigal,f1,f2)
%
% Resampling function based on interpolating data with an spline-like
% function and reconstructing signal at the new sampling frequency.
%
% NOTE: The function works only on vectors with length grater than 80,000!
% 
% Inputs : 
% original_sigal: the original signal.
% f1: original sampling frequency
% f2: new sampling frequency
%
% Outputs :
% sig: The signal after resampleing
%

if length(original_sigal) < 80000
    error('Resampling Error: For this function, input signal need to contain at least 80,000 values');
end

sig = original_sigal;
sigIsInColumn = true;
if size(sig,1) >= size(sig,2)
  sig = sig';
else
    sigIsInColumn = false;
end


tail_length=20000;
max_length=round(length(sig)-3*tail_length);
steps=floor(length(sig)/max_length);
if mod(length(sig),max_length)>2*tail_length
    steps=steps+1;
end
index=1;
sig_=[];
for i=1:steps-1
    sig_{i}=sig(index:index+max_length);
    index=index+max_length-tail_length;
end
sig_{steps}=sig(index:length(sig));
sig_end=sig(end-tail_length:end);
for i=1:3
    sig_end=[sig_end,sig_end];
end
for i=1:length(sig_)
    sig_app=sig_{i};
    t=1/f1:1/f1:length(sig_app)/f1;
    pp=pchip(t,sig_app);
    t2=1/f2:1/f2:max(t);
    clear t;
    sig_{i}=ppval(pp,t2);
    clear t2;
end
t=1/f1:1/f1:length(sig_end)/f1;
t2=1/f2:1/f2:max(t);
pp=pchip(t,sig_end);
sig_end=ppval(pp,t2);
index=1;
sig=[];
for i=1:steps-1
    app=sig_{i};
    sig(index:index+round(max_length*f2/f1)-1)=app(1:round(max_length*f2/f1));
    index=index+round((max_length-tail_length)*f2/f1);
end
sig(index:index+length(sig_{steps})-1)=sig_{steps};
sig(end-round(tail_length*f2/f1):end)=[];
sig(end+1:end+round(tail_length*f2/f1))=sig_end(1:round(tail_length*f2/f1));

if sigIsInColumn, sig = sig'; end;