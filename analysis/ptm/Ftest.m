
% clear all;
% close all;

% F-testing

% vs. full model ----------------------

% enter the parameters of the model to be tested against the full model 
kreduced = 5; % number of parameters in reduced model
Rsquare_reduced = 0.977604;

% enter the full model parameters 
kfull = 6; % number of parameters in full model 
Rsquare_full = 0.98516;

N = 8*4; % number of observations

df1 = kfull - kreduced;
df2 = N - (kfull);


F = ( (Rsquare_full-Rsquare_reduced)/df1 ) / ( (1-Rsquare_full)/df2 );
Ftable = finv(0.95,df1,df2);
p = 1-fcdf(F,df1,df2);

if F >= Ftable
    fprintf(['The model is significantly different from the full model \n' ...
        'F(', num2str(df1), ',', num2str(df2), ') = ', num2str(F), ', p = ', num2str(p), '\n \n']); 
else
    fprintf(['The models are not significantly different \n' ...
        'F(', num2str(df1), ',', num2str(df2), ') = ', num2str(F), ', p = ', num2str(p), '\n \n']); 
end


% vs. default model --------------------

kfull2 = kreduced;
Rsquare_full2 = Rsquare_reduced;

% enter the parameters of the default model 
kreduced2 = 4;
Rsquare_reduced2 = 0.9511;

df1_2 = kfull2 - kreduced2;
df2_2 = N - kfull2;

F2 = ( (Rsquare_full2-Rsquare_reduced2)/df1_2 ) / ( (1-Rsquare_full2)/df2_2 );
Ftable2 = finv(0.95,df1_2,df2_2);
p2 = 1-fcdf(F2,df1_2,df2_2);

if F2 >= Ftable2
    fprintf(['The model is significantly different from the default model \n' ...
        'F(', num2str(df1_2), ',', num2str(df2_2), ') = ', num2str(F2), ', p = ', num2str(p2), '\n \n']); 
else
    fprintf(['The models are not significantly different \n' ...
        'F(', num2str(df1_2), ',', num2str(df2_2), ') = ', num2str(F2), ', p = ', num2str(p2), '\n \n']); 
end

