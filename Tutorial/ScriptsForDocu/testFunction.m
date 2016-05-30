%% Test Function
% The test function as follows:
%
% $$f_{\rm{Test}} = \frac{10x_{1}}{5 + x_{1}}\sin\left( x_{2} \right) +
% x_{3}$$
%
% That is, the output is related to the first input variable by a Michaelis
% Menten curve(steep increase for small values followed by a plateau). The
% second input variable leads to oscillation and the third input variable
% leads to an linear monotonically increase. The effect of the
% first and second input are coupled with each other.
%
% For a better visualization, pair-wise combination of input variables are
% plotted against the output. Remaining input variable was hold constant at
% the value of 2.
visualizeTestFunction();

