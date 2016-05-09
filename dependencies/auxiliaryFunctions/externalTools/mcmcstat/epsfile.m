function epsfile(file,stretch,varargin)
%EPSFILE  dumps current image into eps file
% epsfile(filename)
% epsfile(filename,1) sets 'PaperPositionMode' to 'auto'

% $Revision: 1.3 $  $Date: 2013/07/11 09:32:19 $

if length(file)==0
  error 'usage epsfile file'
end
if nargin > 1 & stretch~=0
  set(gcf,'PaperPositionMode','auto');
end
print('-depsc2','-noui',varargin{:},file)
