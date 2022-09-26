function [str,ndis] = normdis(dis)

% dis=varargin(1);
bdis=dis-dis(1);

ndis=bdis/dis(1);
str=max(ndis)-min(ndis);

