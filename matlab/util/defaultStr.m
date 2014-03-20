function s = defaultStr(arg1, arg2, arg3)
% function defaultStr(varname,val)
% function default(varname, fieldname, defvalue)
%
% in two argument form, makes varname default to var in callers space
% in three argument form, 
%   for a structure varname, make defvalue the default for
%   field fieldname
%
% In this version, the default value is passed as a string
% that is only evaluated if necessary.
%


if (nargin == 2),
    varname = arg1;
    val = arg2;
    
    str = ['exist(''', varname, ''', ''var'')'];

    if ~(evalin('caller',str))
        assignin('caller',varname,evalin('caller',val));
    else
        str = ['isempty(', varname ')'];
        if (evalin('caller',str))
            assignin('caller',varname,evalin('caller',val));
        end
    end

end

if (nargin == 3),
    
    varname = arg1;
    fieldname = arg2;
    val = arg3;

    str = ['exist(''', varname, ''', ''var'')'];
    if ~(evalin('caller',str))
        assignin('caller',varname,setfield([],fieldname,evalin('caller',val)));
    else
        str = ['isfield(',varname,',''', fieldname,''')'];
        if ~(evalin('caller',str))
            prev = evalin('caller',varname);
            assignin('caller',varname,setfield(prev,fieldname,evalin('caller',val)));;
        end
    end

end
