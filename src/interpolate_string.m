function str = interpolate_string(str)
%INTERPOLATE_STRING Apply string interpolation
%   STR = INTERPOLATE_STRING(TEMPLATE) replaces all expressions surrounded
%   by curly braces in TEMPLATE with their value, and returns STR.
%   TEMPLATE must be a scalar string (or char array).
%
%   Expressions inside the curly braces cannot contain curly braces
%   themselves, for example, you can't use "c{1}" as an expression.
%
%   You can include curly braces in the output string by escaping the
%   opening brace with a backslash.
%
%   Examples:
%
%   >> interpolate_string('cos(4) = {cos(4)}')
%   ans =
%       'cos(4) = -0.6536'
%
%   >> s = 3.4; p = 'foo'
%   >> interpolate_string("The value of s is {s}, and that of p is {p}!")
%   ans =
%       "The value of s is 3.4000, and that of p is foo!"
%
%   >> interpolate_string("The value of s is {s}, \{this is not replaced}.")
%   ans = 
%       "The value of s is 3.4000, {this is not replaced}."
 
% (c)2022 by Cris Luengo
 
% Validate and normalize input string
return_char = false;
if ischar(str)
    str = string(str);
    return_char = true;
end
if ~isstring(str) || ~isscalar(str)
    error("Input must be a scalar string")
end
 
% Find all expressions to be evaluated
tokens = regexp(str, "(?<!\\)\{(.+?)\}", 'tokens');
tokens = unique([tokens{:}]);
 
% Evaluate the expressions
values = strings(size(tokens));
for ii = 1:numel(tokens)
    v = evalin('caller', tokens(ii));
    values(ii) = mlreportgen.utils.toString(v);
end
 
% Substitute values into string
for ii = 1:numel(tokens)
    t = regexptranslate('escape', tokens(ii));
    str = regexprep(str, "(?<!\\)\{(??@t)\}", values(ii));
end
 
% Remove escaped curly braces
str = regexprep(str, '\\\{', '{');
 
% Output should have same type as input
if return_char
    str = char(str);
end