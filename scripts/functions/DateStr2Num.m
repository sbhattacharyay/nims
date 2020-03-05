function D = DateStr2Num(S, F)
% DATESTR2NUM - Fast conversion of DATESTR to DATENUM
% The builtin DATENUM command is very powerful, but if you know the input
% format exactly, a specific MEX can be much faster: e.g. for single strings
% DateStr2Num is 50 to 100 times faster than DATENUM, for a {1 x 10000} cell
% string the speed up factor is 300 to 700 (Matlab 2009a/2011b/64, MSVC 2008).
%
% D = DateStr2Num(S, F)
% INPUT:
%   S: String or cell string in the format specified by F.
%      In opposite to DATENUM the validity of the input string is not checked
%      in any way (e.g. 1 <= month <= 12).
%   F: Integer number defining the input format. Accepted:
%           0: 'dd-mmm-yyyy HH:MM:SS'      01-Mar-2000 15:45:17
%           1: 'dd-mmm-yyyy'               01-Mar-2000
%          29: 'yyyy-mm-dd'                2000-03-01
%          30: 'yyyymmddTHHMMSS'           20000301T154517
%          31: 'yyyy-mm-dd HH:MM:SS'       2000-03-01 15:45:17
%      Not DATEFROM numbers of DATESTR:
%         230: 'mm/dd/yyyyHH:MM:SS'        12/24/201515:45:17
%         231: 'mm/dd/yyyy HH:MM:SS'       12/24/2015 15:45:17
%         240: 'dd/mm/yyyyHH:MM:SS'        24/12/201515:45:17
%         241: 'dd/mm/yyyy HH:MM:SS'       24/12/2015 15:45:17
%      Including the milliseconds:
%        1000: 'dd-mmm-yyyy HH:MM:SS.FFF'  01-Mar-2000 15:45:17.123
%        1030: 'yyyymmddTHHMMSS.FFF'       20000301T154517.123
%      The '-', ':', '/' and space characters are ignored and can be chosen
%      freely. Additional trailing characters are allowed.
%      Optional, default: 0.
%
% OUTPUT:
%   D: Serial date number. If S is a cell string, D has is same size.
%
% EXAMPLES:
%   C = {'2010-06-29 21:59:13', '2010-06-29 21:59:13'};
%   D = DateStr2Num(C, 31)
%   >> [734318.916122685, 734318.916122685]
% Equivalent Matlab command (but a column vector is replied):
%   D = datenum(C, 'yyyy-mm-dd HH:MM:SS')
%
% NOTES: The parsing of the strings works for clean ASCII characters only:
%   '0' MUST have the key code 48!
%   Month names must be English with the 2nd and 3rd charatcer in lower case.
%   Sorry for ugly numerical codes to define the format! This helps to squeeze
%   out some speed. If you need a 'nice' function, use DATENUM.
%   Trailing characters are ignore, such that you can ignore the time of
%   'dd-mmm-yyyy HH:MM:SS' by choosing F=1.
%
% COMPILATION:
%   This function must be compiled before using.
%   Windows: mex -O DateStr2Num.c
%   Linux:   mex -O CFLAGS="\$CFLAGS -std=c99" DateStr2Num.c
%   Download precompiled Mex: http://www.n-simon.de/mex
%
% Tested: Matlab/64 7.8, 7.13, 8.6, 9.1, Win7/64
% Compiler: LCC 2.4/3.8, BCC 5.5, Open Watcom 1.8, MSVC 2008
%         Compatibility to MacOS, Linux, 64 bit is assumed, but not tested.
% Author: Jan Simon, Heidelberg, (C) 2010-2018 matlab.2010(a)n(MINUS)simon.de
%
% See also DATESTR, DATENUM, DATEVEC.
% FEX: DateConvert 25594 (Jan Simon)

% $JRev: R-k V:011 Sum:OZg+ja45HeHv Date:27-May-2018 21:45:10 $
% $License: BSD (use/copy/change/redistribute on own risk, mention the author) $
% $UnitTest: uTest_DateStr2Num $
% $File: Tools\GLString\DateStr2Num.m $
% History:
% 001: 29-Jun-2010 22:07, First version.
% 002: 30-Jun-2010 16:08, Accept formats 0, 1, 29, 30, 31.
% 008: 28-Jul-2013 12:48, 300 -> 1030, format 1000 added.
% 009: 02-Dec-2014 08:08, New formats: 230,213,240,241, installation routine.
% 011: 27-May-2018 21:40, Upper/lower case in month names accepted in C-Mex.

% ==============================================================================
% The MEX is not found:
persistent WarnOnce     % Uncomment this to warn once only
if isempty(WarnOnce)
   warning('JSimon:DateStr2Num:NoMex', ...
      ['The slow M-version of DateStr2Num is used.\n', ...
      'Start the automatic compilation:']);
   InstallMex('DateStr2Num.c', 'uTest_DateStr2Num');
   WarnOnce = true;
end

% Compared to DATENUM, this M-version is 20 times faster than DATENUM(S, 31) for
% a string and 1.1 to 1.4 times slower for a {1 x 10000} cell string.
if ischar(S)
   S = {S};
elseif ~iscellstr(S)
   error('JSimon:DateStr2Num:BadInputType', ...
      'Input must be a string or cell string.');
end

if nargin == 1
   F = 0;
end

nS = numel(S);
D  = zeros(size(S));
switch F
   case 0  % 'dd-mmm-yyyy HH:MM:SS'
      for iS = 1:nS
         aS    = S{iS};
         x     = aS - '0';
         month = (strfind('janfebmaraprmayjunjulaugsepoctnovdec', ...
            lower(aS(4:6))) + 2) / 3;
         year  = x(8) * 1000 + x(9) * 100 + x(10) * 10 + x(11);
         
         D(iS) = datenummx(year, month, x(1) * 10 + x(2), ...
            x(13) * 10 + x(14), x(16) * 10 + x(17), ...
            x(19) * 10 + x(20));
      end
      
   case 1  % 'dd-mmm-yyyy'
      for iS = 1:nS
         aS    = S{iS};
         x     = aS - '0';
         month = (strfind('janfebmaraprmayjunjulaugsepoctnovdec', ...
            lower(aS(4:6))) + 2) / 3;
         year  = x(8) * 1000 + x(9) * 100 + x(10) * 10 + x(11);
         
         D(iS) = datenummx(year, month, x(1) * 10 + x(2));
      end
      
   case 29  % 'yyyy-mm-dd'
      for iS = 1:nS
         x     = S{iS} - '0';
         D(iS) = datenummx( ...
            x(1)*1000 + x(2)*100 + x(3)*10 + x(4), ... % Year
            x(6)*10 + x(7), ...  % Month
            x(9)*10 + x(10));     % Day
      end
      
   case 30  % 'yyyymmddTHHMMSS'
      for iS = 1:nS
         x     = S{iS} - '0';
         D(iS) = datenummx( ...
            x(1)*1000 + x(2)*100 + x(3)*10 + x(4), ... % Year
            x(5)*10 + x(6), ...  % Month
            x(7)*10 + x(8), ...  % Day
            x(10)*10 + x(11), x(12)*10 + x(13), x(14)*10 + x(15));
      end
      
   case 31  % 'yyyy-mm-dd HH:MM:SS'
      for iS = 1:nS
         x     = S{iS} - '0';
         D(iS) = datenummx( ...
            x(1)*1000 + x(2)*100 + x(3)*10 + x(4), ... % Year
            x(6)*10 + x(7), ...  % Month
            x(9)*10 + x(10), ... % Day
            x(12)*10 + x(13), x(15)*10 + x(16), x(18)*10 + x(19));
      end
      
   case 230  % 'mm/dd/yyyyHH:MM:SS'
      for iS = 1:nS
         x     = S{iS} - '0';
         D(iS) = datenummx( ...
            x(7)*1000 + x(8)*100 + x(9)*10 + x(10), ... % Year
            x(1)*10 + x(2), ...  % Month
            x(4)*10 + x(5), ...  % Day
            x(11)*10 + x(12), x(14)*10 + x(15), x(17)*10 + x(18));
      end
      
   case 231  % 'mm/dd/yyyy HH:MM:SS'
      for iS = 1:nS
         x     = S{iS} - '0';
         D(iS) = datenummx( ...
            x(7)*1000 + x(8)*100 + x(9)*10 + x(10), ... % Year
            x(1)*10 + x(2), ...  % Month
            x(4)*10 + x(5), ...  % Day
            x(12)*10 + x(13), x(15)*10 + x(16), x(18)*10 + x(19));
      end
      
   case 240  % 'dd/mm/yyyyHH:MM:SS'
      for iS = 1:nS
         x     = S{iS} - '0';
         D(iS) = datenummx( ...
            x(7)*1000 + x(8)*100 + x(9)*10 + x(10), ... % Year
            x(4)*10 + x(5), ...  % Month
            x(1)*10 + x(2), ...  % Day
            x(11)*10 + x(12), x(14)*10 + x(15), x(17)*10 + x(18));
      end
      
   case 241  % 'dd/mm/yyyy HH:MM:SS'
      for iS = 1:nS
         x     = S{iS} - '0';
         D(iS) = datenummx( ...
            x(7)*1000 + x(8)*100 + x(9)*10 + x(10), ... % Year
            x(4)*10 + x(5), ...  % Month
            x(1)*10 + x(2), ...  % Day
            x(12)*10 + x(13), x(15)*10 + x(16), x(18)*10 + x(19));
      end
      
   case 1000  % 'dd-mmm-yyyy HH:MM:SS.FFF'
      for iS = 1:nS
         aS    = S{iS};
         x     = aS - '0';
         month = (strfind('janfebmaraprmayjunjulaugsepoctnovdec', ...
            lower(aS(4:6))) + 2) / 3;
         year  = x(8) * 1000 + x(9) * 100 + x(10) * 10 + x(11);
         
         D(iS) = datenummx(year, month, x(1) * 10 + x(2), ...
            x(13) * 10 + x(14), x(16) * 10 + x(17), ...
            x(19) * 10 + x(20) + ...
            x(22) / 10 + x(23) / 100 + x(24) / 1000);
      end
      
   case 1030  % yyyymmddTHHMMSS.FFF
      for iS = 1:nS
         x     = S{iS} - '0';
         D(iS) = datenummx( ...
            x(1)*1000 + x(2)*100 + x(3)*10 + x(4), ... % Year
            x(5)*10 + x(6), ...  % Month
            x(7)*10 + x(8), ...  % Day
            x(10)*10 + x(11), x(12)*10 + x(13), x(14)*10 + x(15) + ...
            x(17) / 10 + x(18) / 100 + x(19) / 1000);
      end
      
   otherwise
      error('JSimon:DateStr2Num:BadFormat', 'Bad format');
end

% return;
