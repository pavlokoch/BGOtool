function counts = trigger_lookups(bg, window) 

% Look-up tables from ASIM-TER-MXGS-ICD-001_2E MXGS Software Interface Control Document.pdf

% bg - bakground count rate
% window 1 - short 1
% window 2 - short 2
% window 3 - short 3
% window 4 - long

% check range
if bg<0.1
%     warning('Background count rate below minimal 0.1 counts/ms. Value 0.1 is used instead.');
    bg = 0.1;
end
if bg>12.8
%     warning('Background count rate above maximal 12.8 counts/ms. Value 12.8 is used instead.');
    bg = 12.8;
end

% round bg to 0.1 precision 
bg = round(bg,1);

% energies
en = round(0.1:0.1:12.8,1);

% HED short trigger window 1 variable thresholds look-up table
% (HED_STW_VAR_THR_1)
% A 1-dimensional array containing thresholds to be used for
% HED short windows calculated background rates in range
% 0.1 – 12.8 counts/ms (in increments of 0.1 counts/ms)
% Values range: 0 – 255 Units: Counts

HED_STW_VAR_THR_1 = [4, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 8, 8, 9, 9, 9, 9, ...
9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, ...
12, 12, 12, 12, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 13, ...
13, 13, 13, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 15, 15, 15, ...
15, 15, 15, 15, 15, 15, 15, 15, 16, 16, 16, 16, 16, 16, 16, 16, ...
16, 16, 16, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, ...
18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 19, 19, 19, 19, ...
19, 19, 19, 19, 19, 19, 19, 19, 19, 20, 20, 20, 20, 20];

% HED short trigger window 2 variable thresholds look-up table
% (HED_STW_VAR_THR_2)
% A 1-dimensional array containing thresholds to be used for
% HED short windows calculated background rates in range
% 0.1 – 12.8 counts/ms (in increments of 0.1 counts/ms)
% Values range: 0 – 255 Units: Counts

HED_STW_VAR_THR_2 = [6, 7, 7, 8, 9, 9, 10, 10, 10, 11, 11, 12, 12, ...
12, 13, 13, 13, 14, 14, 14, 14, 15, 15, 15, 16, 16, 16, 16, 17, ...
17, 17, 17, 18, 18, 18, 18, 19, 19, 19, 19, 20, 20, 20, 20, 21, ...
21, 21, 21, 22, 22, 22, 22, 22, 23, 23, 23, 23, 24, 24, 24, 24, ...
24, 25, 25, 25, 25, 25, 26, 26, 26, 26, 26, 27, 27, 27, 27, 27, ...
28, 28, 28, 28, 28, 29, 29, 29, 29, 29, 30, 30, 30, 30, 30, 31, ...
31, 31, 31, 31, 31, 32, 32, 32, 32, 32, 33, 33, 33, 33, 33, 33, ...
34, 34, 34, 34, 34, 35, 35, 35, 35, 35, 35, 36, 36, 36, 36, 36, ...
37, 37, 37];

% HED short trigger window 3 variable thresholds look-up table
% (HED_STW_VAR_THR_3)
% A 1-dimensional array containing thresholds to be used for
% HED short windows calculated background rates in range
% 0.1 – 12.8 counts/ms (in increments of 0.1 counts/ms)
% Values range: 0 – 255 Units: Counts

HED_STW_VAR_THR_3 = [7, 9, 10, 11, 12, 13, 14, 15, 16, 16, 17, 18, ...
19, 19, 20, 21, 21, 22, 23, 23, 24, 24, 25, 26, 26, 27, 27, 28, ...
28, 29, 30, 30, 31, 31, 32, 32, 33, 33, 34, 34, 35, 35, 36, 37, ...
37, 38, 38, 39, 39, 40, 40, 41, 41, 42, 42, 43, 43, 44, 44, 44, ...
45, 45, 46, 46, 47, 47, 48, 48, 49, 49, 50, 50, 51, 51, 52, 52, ...
52, 53, 53, 54, 54, 55, 55, 56, 56, 57, 57, 57, 58, 58, 59, 59, ...
60, 60, 61, 61, 61, 62, 62, 63, 63, 64, 64, 65, 65, 65, 66, 66, ...
67, 67, 68, 68, 68, 69, 69, 70, 70, 71, 71, 71, 72, 72, 73, 73, ...
74, 74, 74, 75];

% HED long trigger window variable thresholds look-up table
% (HED_LTW_VAR_THR)
% A 1-dimensional array containing thresholds to be used for
% LED long window calculated background rates in range 0.1 –
% 12.8 counts/ms (in increments of 0.1 counts/ms)
% Values range: 0 – 1023 Units: Counts

HED_LTW_VAR_THR = [14, 20, 24, 29, 33, 38, 42, 45, 49, 53, 57, ...
61, 64, 67, 71, 74, 78, 81, 85, 88, 92, 95, 98, 102, 105, 109, ...
112, 115, 119, 122, 125, 128, 132, 135, 138, 141, 145, 148, ...
151, 154, 157, 161, 164, 167, 170, 173, 176, 179, 183, 186, ...
189, 192, 195, 198, 201, 204, 207, 210, 213, 216, 219, 223, ...
226, 229, 232, 235, 238, 241, 244, 247, 250, 253, 256, 259, ...
262, 265, 268, 271, 274, 277, 280, 283, 286, 289, 292, 295, ...
297, 300, 303, 306, 309, 312, 315, 318, 321, 324, 327, 330, ...
333, 336, 339, 342, 345, 347, 350, 353, 356, 359, 362, 365, ...
368, 371, 374, 377, 379, 382, 385, 388, 391, 394, 397, 400, ...
403, 405, 408, 411, 414, 417];

switch window
    case 1
        counts = HED_STW_VAR_THR_1(en==bg);
    case 2
        counts = HED_STW_VAR_THR_2(en==bg);
    case 3
        counts = HED_STW_VAR_THR_3(en==bg);
    case 4
        counts = HED_LTW_VAR_THR(en==bg);
    otherwise
        error('Window must be 1,2,3 for short windows or 4 for long.')
end
        
          






