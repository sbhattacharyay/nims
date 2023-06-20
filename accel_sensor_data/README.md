# Unfiltered triaxial accelerometry (10 Hz) for NIMS

## Data access instructions
The data corresponding to this directory can be found at the following citation link:
> Bhattacharyay, Shubhayu; Wang, Matthew; Rattray, John; Etienne-Cummings, Ralph; Kudela, Pawel; Stevens, Robert David, 2021, "Unfiltered triaxial accelerometry (10 Hz) for NIMS", https://doi.org/10.7910/DVN/SXNBGB, Harvard Dataverse, V1

You will be asked to submit a request to access the data.

## Description
Raw accelerometry (10 Hz, unit: g) file names are formatted as: `‘accel_’ + unique patient identifier (UPI) + ‘.mat’`. For patients for whom we collected data too large to fit into one file, data is split into one file per sensor (7 sensors total), and file names are formatted as follows: `‘accel_’ + unique patient identifier (UPI) + ‘_C’ + sensor number (1 – 7) + ‘.mat’`.

For complete `.mat` files, each file contains a MATLAB cell array (12 x 7) called `data`. The columns of this cell array correspond to the different sensors (in order: Bed [control sensor], Left Ankle, Left Elbow, Left Wrist, Right Ankle, Right Elbow, Right Wrist) and the rows of the cell array correspond to the different channels of the sensors. For this study, only rows 5-7 (x-,y-, and z-axes respectively) and row 12 (time of day stamps) are used, and thus rows 1-4 and 8-11 can be removed.

For sensor-specific `.mat` files (ending in `_C#.mat`), each file contains a MATLAB cell array (1 x 12) called `C#`. This corresponds to the different channels of the sensor. For this study, only cells 5-7 (x-,y-, and z-axes respectively) and cell 12 (time of day stamps) are used, and thus cells 1-4 and 8-11 can be removed. The C# correspond to sensors as follows: C1 - Bed [control sensor], C2 - Left Ankle, C3 - Left Elbow, C4 - Left Wrist, C5 - Right Ankle, C6 - Right Elbow, and C7 - Right Wrist. To concatenate these channels into a 'data' cell array, as seen in complete `.mat` files, load all the `_C#.mat` for a specific patient, and run the following line of MATLAB code: `data = [C1;C2;C3;C4;C5;C6;C7]';`

## Data format
### For complete `.mat` files:
`data` is a 12 x 7 cell array that can be labeled as follows:
| Bed             | LA             | LE             | LW             | RA             | RE             | RW             |
|-----------------|----------------|----------------|----------------|----------------|----------------|----------------|
| Blank           | Blank          | Blank          | Blank          | Blank          | Blank          | Blank          |
| Blank           | Blank          | Blank          | Blank          | Blank          | Blank          | Blank          |
| Blank           | Blank          | Blank          | Blank          | Blank          | Blank          | Blank          |
| Blank           | Blank          | Blank          | Blank          | Blank          | Blank          | Blank          |
| x-axis (Bed)    | x-axis (LA)    | x-axis (LE)    | x-axis (LW)    | x-axis (RA)    | x-axis (RE)    | x-axis (RW)    |
| y-axis (Bed)    | y-axis (LA)    | y-axis (LE)    | y-axis (LW)    | y-axis (RA)    | y-axis (RE)    | y-axis (RW)    |
| z-axis (Bed)    | z-axis (LA)    | z-axis (LE)    | z-axis (LW)    | z-axis (RA)    | z-axis (RE)    | z-axis (RW)    |
| Blank           | Blank          | Blank          | Blank          | Blank          | Blank          | Blank          |
| Blank           | Blank          | Blank          | Blank          | Blank          | Blank          | Blank          |
| Blank           | Blank          | Blank          | Blank          | Blank          | Blank          | Blank          |
| Blank           | Blank          | Blank          | Blank          | Blank          | Blank          | Blank          |
| TimeOfDay (Bed) | TimeOfDay (LA) | TimeOfDay (LE) | TimeOfDay (LW) | TimeOfDay (RA) | TimeOfDay (RE) | TimeOfDay (RW) |

### For sensor-specific `.mat` files (ending in `_C#.mat`):
`C#` is a 1 x 12 cell array that can be labeled as follows:
| Blank             | Blank             | Blank             | Blank             | x-axis            | y-axis             | z-axis             | Blank             | Blank             | Blank             | Blank             | TimeOfDay             |
|-----------------|----------------|----------------|----------------|----------------|----------------|----------------|----------------|----------------|----------------|----------------|----------------|
