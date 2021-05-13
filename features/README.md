# Extracted motion features (0.2 Hz) from multisegmental, triaxial accelerometry for NIMS

## Data access instructions
The data corresponding to this directory can be found at the following citation link:
> Bhattacharyay, Shubhayu, 2021, "Extracted motion features (0.2 Hz) from multisegmental, triaxial accelerometry for NIMS", https://doi.org/10.7910/DVN/FMGE6E, Harvard Dataverse, V2, UNF:6:fwU1h/kxnsd45tYZZ5qHng== [fileUNF]

You will be asked to submit a request to access the data.

## Description
This dataset contains the motion features, extracted from non-overlapping non-overlapping 5-second windows (0.2 Hz) of filtered accelerometry, for each unique patient identifier (UPI). The format of the files are: `'features_' + UPI + '.csv'`.

The contents of each file include the UPI, the Recoding Index (which identifies unique timepoints in the recording, starting from 1), the hours from ICU admission of each time point, the time of day of the corresponding timepoint, the feature type (see below), and the columns containing corresponding feature values for each of the 7 sensors: Bed (control sensor), left ankle (LA), left elbow (LE), left wrist (LW), right ankle (RA), right elbow (RE), and right wrist (RW).

The feature types, encoded in the 'Feature' column, are as follows: `SMA` - signal magnitude area, `HLF_h` and `HLF_l` - a pair of high-frequency component (`HLF_h`) and low-frequency component (`HLF_l`) time-domain medians, `MFR` - median frequency, `FDE` - frequency-domain entropy, `BPW` - band power between 0.3 and 3.5 Hz, and and `WVL` - level 2 – 6 detail coefficients of the 5th-order Daubechies wavelet transform.

## Data format
### For patient-specific `.csv`:
Extracted features data is stored in an N x 12 table with the following column variables:
| UPI             | RecordingIdx             | HoursFromICUAdmission             | TimeOfDay             | Feature            | Bed             | LA             | LE             | LW             | RA             | RE             | RW             |
|-----------------|----------------|----------------|----------------|----------------|----------------|----------------|----------------|----------------|----------------|----------------|----------------|
