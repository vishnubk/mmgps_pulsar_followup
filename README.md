# mmgps_pulsar_followup
**A Nextflow pipeline for followup of Pulsar Discoveries in the MPIfR MeerKAT Galatic Plane Survey (MMGPS)**

INSTALLATION

1. Install JAVA and Nextflow. Follow the instructions here https://www.nextflow.io/docs/latest/getstarted.html
2. Create a .env file within this repo. Fill in the MMGPS Database credentials (username, password, host etc).
3. Fill in the nextflow.config file with your target_name (eg J0737-3039), beam_name, spin period range, utc start and utc end to cover. If no UTC is specified it will analyse all observations.
4. If you do not have an Ephemeris for your pulsar yet, I recommend only folding the APSUSE filterbank files first i.e APSUSE_FOLDS = 1 and PTUSE_FOLDS = 0. This will run filtool, peasoup and fold all candidates within your spin period range using dspsr automatically.
5. Once you create TOAs, and get a preliminary ephemeries for the APSUSE observations, set PTUSE_FOLDS = 1, and fold all your high resolution PTUSE files.

