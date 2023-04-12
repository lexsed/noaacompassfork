# noaacompassfork

This is a fork of NOAA magnetic compass Web App. The original app is available at https://www.ngdc.noaa.gov/geomag/calculators/mobileCompass.shtml

Per https://www.ngdc.noaa.gov/ngdcinfo/privacy.html#disclaimer, the original autor of the app is the National Geophysical Data Center (NGDC) of the National Oceanic and Atmospheric Administration (NOAA). Being published by NOAA, the app is presumed to be in the public domain. The app does depend on external libraries and APIs, which may be subject to different licenses, however references to such dependencies have been removed in this fork.

The modifications made to the original app are:
- removal of depedencies and trackers
- removal of the magnetic storm indicator (and corresponding API requests)
- removal of the help button and corresponding dependencies
- addition of location indicator (only updates when the page is reloaded)