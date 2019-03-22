# validationPCR

You run the validation Script by:
python runValidation.py config.ini

This will store a object with a new something like validationPool_"run1Name"_"run2Name".obj

This can then be plotted using:
python plotValidation.py config

The results will be stored in a PDF files

For individual river plotting:
python plotTimeserie.py config.ini GRDC#1 GRDC#2 GRDC#x

Please not this is only on option when fullAnalysis is enabled in the config.ini
