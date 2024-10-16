#!/bin/bash
# Submit read manifest files via ENA Webin Client
for filename in manifest-ELA-L*.txt
do
	java -jar webin-cli-6.7.2.jar -ascp -context reads -userName Webin-50792 -password jTR8EEI9rR -manifest $filename -submit
done
