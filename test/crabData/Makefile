# This script is used to create crab config files for jobs.
# Different templates are used for different era definition.

# The different eras are processedn as different classes of targets.
# This pipeline is driven using the submit shell scripts as targets.

.PHONY: all clean 2022 2023 2024

# Define the different eras
all: submit_2022CD.sh submit_2022E.sh submit_2022FG.sh submit_2023.sh submit_2024.sh

# Phony targets for preparing eras
2022: submit_2022CD.sh submit_2022E.sh submit_2022FG.sh

2023: submit_2023.sh

2024: submit_2024.sh

# For each era, use corresponding CMSSW, crab config template and data list to create crab config file.
submit_%.sh: crab3_template_%.py RunDatalist_%.txt ../runMultiLepPAT_dataRun3_miniAOD_Run%.py registerData.sh
	@echo "Creating crab config file for era $*"
	@./registerData.sh -l RunDatalist_$*.txt -t crab3_template_$*.py -s $@
	@chmod +x $@
# Clean up the crab config files
clean:
	@rm -f submit_*.sh
	@rm -f crab3_*_Run*_MINIAOD.py
	@rm -f crab3_*_Run*_MINIAOD.py-e

	
