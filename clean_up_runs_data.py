#!/usr/bin/env python

###

import subprocess
import os
import numpy as np

###

for x in range (1, 11):
	runs_data_in = open(str(x) + "clade_support/clade_support.txt", "r+")
	runs_data_out = open(str(x) + "clade_support/clade_support_out.txt", "w+")
	runs_data_in_base = runs_data_in.read()
	numbers = range(0, 2000)
	runs_data_in_base = runs_data_in_base.replace("[1]", "clade_support<-c(")
	for q in range(0, 2000):
		runs_data_in_base = runs_data_in_base.replace("[" + str(numbers[q]) + "]", "  ")
		runs_data_in_base = runs_data_in_base.replace("[" + str(numbers[q]) + "]", "  ")
	runs_data_in_base = runs_data_in_base.replace("      ", " ")
	runs_data_in_base = runs_data_in_base.replace("     ", " ")
	runs_data_in_base = runs_data_in_base.replace("   ", " ")
	runs_data_in_base = runs_data_in_base.replace("  ", " ")
	runs_data_in_base = runs_data_in_base.replace(" ", ", ")
	runs_data_in_base = runs_data_in_base + ")"
	runs_data_in_base = runs_data_in_base.replace(", clade", "clade")
	runs_data_in_base = runs_data_in_base.replace("(, ", "(")
	runs_data_out.write(runs_data_in_base)
