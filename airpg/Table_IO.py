import os, logging
import pandas as pd

class TableIO:

	def __init__(self, fp_entry_table, fp_ir_table = None, fp_blacklist = None, fp_duplicates = None, logger = None):
		self.log = logger or logging.getLogger(__name__ + ".TableIO")
		self.entry_table = None
		self.duplicates = {}
		self.ir_table = None
		self.blacklist = []

		self.read_entry_table(os.path.abspath(fp_entry_table))

		if fp_ir_table:
			self.read_ir_table(os.path.abspath(fp_ir_table))
		if fp_blacklist:
			self.read_blacklist(os.path.abspath(fp_blacklist))
		if fp_duplicates:
			self.read_duplicates(os.path.abspath(fp_duplicates))

	###############
	# I/O methods #
	###############

	def read_entry_table(self, fp_entry_table):
		'''
		Read a tab-separated file of GenBank entry information.
		If the file doesn't exist yet, create it and write column headers
		Params:
		 - fp_entry_table: file path to input file
		'''
		if os.path.isfile(fp_entry_table):
			self.entry_table = pd.read_csv(fp_entry_table, sep = '\t', index_col = 0, encoding = 'utf-8')
		else:
			columns = ["UID", "ACCESSION", "VERSION", "ORGANISM", "SEQ_LEN", "CREATE_DATE", "AUTHORS", "TITLE", "REFERENCE", "NOTE", "TAXONOMY"]
			self.entry_table = pd.DataFrame(columns = columns)
			self.entry_table = self.entry_table.set_index("UID", drop = True)
			self.write_entry_table(fp_entry_table)

	def write_entry_table(self, fp_entry_table, append = False):
		'''
		Write a list of GenBank entry information to tab-separated file.
		Params:
		 - fp_entry_table: file path to output file
		'''
		if append:
			self.entry_table.to_csv(fp_entry_table, sep = '\t', encoding = 'utf-8', header = False, mode = "a")
		else:
			self.entry_table.to_csv(fp_entry_table, sep = '\t', encoding = 'utf-8', header = True)

	def append_entry_to_table(self, entry, uid, fp_entry_table):
		'''
		Write information on one GenBank entry to tab-separated file
		Params:
		 - entry: dict. Keys are column names
		 - uid: Unique identifier for this GenBank entry
		 - fp_entry_table: file path to output file
		'''
		if os.path.isfile(fp_entry_table):
			for key, value in entry.items():
				entry[key] = [value]
			entry["UID"] = [uid]
			temp_df = pd.DataFrame(entry)
			temp_df = temp_df.set_index("ACCESSION", drop = True)
			temp_df.to_csv(fp_entry_table, sep = '\t', header = False, encoding = 'utf-8', mode = "a")
		else:
			raise Exception("Error trying to append GenBank entry to file '%s': File does not exist!" % (fp_entry_table))

	def read_ir_table(self, fp_ir_table):
		'''
		Read a tab-separated file of information on inverted repeats per GenBank accession
		If the file doesn't exist yet, create it and write column headers
		Params:
		 - fp_ir_table: file path to input file
		'''
		if os.path.isfile(fp_ir_table):
			self.ir_table = pd.read_csv(fp_ir_table, sep = '\t', index_col = 0, encoding = 'utf-8')
		else:
			columns = ["ACCESSION", "IRa_REPORTED", "IRa_REPORTED_START", "IRa_REPORTED_END", "IRa_REPORTED_LENGTH", "IRb_REPORTED", "IRb_REPORTED_START", "IRb_REPORTED_END", "IRb_REPORTED_LENGTH"]
			self.ir_table = pd.DataFrame(columns = columns)
			self.ir_table = self.ir_table.set_index("ACCESSION", drop = True)
			self.write_ir_table(fp_ir_table)

	def write_ir_table(self, fp_ir_table, append = False):
		'''
		Write a list of per-accession inverted repeat information to tab-separated file
		Params:
		 - fp_ir_table: file path to output file
		'''
		if append:
			self.ir_table.to_csv(fp_ir_table, sep = '\t', encoding = 'utf-8', header = False, mode = "a")
		else:
			self.ir_table.to_csv(fp_ir_table, sep = '\t', encoding = 'utf-8', header = True)

	def append_ir_info_to_table(self, ir_info, accession, fp_ir_table):
		'''
		Write information on one accession's inverted repeats to tab-separated file
		Params:
		 - ir_info: dict. Keys are column names
		 - accession: accession number of this record
		 - fp_ir_table: file path to output file
		'''
		if os.path.isfile(fp_ir_table):
			for key, value in ir_info.items():
				ir_info[key] = [value]
			ir_info["ACCESSION"] = [accession]
			temp_df = pd.DataFrame(ir_info)
			temp_df = temp_df.set_index("ACCESSION", drop = True)
			temp_df.to_csv(fp_ir_table, sep = '\t', header = False, encoding = 'utf-8', mode = "a")
		else:
			raise Exception("Error trying to append IR info to file '%s': File does not exist!" % (fp_ir_table))

	def read_blacklist(self, fp_blacklist):
		'''
		Read a file of blacklisted genera.
		Params:
		 - fp_blacklist: file path to input file
		'''
		with open(fp_blacklist, "r") as fh_blacklist:
			for line in [l.strip() for l in fh_blacklist.readlines()]:
				if not line.startswith("#"):
					self.blacklist.append(line)

	def read_duplicates(self, fp_duplicates):
		'''
		Read a tab-separated file of UIDs, corresponding RefSeq accession numbers and duplicate accession numbers.
		Params:
		 - fp_duplicates: file path to input file
		'''
		with open(fp_duplicates, "r") as fh_duplicates:
			for dup_tup in [line.rstrip().split('\t') for line in fh_duplicates.readlines()]:
				self.duplicates[dup_tup[0]] = [dup_tup[1], dup_tup[2]]

	def write_duplicates(self, fp_duplicates):
		'''
		Write a list of UIDs, corresponding RefSeq accession numbers and duplicate accession numbers to tab-separated file.
		Params:
		 - fp_duplicates: file path to output file
		'''
		with open(fp_duplicates, "w") as fh_duplicates:
			for d_key in self.duplicates.keys():
				fh_duplicates.write("%s\t%s\t%s\n" % (str(d_key), str(self.duplicates[d_key][0]), str(self.duplicates[d_key][1])))

	def append_duplicates(self, fp_duplicates):
		'''
		Append a list of UIDs, corresponding RefSeq accession numbers and duplicate accession numbers to tab-separated file.
		Params:
		 - fp_duplicates: file path to output file
		'''
		with open(fp_duplicates, "a") as fh_duplicates:
			for d_key in self.duplicates.keys():
				fh_duplicates.write("%s\t%s\t%s\n" % (str(d_key), str(self.duplicates[d_key][0]), str(self.duplicates[d_key][1])))

	#####################
	# List edit methods #
	#####################

	def remove_naturally_irl_genera(self, genera_list):
		'''
		Remove entries from ir_table that belong to a genus that naturally lacks one or both IRs.
		'''
		genus_accessions = {}
		for genus in genera_list:
			# get accessions from entry_table where ORGANISM starts with genus
			# check if these accession have two reported IRs in ir_table
			# if not, remove it from ir_table
			genus_accessions[genus] = list(self.entry_table.loc[[(entry[-1].rstrip('.') == genus) for entry in self.entry_table["TAXONOMY"].str.split(';')]]["ACCESSION"])
		self.log.debug("Found %s accessions that belong to potential IR-lacking genera." % str(sum(len(genus_accessions[x]) for x in genus_accessions)))
		for genus, accessions in genus_accessions.items():
			lacks_irs = True
			i = 0
			while lacks_irs == True and i < len(accessions):
				if self.ir_table[accessions[i]]["IRa_REPORTED"] == "yes" and self.ir_table[accessions[i]]["IRb_REPORTED"] == "yes":
					lacks_irs = False
					self.log.debug("Found an entry with two reported IRs for " + str(genus))
				i += 1
			if lacks_irs:
				self.log.debug("Dropping accessions for genus " + str(genus))
				for accession in accessions:
					self.log.debug("Dropping accession " + str(accession))
					self.ir_table.drop(accession, inplace=True)

	def remove_blacklisted_entries(self):
		'''
		Remove entries from entry table that match blacklisted genera.
		'''
		for genus in self.blacklist:
			# TM: The next line took a while to figure out, so for the sake of my own and future contributers' sanities, here's a breakdown of what it does:
			# self.entry_table["TAXONOMY"].str provides the whole taxonomy column for elementwise(i.e. rowwise) string operations. Since our TAXONOMY information is semicolon-separated, each row is split.
			# This results in a Series of string lists. The last element (the genus) of each string list is compared to the current genus from the blacklist (entry[-1].rstrip('.') == genus).
			# This in turn results in a list of bools, making self.entry_table.loc return all rows where the list of bools has True (i.e. all entries that match a blacklisted entry)
			# Finally, we want only the index of those rows, to tell the dataframe which ones should get dropped.
			self.entry_table.drop(self.entry_table.loc[[(entry[-1].strip('. ') == genus) for entry in self.entry_table["TAXONOMY"].str.split(';')]].index, inplace = True)
		if len(self.blacklist) == 0:
			self.log.info("Blacklist is empty. No entries removed.")

	def remove_duplicates(self):
		'''
		Remove entries from entry table that match duplicate accession numbers
		'''
		for d_key in self.duplicates.keys():
			self.entry_table.drop(self.entry_table.loc[self.entry_table["ACCESSION"] == self.duplicates[d_key][1]].index, inplace = True)
