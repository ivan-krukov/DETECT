"""Check the dependencies for DETECT"""

import subprocess

def run_which (query):
	p=subprocess.Popen(("which",query),stdout=subprocess.PIPE)
	stdin,stderr=p.communicate()
	return stdin

if __name__=="__main__":
	fail=False
	if not run_which("blastp"):
		print "It appears that BLASTP is not installed or is not in the path.\nPlease consult http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download"
		fail=True
	if not run_which("blastdbcmd"):
		print "It appears that BLASTDBCMD is not installed or is not in the path.\nPlease consult http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download  "
		fail=True
	if not run_which("needle"):
		print "It appears that NEEDLE is not installed or is not in the path.\nPlease consult http://emboss.sourceforge.net/download/"
		fail=True
	if not fail:
		print "All dependencies appear to be present"
	else:
		exit(1)
