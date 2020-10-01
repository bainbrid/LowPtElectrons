import pandas as pd
import sys

ntuples = True

f = None
if len(sys.argv) > 1 : f = open(str(sys.argv[1]),'r')
else : f = open('pippo','r')
if f is None : quit()

name = None
if len(sys.argv) > 2 : name = str(sys.argv[2])
else : name = 'inAODSIM_'

issues = []
paths = []
jobids = []
sizes = []
for line in f.readlines() :
    try :
        print line.split(name)
        jobid = int(line.split(name)[-1].split('.root')[0])
        jobids.append(jobid)
        if len(line.split()) == 9 : sizes.append(int(line.split()[4]))
    except ValueError :
        #if 'MINIAOD' not in line : 
        issues.append(line.strip('\n'))

print 
print 'Number of jobids: ',len(jobids)
print 'Number of issues: ',len(issues)
print '# unique jobids:  ',len(list(set(jobids)))
if len(jobids) != len(list(set(jobids))) : 
    print "Duplicates!", list(set([ x for x in jobids if jobids.count(x) > 1 ]))
print 'Minimum jobid:    ',min(jobids)
print 'Maximum jobid:    ',max(jobids)
print 'Minimum size:     ',min(sizes)
print 'Maximum size:     ',max(sizes)
missing = []
for jobid in range(min(jobids),max(jobids)+1) :
    if jobid not in jobids : missing.append(jobid)
print 'Missing jobids:   ',[ x for x in missing ]
print

print 'df.describe():'
df = pd.DataFrame.from_dict({'jobid':jobids,'size':sizes})
print df.describe(include='all')
print

print 'Issues parsing following text:'
for index,issue in enumerate(issues) : print "  #{:.0f}: {:s}".format(index,issue)
