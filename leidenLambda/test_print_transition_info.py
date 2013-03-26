from leidenLambda import molData 

specStr  = 'O'
inPath   = 'oatom.dat'
idx      = 0
T_kin    = 50.0
collider = 'p-H2' 
#----------------------------------------------------------------------------------------------
#reading the whole database of line info of species from LAMBDA
lambdaPath = '/home/mher/ism/code/radex/Radex/data/home.strw.leidenuniv.nl/~moldata/datafiles'
reader = molData.reader(dirPath = lambdaPath, species = specStr)
specInfo = reader.get_specie(specStr, inPath)

####picking manually####
specInfo = specInfo[1]
######################

transRad = specInfo.transRad[idx]

u, l = transRad['u'], transRad['l']
nu = transRad['nu']
en = transRad['E']

print 'upper-lower = %d-%d' % (u, l) 
print 'frequency   = %.2f GHz' % nu 
print 'wavelength  = %.2f um' % (3e5/nu)
print 'eneryg      = %.2f K' % en
ncr = specInfo.critical_density(upper=u, lower=l, T_kin=T_kin, collider=collider)

print 'n_critical  = %.2e' % ncr
print 'done'
