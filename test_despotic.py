from despotic import cloud

cld = cloud()
cld.nH      = 1.0e3
cld.colDen  = 1.0e22
cld.sigmaNT = 4.0e5
cld.Tg      = 10.0
cld.xoH2    = 0.1
cld.xpH2    = 0.4
cld.xHe     = 0.0
cld.xHI     = 0.000001
cld.xHplus  = 0.000001
cld.xe      = 0.000001
cld.addEmitter('CO', 1.0e-4)

lines = cld.lineLum('CO')