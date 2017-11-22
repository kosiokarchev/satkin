CREATE TABLE IF NOT EXISTS
cube (
    galaxyId INT PRIMARY KEY,
    haloId INT, subHaloId INT, fofCentralId INT, fofSubhaloId INT,
	centralMVir REAL, centralRvir REAL,
	distanceToCentralGalX REAL, distanceToCentralGalY REAL, distanceToCentralGalZ REAL,
	type INT,
    x REAL, y REAL, z REAL,
    velX REAL, velY REAL, velZ REAL,
    stellarSpinX REAL, stellarSpinY REAL, stellarSpinZ REAL,
    mvir REAL, rvir REAL, stellarMass REAL, sfr REAL,
    redshift REAL
)
WITHOUT ROWID